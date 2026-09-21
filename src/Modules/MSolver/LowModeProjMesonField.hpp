/*
 * LowModeProjMesonField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Michael Lynch <michaellynch628@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 */

/*  END LEGAL */
#ifndef HadronsMILC_MSolver_LowModeProjMesonField_hpp_
#define HadronsMILC_MSolver_LowModeProjMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <A2AMatrix.hpp>
#include <EigenPack.hpp>
#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *         Low mode projector from a precomputed meson field                *
 ******************************************************************************/
/*  Drop-in alternative to TLowModeProjMILC (projector branch) that takes
    the <eta_i|e_j>(t)-type contraction coefficients from a meson-field HDF5
    file (loaded by MIO::LoadMesonField) instead of computing them live.

    The file's eigenvector rows come in |e+o>/|e-o> pairs (rows 2k/2k+1 of
    eigenpair k) -- produced identically by EigenPackFullPairs (full-volume)
    or by MesonField's checkerboarded cbPairs path (on-demand CB pairs; the
    writer's CB 1/sqrt(2) norm reduction absorbs the pair construction's
    1/sqrt(2), so the rows are numerically identical either way). Their sum
    and difference reconstruct the parity-split inner products LowModeProj
    computes live:

      SUM_k(t,j)  = M[t][2k][j] + M[t][2k+1][j] ~ <e_E^k|eta_j,E>(t)
      DIFF_k(t,j) = M[t][2k][j] - M[t][2k+1][j] ~ (i/lam_k) <Meooe(e_E^k)|eta_j,O>(t)

    giving (see the design artifact, decision D8, for the full derivation):

      sol = (norm_cb / pairScale) * [ SUM_k e_k  -  i * Meooe( SUM_k (DIFF_k/lam_k) e_k ) ]

    with norm_cb = 1/sqrt(norm2(evec[0])) computed live on the checkerboarded
    pack, exactly as LowModeProj. Everything that cannot come from the file
    (the Meooe parity move, the 1/lam_D^2 weighting folded into DIFF/lam,
    the final normalization and the axpy onto live checkerboarded
    eigenvectors) stays live.

    action      Staggered action module (Meooe parity move)
    lowModes    MassShiftEigenPack module with the CHECKERBOARDED
                eigenvectors/eigenvalues (row pair k <-> evec[k], eval[k])
    gammas      "" (default) or a standard spin-taste pair list
                "(G1 G1) (G5 G5)" (StagGamma::ParseSpinTasteString, the
                same parser the meson-field writer uses): one solver
                family per gamma; empty means the legacy single (G1,G1)
                family with bare names
    mesonField  whitespace-separated LoadMesonField module names, one per
                gamma (positional parallel list, gammas[i] <-> entry i),
                each holding a [nt, 2*nEvec, nNoise] table whose
                coefficients already carry that gamma
    noiseIndex  j: which noise vector (file column) this instance projects
    eigStart    first eigenpair to include (pair space)
    nEigs       number of eigenpairs (< 1: all)
    negFirst    ""/"false" (default): row 2k is |e+o>; "true": |e-o> comes
                first (flips the DIFF sign)
    pairScale   production normalization constant P (default sqrt(2), i.e.
                files whose |e+o> vectors were built as
                (e_E +/- Meooe(e_E)/(i*lam))/sqrt(2) from a unit-checkerboard
                pack; use 1.0 for unit-norm full-volume Lanczos files)
    noise       optional name of the noise-vector environment object
                (std::vector<FermionField>, e.g. "<noise module>_vec") used
                for a one-shot pairing/normalization self-check at execute
                time (only the (G1,G1) table is checkable against the live,
                gamma-independent reference; other gammas log a skip)

    Instead of a fixed timeslice parameter, the module creates one solver
    pair per lattice timeslice t in [0, nt) PER GAMMA g: "<name>_t<t>" and
    "<name>_t<t>_subtract" for a single gamma (or an empty gammas list --
    bit-identical to the previous fixed-G1_G1 behavior), and
    "<name>_<spin>_<taste>_t<t>" / "..._subtract" when the list holds more
    than one gamma (GaugeProp/Meson single-gamma-no-suffix convention),
    each a closure over the same solver body with {table, timeslice} bound
    by value (the "_subtract" suffix stays last, following LowModeProj's
    name/name_subtract convention, so appending "_subtract" to any
    wrapper name yields its subtract variant). The full 2*nt*nGamma name
    family is enumerated in getOutput() from env().getDim().back() (the
    first HadronsMILC module with a lattice-dimension-dependent output
    list -- FermionFlow/ScalarVP core precedents). A solver call ignores
    the source passed to it (the coefficients come from the file); the
    "<name>..._subtract" variant returns source - projection, mirroring
    LowModeProj. The reconstruction body is gamma-blind: the file rows
    already carry Gamma (see the meson-field writer), so each family
    applies the identical SUM/DIFF algebra to its own table.
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class LowModeProjMesonFieldPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LowModeProjMesonFieldPar,
                                  std::string,   action,
                                  std::string,   lowModes,
                                  std::string,   mesonField,
                                  std::string,   gammas,
                                  unsigned int,  noiseIndex,
                                  unsigned int,  eigStart,
                                  int,           nEigs,
                                  std::string,   negFirst,
                                  std::string,   pairScale,
                                  std::string,   noise);
  LowModeProjMesonFieldPar(void)
      : gammas(""), negFirst(""), pairScale("") {}
};

template <typename FImpl, typename Pack>
class TLowModeProjMesonField : public Module<LowModeProjMesonFieldPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );
  SOLVER_TYPE_ALIASES(FImpl, );

private:
  // gammas "(G1 G1) (G5 G5)" parsed with the standard spin-taste pair
  // parser (empty -> the legacy single (G1,G1) default, bare names);
  // duplicate gamma names are fatal (solver names would collide)
  std::vector<StagGamma::SpinTastePair> gammaList(void) const;
  // one LoadMesonField module name per gamma (positional parallel list,
  // strToVec<std::string>); count mismatch vs the gamma list is fatal
  // before any positional access
  std::vector<std::string> mesonFieldList(
      const std::vector<StagGamma::SpinTastePair> &gammas) const;

public:
  // constructor
  TLowModeProjMesonField(const std::string name);
  // destructor
  virtual ~TLowModeProjMesonField(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  virtual DependencyMap getObjectDependencies(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(StagLMAMesonField,
                    ARG(TLowModeProjMesonField<STAGIMPL,
                                               MassShiftEigenPack<STAGIMPL>>),
                    MSolver);

/******************************************************************************
 *                  TLowModeProjMesonField implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TLowModeProjMesonField<FImpl, Pack>::TLowModeProjMesonField(
    const std::string name)
    : Module<LowModeProjMesonFieldPar>(name) {}

// gamma/meson-field parallel lists //////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<StagGamma::SpinTastePair>
TLowModeProjMesonField<FImpl, Pack>::gammaList(void) const {
  std::vector<StagGamma::SpinTastePair> gammas;
  if (par().gammas.empty()) {
    gammas.push_back(std::make_pair(StagGamma::StagAlgebra::G1,
                                    StagGamma::StagAlgebra::G1));
  } else {
    gammas = StagGamma::ParseSpinTasteString(par().gammas);
  }
  std::vector<std::string> names;
  for (auto &g : gammas) {
    std::string name = StagGamma::GetName(g);
    for (auto &n : names) {
      if (n == name) {
        HADRONS_ERROR(Argument, "duplicate gamma '" + name +
                                    "' in gammas list (solver names would "
                                    "collide)");
      }
    }
    names.push_back(name);
  }

  return gammas;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLowModeProjMesonField<FImpl, Pack>::mesonFieldList(
    const std::vector<StagGamma::SpinTastePair> &gammas) const {
  auto mfs = strToVec<std::string>(par().mesonField);
  if (mfs.size() != gammas.size()) {
    HADRONS_ERROR(Argument,
                  "gammas list has " + std::to_string(gammas.size()) +
                      " entries but mesonField names " +
                      std::to_string(mfs.size()) +
                      " module(s): one LoadMesonField module per gamma "
                      "(positional parallel list)");
  }

  return mfs;
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TLowModeProjMesonField<FImpl, Pack>::getInput(void) {
  auto mfs = mesonFieldList(gammaList());
  std::vector<std::string> in{par().action, par().lowModes};
  for (auto &mf : mfs) {
    in.push_back(mf);
  }
  if (!par().noise.empty()) {
    in.push_back(par().noise);
  }

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLowModeProjMesonField<FImpl, Pack>::getOutput(void) {
  std::vector<std::string> out;
  auto gammas = gammaList();
  int nt = env().getDim().back();

  // single gamma -> bare names (legacy, bit-identical); multiple gammas ->
  // "_<spin>_<taste>" segment before the timeslice (GaugeProp.hpp naming
  // convention), "_subtract" stays last
  for (unsigned int g = 0; g < gammas.size(); ++g) {
    std::string suffix =
        (gammas.size() > 1) ? ("_" + StagGamma::GetName(gammas[g])) : "";
    for (int t = 0; t < nt; ++t) {
      std::string name = getName() + suffix + "_t" + std::to_string(t);
      out.push_back(name);
      out.push_back(name + "_subtract");
    }
  }

  return out;
}

template <typename FImpl, typename Pack>
DependencyMap
TLowModeProjMesonField<FImpl, Pack>::getObjectDependencies(void) {
  DependencyMap dep;
  auto gammas = gammaList();
  auto mfs = mesonFieldList(gammas);
  int nt = env().getDim().back();

  // every wrapper name carries its own {action, lowModes, mesonField[g]}
  // (and noise) edges so the GC never frees a loaded table before
  // late-timeslice consumers run; each family references only its own
  // gamma's loader
  for (unsigned int g = 0; g < gammas.size(); ++g) {
    std::string suffix =
        (gammas.size() > 1) ? ("_" + StagGamma::GetName(gammas[g])) : "";
    for (int t = 0; t < nt; ++t) {
      for (auto &s : {"", "_subtract"}) {
        std::string name = getName() + suffix + "_t" + std::to_string(t) + s;
        dep.insert({par().action, name});
        dep.insert({par().lowModes, name});
        dep.insert({mfs[g], name});
        if (!par().noise.empty()) {
          dep.insert({par().noise, name});
        }
      }
    }
  }

  return dep;
}

/******************************************************************************
 *              TLowModeProjMesonField setup                                  *
 ******************************************************************************/

template <typename FImpl, typename Pack>
void TLowModeProjMesonField<FImpl, Pack>::setup(void) {
  int Ls = env().getObjectLs(par().action);
  if (Ls > 1) {
    HADRONS_ERROR(Argument, "Ls > 1 not implemented");
  }

  // optional string parameters (missing XML nodes are tolerated for strings
  // only -- useStencil precedent, commit 297b724)
  bool negFirst = false;
  if (par().negFirst == "true") {
    negFirst = true;
  } else if ((par().negFirst != "") && (par().negFirst != "false")) {
    HADRONS_ERROR(Argument, "negFirst must be '', 'false' or 'true' (got '" +
                                par().negFirst + "')");
  }
  RealD pairScale = std::sqrt(2.0);
  if (!par().pairScale.empty()) {
    auto scale = strToVec<RealD>(par().pairScale);
    if (scale.size() != 1) {
      HADRONS_ERROR(Argument, "pairScale must be a single number (got '" +
                                  par().pairScale + "')");
    }
    pairScale = scale[0];
  }

  auto gammas = gammaList();
  auto mfs = mesonFieldList(gammas);

  LOG(Message) << "Setting up meson-field driven low mode projector "
               << "for action '" << par().action
               << "' using eigenvectors from '" << par().lowModes
               << "' and meson fields (noise index " << par().noiseIndex
               << ", one solver pair per timeslice per gamma):"
               << std::endl;
  for (unsigned int g = 0; g < gammas.size(); ++g) {
    LOG(Message) << "  gamma '" << StagGamma::GetName(gammas[g])
                 << "' from '" << mfs[g] << "'" << std::endl;
  }

  auto &mat = envGet(FMat, par().action);
  auto &epack = envGet(Pack, par().lowModes);

  // eigenpair range: the fixed three-clause bounds check of LowModeProj
  // (commit a64d61b)
  unsigned int eigStart = par().eigStart;
  int nEigs = par().nEigs;
  if (nEigs < 1) {
    nEigs = epack.evec.size();
  }
  if (eigStart > nEigs || eigStart > epack.evec.size() ||
      nEigs - eigStart > epack.evec.size() - eigStart) {
    HADRONS_ERROR(Argument,
                  "Requested eigs (parameters eigStart and nEigs) out of "
                  "bounds.");
  }

  // meson-field tables: one per gamma; rows must be |e+o>/|e-o> pairs of
  // the pack eigenvectors (2 rows per eigenpair)
  int nt = env().getDim().back();
  std::vector<std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>> *> mfTables;
  for (auto &mfName : mfs) {
    auto &mf =
        envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, mfName);
    if (static_cast<int>(mf.size()) != nt) {
      HADRONS_ERROR(Size, "meson field '" + mfName + "' has " +
                              std::to_string(mf.size()) +
                              " timeslices, expected " +
                              std::to_string(nt));
    }
    mfTables.push_back(&mf);
  }
  // NOTE: the row-count (mf[0].rows() == 2*evec.size()), column-count
  // (noiseIndex < mf[0].cols()) and D9 self-check validations live in
  // execute(), not here: the scheduler dry-runs every module's setup()
  // for memory profiling BEFORE any module executes, when the loader's
  // container still holds nt empty 0x0 matrices

  envCache(FermionField, "rbTemp", 1, envGetRbGrid(FermionField));
  envCache(FermionField, "rbTempNeg", 1, envGetRbGrid(FermionField));
  envCache(FermionField, "rbFermNeg", 1, envGetRbGrid(FermionField));

  // bind the fixed instance index NOW (LowModeProj captures bound values
  // by value, not par() at call time); the table index, timeslice and
  // subtract flag enter as factory arguments so each wrapper closes over
  // its own copy
  const unsigned int jIdx = par().noiseIndex;

  // the solver body is gamma-blind: the file coefficients already carry
  // Gamma, so every family applies the identical SUM/DIFF algebra to its
  // own table
  auto makeSolver = [&mat, &epack, mfTables, eigStart, nEigs, negFirst,
                     pairScale, jIdx, this](const unsigned int g,
                                            const unsigned int tSlice,
                                            bool subGuess) {
    auto *mf = mfTables[g];
    return [&mat, &epack, mf, subGuess, eigStart, nEigs, negFirst, pairScale,
            tSlice, jIdx, this](FermionField &sol,
                                const FermionField &source) {
      auto &rbTemp = envGet(FermionField, "rbTemp");
      auto &rbTempNeg = envGet(FermionField, "rbTempNeg");
      auto &rbFermNeg = envGet(FermionField, "rbFermNeg");

      int cb = epack.evec[0].Checkerboard();

      RealD norm = 1. / ::sqrt(norm2(epack.evec[0]));

      const A2AMatrix<HADRONS_A2AM_IO_TYPE> &mft = (*mf)[tSlice];
      const unsigned int j = jIdx;

      rbTemp = Zero();
      rbTemp.Checkerboard() = cb;
      rbTempNeg = Zero();
      rbTempNeg.Checkerboard() = cb;
      rbFermNeg = Zero();
      rbFermNeg.Checkerboard() = (cb == Even) ? Odd : Even;

      // accumulate the two parity channels from the file row pairs; the
      // subtraction channel is accumulated as DIFF/lam so that the live
      // Meooe below completes the 1/lam_D^2 weighting of LowModeProj
      for (int k = (eigStart + nEigs - 1); k >= int(eigStart); k--) {
        const FermionField &e = epack.evec[k];
        const RealD lam_D = epack.eval[k].imag();

        ComplexD sum = ComplexD(mft(2 * k, j)) +
                       ComplexD(mft(2 * k + 1, j));
        ComplexD diff = ComplexD(mft(2 * k, j)) -
                        ComplexD(mft(2 * k + 1, j));
        if (negFirst) {
          diff = -diff;
        }

        axpy(rbTemp, sum, e, rbTemp);
        axpy(rbTempNeg, diff / lam_D, e, rbTempNeg);
      }

      // sol = (norm/pairScale) * [ SUM-channel - i * Meooe(DIFF/lam-channel) ]
      mat.Meooe(rbTempNeg, rbFermNeg);
      rbFermNeg = ComplexD(0., -1.) * rbFermNeg;
      {
        setCheckerboard(sol, rbTemp);
        setCheckerboard(sol, rbFermNeg);
      }
      sol *= norm / pairScale;
      if (subGuess) {
        sol = source - sol;
      }
    };
  };

  // one <name>[_<gamma>]_t<t> / ..._subtract pair per timeslice per
  // gamma; the loop recomputes exactly the name family getOutput()
  // declares at push time (FermionFlow pattern -- nothing is cached
  // between calls)
  for (unsigned int g = 0; g < gammas.size(); ++g) {
    std::string suffix =
        (gammas.size() > 1) ? ("_" + StagGamma::GetName(gammas[g])) : "";
    for (int t = 0; t < nt; ++t) {
      std::string name = getName() + suffix + "_t" + std::to_string(t);
      envCreate(Solver, name, Ls, makeSolver(g, t, false), mat);
      envCreate(Solver, name + "_subtract", Ls, makeSolver(g, t, true), mat);
    }
  }
}

/******************************************************************************
 *              TLowModeProjMesonField execution                              *
 ******************************************************************************/

template <typename FImpl, typename Pack>
void TLowModeProjMesonField<FImpl, Pack>::execute(void) {
  // Table-dimension checks and the self-check run here, not in setup():
  // the loader executes earlier in program order, so its table is filled
  // by the time this module executes -- but the scheduler's memory-
  // profiling pass dry-runs setup() on every module before anything
  // executes, when the table is still empty
  auto gammas = gammaList();
  auto mfs = mesonFieldList(gammas);
  auto &epack = envGet(Pack, par().lowModes);
  int nt = env().getDim().back();

  for (unsigned int g = 0; g < gammas.size(); ++g) {
    auto &mf =
        envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, mfs[g]);
    if (static_cast<unsigned int>(mf[0].rows()) !=
        2 * static_cast<unsigned int>(epack.evec.size())) {
      HADRONS_ERROR(Size, "meson field '" + mfs[g] + "' has " +
                              std::to_string(mf[0].rows()) +
                              " rows, expected 2*" +
                              std::to_string(epack.evec.size()) +
                              " (|e+o>/|e-o> pair layout)");
    }
    if (par().noiseIndex >= static_cast<unsigned int>(mf[0].cols())) {
      HADRONS_ERROR(Argument, "noiseIndex " +
                                  std::to_string(par().noiseIndex) +
                                  " out of range for '" + mfs[g] +
                                  "' (Nj = " + std::to_string(mf[0].cols()) +
                                  ")");
    }
  }

  // optional pairing/normalization self-check: the file pair sum of the
  // first eigenpair, summed over ALL timeslices, equals P * <e|eta_j,E>
  // with the live checkerboarded eigenvector; their ratio exposes the
  // production constants and catches pairing/order/normalization
  // mistakes. The live reference is gamma-independent, so the comparison
  // is only meaningful for the (G1,G1) table; other gammas get a skip
  // notice
  if (!par().noise.empty()) {
    auto &noise = envGet(std::vector<FermionField>, par().noise);
    if (par().noiseIndex >= noise.size()) {
      HADRONS_ERROR(Size, "noiseIndex out of range for noise object '" +
                              par().noise + "'");
    }
    int gIdentity = -1;
    for (unsigned int g = 0; g < gammas.size(); ++g) {
      if ((gammas[g].first == StagGamma::StagAlgebra::G1) &&
          (gammas[g].second == StagGamma::StagAlgebra::G1)) {
        gIdentity = g;
        break;
      }
    }
    if (gIdentity < 0) {
      LOG(Message) << "Self-check skipped: no (G1,G1) gamma in the list "
                      "(the live reference <e|eta> is gamma-independent)"
                   << std::endl;
    } else {
      auto &mf = envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>,
                        mfs[gIdentity]);
      unsigned int eigStart = par().eigStart;
      RealD pairScale = std::sqrt(2.0);
      if (!par().pairScale.empty()) {
        pairScale = strToVec<RealD>(par().pairScale)[0];
      }

      FermionField rbNoise(envGetRbGrid(FermionField));
      int cb = epack.evec[0].Checkerboard();

      rbNoise = Zero();
      rbNoise.Checkerboard() = cb;
      pickCheckerboard(cb, rbNoise, noise[par().noiseIndex]);

      ComplexD ipFull =
          TensorRemove(innerProduct(epack.evec[eigStart], rbNoise));
      ComplexD sumFile = 0.;
      for (int t = 0; t < nt; ++t) {
        sumFile += ComplexD(mf[t](2 * eigStart, par().noiseIndex)) +
                   ComplexD(mf[t](2 * eigStart + 1, par().noiseIndex));
      }
      if (std::abs(ipFull) > 1.e-12) {
        ComplexD pLive = sumFile / ipFull;
        LOG(Message) << "Self-check (gamma '"
                     << StagGamma::GetName(gammas[gIdentity])
                     << "'): file-derived production constant P = " << pLive
                     << " (configured pairScale = " << pairScale
                     << ")" << std::endl;
        if (std::abs(pLive - pairScale) > 0.05 * std::abs(pairScale)) {
          LOG(Warning) << "Meson-field pair normalization mismatch: "
                          "derived P = " << pLive << " but pairScale = "
                       << pairScale
                       << " -- check the |e+o>/|e-o> pair ordering and the "
                          "production normalization"
                       << std::endl;
        }
      } else {
        LOG(Warning) << "Self-check skipped: live inner product "
                        "<e_eigStart|eta_noiseIndex> vanishes"
                     << std::endl;
      }
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MSolver_LowModeProjMesonField_hpp_
