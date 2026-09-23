/*
 * LMAMesonFieldProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MFermion_LMAMesonFieldProp_hpp_
#define HadronsMILC_MFermion_LMAMesonFieldProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <A2AMatrix.hpp>
#include <EigenPack.hpp>
#include <Modules/MContraction/MesonField.hpp>
#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *        Low mode propagators from precomputed meson fields                *
 ******************************************************************************/
/*  GaugeProp-style replacement of the former file-driven LMA solver
    family: instead of registering lazy Solver closures consumed
    one at a time through GaugeProp middlemen, this module EAGERLY produces
    one scalar PropagatorField per (timeslice, gamma) output name during
    its own execute(); downstream modules fetch the field objects directly
    by name (no middleman, no solver indirection).

    The file's eigenvector rows come in |e+o>/|e-o> pairs (rows 2k/2k+1 of
    eigenpair k) -- produced identically by EigenPackFullPairs (full-volume)
    or by MesonField's checkerboarded cbPairs path. Their sum and
    difference reconstruct the parity-split inner products LowModeProj
    computes live:

      SUM_k(t,j)  = M[t][2k][j] + M[t][2k+1][j] ~ <e_E^k|eta_j,E>(t)
      DIFF_k(t,j) = M[t][2k][j] - M[t][2k+1][j]
                    ~ (i/lam_k) <Meooe(e_E^k)|eta_j,O>(t)

    giving (see the 2026-09-17 design artifact, decision D8, for the full
    derivation):

      ferm_c = (norm_cb / pairScale) *
               [ SUM_k e_k  -  i * Meooe( SUM_k (DIFF_k/lam_k) e_k ) ]

    with ferm_c placed in COLOR SLOT c of noise n's output propagator:
    each color-diluted noise source occupies 3 adjacent table columns
    (one per color index, TimeDilutedSpinColorDiagonal color-fast
    layout), global noise g = noiseIndex + n in window [3g, 3g + 3);
    the module reconstructs one FermionField per column j = 3g + c and
    assembles them FermToProp-style -- the GaugeProp solveField
    color-loop pattern applied to meson-field columns, run once per
    noise window (per-noise eigenpass loops: each noise's propagator is
    isolated; averaging over noises is deferred to the
    MContraction::Meson vector contract, whose srcCorrs hold the
    per-source correlators and corr their average).

    action      Staggered action module (Meooe parity move)
    lowModes    MassShiftEigenPack module with the CHECKERBOARDED
                eigenvectors/eigenvalues (row pair k <-> evec[k], eval[k])
    spinTaste   SpinTasteParams (the standard gammas/gauge/applyG5 struct
                shared with GaugeProp, Meson and MesonField). gammas ""
                (default) or a standard spin-taste pair list
                "(G1 G1) (G5 G5)": one output family per gamma; empty
                means the legacy single (G1,G1) family with bare names.
                Parsed TWICE the GaugeProp/Meson way (GaugeProp.hpp
                parseGammas): the RAW pairs are the output LABELS
                (GaugeProp's per-gamma guess grammar -- a StagGaugeProp
                with the same gammas and guess "<name>_t<t>" derives
                exactly these suffixes from its own parameters); the
                applyG5-conjugated pairs are the GAMMA KEYS that must
                match each loaded file's metadata (the producing
                StagA2AMesonField writes gamma_spin/gamma_taste under
                its own conjugated gammas, MesonField.hpp:532). With
                applyG5=true the file pairing permutes while the labels
                stay on the raw list. gauge must be empty: this module
                applies no spin-taste operator, so a gauge field would
                be unused
    mesonField  whitespace-separated LoadMesonField module names, one per
                gamma (positional parallel list, gammas[i] <-> entry i);
                each loader's published MesonFieldMILCMetadata side object
                ("<loader>_metadata") is cross-checked against the
                applyG5-conjugated VALUE of gammas[i] at execute time
                (mesonField[i] must load the file produced under the
                conjugated i-th gamma -- with applyG5=true the file
                pairing permutes while the labels stay on the raw list),
                so a miswired gamma/file pairing fails loudly instead of
                silently mislabeling output names
    noiseIndex  index of the FIRST noise window (unit = NOISES, not
                columns): the module reconstructs noises
                noiseIndex..noiseIndex+nNoise-1, window n reading
                columns (noiseIndex+n)*3 .. (noiseIndex+n)*3+2 (one
                column per color). The table must COVER the windows:
                cols >= (noiseIndex+nNoise)*3 (execute-time check)
    nNoise      number of noise sources reconstructed per output name
                (>= 1). nNoise == 1 publishes a scalar PropagatorField
                per output name (legacy grammar); nNoise > 1 publishes
                a std::vector<PropagatorField> of length nNoise,
                noise-major (entry n = noise n's propagator) -- the
                RandomWall scalar-collapse precedent (ef6d6d0)
    tA          first timeslice to produce (inclusive)
    tB          last timeslice to produce (inclusive, must be < nt)
    tStep       timeslice stride (>= 1); one PropagatorField per gamma per
                t in [tA, tB] with stride tStep, named "<name>_t<t>" for
                a single gamma (or an empty gammas list) and
                "<name>_t<t>_<spin>_<taste>" when the list holds more
                than one gamma (gamma segment LAST, mirroring GaugeProp's
                per-gamma suffix convention; built from the RAW gamma
                label, independent of applyG5, so GaugeProp guess
                lookups match)
    eigStart    first eigenpair to include (pair space)
    nEigs       number of eigenpairs (< 1: all)
    negFirst    ""/"false" (default): row 2k is |e+o>; "true": |e-o> comes
                first (flips the DIFF sign)
    pairScale   production normalization constant P (default sqrt(2); use
                1.0 for unit-norm full-volume Lanczos files)
    noise       optional name of the noise-vector environment object
                (std::vector<FermionField>) used for a one-shot pairing/
                normalization self-check at execute time, run once per
                noise window n against noise[(noiseIndex + n)*3]
                (full-time-extent sum; only the (G1,G1) table is
                checkable against the live, gamma-independent reference;
                other gammas log a skip)

    setup() only parses parameters, checks container sizes and allocates
    and zeroes the outputs (the scheduler dry-runs it for memory profiling
    while the loader tables are still empty 0x0 matrices); ALL table-content
    checks (row layout, color-window bounds, metadata cross-check,
    self-check) live at the top of execute(). The _subtract variants of
    the former solver family are dropped (they were never end-to-end
    validated and the producer contract has no caller-supplied source).
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class LMAMesonFieldPropMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LMAMesonFieldPropMILCPar,
                                  std::string,   action,
                                  std::string,   lowModes,
                                  std::string,   mesonField,
                                  SpinTasteParams, spinTaste,
                                  unsigned int,  noiseIndex,
                                  unsigned int,  nNoise,
                                  unsigned int,  tA,
                                  unsigned int,  tB,
                                  unsigned int,  tStep,
                                  unsigned int,  eigStart,
                                  int,           nEigs,
                                  std::string,   negFirst,
                                  std::string,   pairScale,
                                  std::string,   noise);
  LMAMesonFieldPropMILCPar(void)
      : tStep(1), nNoise(1), negFirst(""), pairScale("") {}
};

template <typename FImpl, typename Pack>
class TLMAMesonFieldPropMILC : public Module<LMAMesonFieldPropMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

private:
  // spinTaste.gammas parsed TWICE, the GaugeProp/Meson way
  // (GaugeProp.hpp parseGammas, Meson.hpp _mapSinkGammas): the RAW
  // parse gives the LABELS (output-name suffixes; GaugeProp guess
  // objects are looked up by exactly these names), the applyG5-
  // conjugated parse gives the VALUES (the gamma keys that must match
  // the loaded file's metadata). Empty gammas -> the legacy single
  // (G1,G1) default, bare names. Order preserved: mesonField is a
  // positional parallel list over the raw order; duplicate LABELS are
  // fatal (output names would collide)
  std::vector<std::pair<std::string, StagGamma::SpinTastePair>>
  parseGammas(void) const;
  // one LoadMesonField module name per gamma (positional parallel list,
  // strToVec<std::string>); count mismatch vs the gamma list is fatal
  // before any positional access
  std::vector<std::string> mesonFieldList(
      const std::vector<std::pair<std::string,
                                  StagGamma::SpinTastePair>> &gammas) const;
  // the timeslices this instance materializes: t in [tA, tB] stride
  // tStep, clamped to the lattice time extent (single source shared by
  // getOutput()/setup()/execute() -- the name family cannot diverge
  // between the three sites)
  std::vector<unsigned int> sliceTimes(void) const;
  // canonical output object name for gamma index g at timeslice t
  std::string outputName(const unsigned int g, const unsigned int t) const;
  // the full (gamma-major) output name family
  std::vector<std::string> outputNames(void) const;

public:
  // constructor
  TLMAMesonFieldPropMILC(const std::string name);
  // destructor
  virtual ~TLMAMesonFieldPropMILC(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(StagLMAMesonFieldProp,
                    ARG(TLMAMesonFieldPropMILC<STAGIMPL,
                                               MassShiftEigenPack<STAGIMPL>>),
                    MFermion);

/******************************************************************************
 *                TLMAMesonFieldPropMILC implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TLMAMesonFieldPropMILC<FImpl, Pack>::TLMAMesonFieldPropMILC(
    const std::string name)
    : Module<LMAMesonFieldPropMILCPar>(name) {}

// gamma/meson-field parallel lists //////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::pair<std::string, StagGamma::SpinTastePair>>
TLMAMesonFieldPropMILC<FImpl, Pack>::parseGammas(void) const {
  std::vector<StagGamma::SpinTastePair> keys, vals;
  if (par().spinTaste.gammas.empty()) {
    keys.push_back(std::make_pair(StagGamma::StagAlgebra::G1,
                                  StagGamma::StagAlgebra::G1));
    vals = keys;
  } else {
    vals = StagGamma::ParseSpinTasteString(par().spinTaste.gammas,
                                           par().spinTaste.applyG5);
    keys = StagGamma::ParseSpinTasteString(par().spinTaste.gammas);
  }
  std::vector<std::pair<std::string, StagGamma::SpinTastePair>> gammas;
  std::vector<std::string> names;
  for (unsigned int i = 0; i < keys.size(); ++i) {
    const std::string label = StagGamma::GetName(keys[i]);
    for (auto &n : names) {
      if (n == label) {
        HADRONS_ERROR(Argument, "duplicate gamma '" + label +
                                    "' in spinTaste.gammas (output names "
                                    "would collide)");
      }
    }
    names.push_back(label);
    gammas.push_back(std::make_pair(label, vals[i]));
  }

  return gammas;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLMAMesonFieldPropMILC<FImpl, Pack>::mesonFieldList(
    const std::vector<std::pair<std::string, StagGamma::SpinTastePair>>
        &gammas) const {
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

// output name family /////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<unsigned int> TLMAMesonFieldPropMILC<FImpl, Pack>::sliceTimes(
    void) const {
  std::vector<unsigned int> ts;
  // defensive: a zero tStep would loop forever below; setup() rejects
  // it loudly, but getOutput() may run before setup(), so guard here
  // too and simply enumerate nothing
  if (par().tStep < 1) {
    return ts;
  }
  const unsigned int nt = static_cast<unsigned int>(env().getDim().back());
  // clamp to the time extent: getOutput() runs before setup() validates
  // the range
  const unsigned int tLast = std::min(par().tB, nt - 1);
  for (unsigned int i = 0;; ++i) {
    const unsigned int t = par().tA + i * par().tStep;
    if ((t < par().tA) || (t > tLast)) {
      break;
    }
    ts.push_back(t);
  }

  return ts;
}

template <typename FImpl, typename Pack>
std::string TLMAMesonFieldPropMILC<FImpl, Pack>::outputName(
    const unsigned int g, const unsigned int t) const {
  auto gammas = parseGammas();
  // suffix from the RAW gamma LABEL (GaugeProp's parseGammas key
  // grammar): a StagGaugeProp with the same spinTaste.gammas and guess
  // "<name>_t<t>" derives exactly this name from its own parameters.
  // Single gamma -> bare names (legacy grammar); multiple gammas ->
  // trailing "_<spin>_<taste>" segment AFTER the timeslice (GaugeProp
  // appends gamma keys last when naming per-gamma objects)
  std::string suffix = (gammas.size() > 1) ? ("_" + gammas[g].first) : "";

  return getName() + "_t" + std::to_string(t) + suffix;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLMAMesonFieldPropMILC<FImpl, Pack>::outputNames(
    void) const {
  std::vector<std::string> names;
  auto gammas = parseGammas();
  auto ts = sliceTimes();
  for (unsigned int g = 0; g < gammas.size(); ++g) {
    for (auto &t : ts) {
      names.push_back(outputName(g, t));
    }
  }

  return names;
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TLMAMesonFieldPropMILC<FImpl, Pack>::getInput(void) {
  auto mfs = mesonFieldList(parseGammas());
  std::vector<std::string> in{par().action, par().lowModes};
  for (auto &mf : mfs) {
    in.push_back(mf);
    // explicit edge on the loader's metadata side object (GaugeProp
    // guess-object pattern): naming an env object in getInput() keeps
    // it alive until this module executes, so the execute-time
    // envGet(<loader>_metadata) cross-check cannot hit a GC-freed object
    in.push_back(mf + "_metadata");
  }
  if (!par().noise.empty()) {
    in.push_back(par().noise);
  }

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLMAMesonFieldPropMILC<FImpl, Pack>::getOutput(void) {
  return outputNames();
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TLMAMesonFieldPropMILC<FImpl, Pack>::setup(void) {
  int Ls = env().getObjectLs(par().action);
  if (Ls > 1) {
    HADRONS_ERROR(Argument, "Ls > 1 not implemented");
  }

  // this module applies no spin-taste operator (the gammas only pair
  // files with output labels), so a gauge field would be unused
  if (!par().spinTaste.gauge.empty()) {
    HADRONS_ERROR(Argument,
                  "spinTaste.gauge must be empty (this module applies no "
                  "spin-taste operator; gauge '" + par().spinTaste.gauge +
                  "' would be unused)");
  }

  // optional string parameters (missing XML nodes are tolerated for
  // strings only -- useStencil precedent, commit 297b724)
  if ((par().negFirst != "") && (par().negFirst != "false") &&
      (par().negFirst != "true")) {
    HADRONS_ERROR(Argument, "negFirst must be '', 'false' or 'true' (got '" +
                                par().negFirst + "')");
  }
  if (!par().pairScale.empty()) {
    auto scale = strToVec<RealD>(par().pairScale);
    if (scale.size() != 1) {
      HADRONS_ERROR(Argument, "pairScale must be a single number (got '" +
                                  par().pairScale + "')");
    }
  }

  // timeslice range: inclusive [tA, tB] with stride tStep >= 1 (the
  // SeqGamma tA/tB idiom, plus the first step parameter in this
  // codebase)
  int nt = env().getDim().back();
  if (par().tStep < 1) {
    HADRONS_ERROR(Argument, "tStep must be >= 1 (got " +
                                std::to_string(par().tStep) + ")");
  }
  // noise count: one propagator per noise window; nNoise == 1 keeps the
  // legacy scalar outputs (RandomWall scalar-collapse precedent)
  if (par().nNoise < 1) {
    HADRONS_ERROR(Argument, "nNoise must be >= 1 (got " +
                                std::to_string(par().nNoise) + ")");
  }
  if (par().tA > par().tB) {
    HADRONS_ERROR(Argument, "tA (" + std::to_string(par().tA) +
                                ") must not exceed tB (" +
                                std::to_string(par().tB) + ")");
  }
  if (par().tB >= static_cast<unsigned int>(nt)) {
    HADRONS_ERROR(Argument, "tB (" + std::to_string(par().tB) +
                                ") out of range: the lattice has " +
                                std::to_string(nt) + " timeslices");
  }

  auto gammas = parseGammas();
  auto mfs = mesonFieldList(gammas);

  LOG(Message) << "Setting up meson-field driven low mode propagator '"
               << getName() << "' for action '" << par().action
               << "' using eigenvectors from '" << par().lowModes
               << "' (noise index " << par().noiseIndex << ", nNoise "
               << par().nNoise << ", applyG5 "
               << (par().spinTaste.applyG5 ? "true" : "false")
               << ", one propagator per timeslice per gamma in [tA="
               << par().tA << ", tB=" << par().tB << "] with stride "
               << par().tStep << "):" << std::endl;
  for (unsigned int g = 0; g < gammas.size(); ++g) {
    LOG(Message) << "  gamma '" << gammas[g].first << "' (file key '"
                 << StagGamma::GetName(gammas[g].second) << "') from '"
                 << mfs[g] << "'" << std::endl;
  }

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

  // meson-field tables: one per gamma; container-size check only (the
  // loader fills the nt-length vector in its own setup, so this is
  // dry-run safe). Row/column/metadata CONTENT checks live in execute()
  for (auto &mfName : mfs) {
    auto &mf = envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, mfName);
    if (static_cast<int>(mf.size()) != nt) {
      HADRONS_ERROR(Size, "meson field '" + mfName + "' has " +
                              std::to_string(mf.size()) +
                              " timeslices, expected " +
                              std::to_string(nt));
    }
  }

  // temps: the per-color-column checkerboard accumulators of the
  // reconstruction (SUM channels rbTemp0..2, DIFF channels
  // rbTempNeg0..2), one Meooe target and one full-grid assembly target.
  // envTmp, not the former envCache: the eager module consumes them only
  // inside its own execute(). The six accumulators exist so that ONE
  // element-wise pass per eigenvector can update every color column at
  // once (the columns differ only in their table coefficients)
  static_assert(FImpl::Dimension == 3,
                "color-fused reconstruction assumes three color columns");
  envTmp(FermionField, "rbTemp0", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTemp1", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTemp2", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTempNeg0", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTempNeg1", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTempNeg2", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbFermNeg", 1, envGetRbGrid(FermionField));
  envTmpLat(FermionField, "sol");

  // allocation-only output creation (the GaugeProp setupHelper
  // contract): one zeroed propagator per output name -- a scalar
  // PropagatorField when nNoise == 1 (legacy consumers and guess wiring
  // address scalar fields), a std::vector<PropagatorField> of length
  // nNoise otherwise (MContraction::Meson's element-wise vector contract
  // consumes it: srcCorrs per source, corr the average). RandomWall
  // scalar-collapse precedent (ef6d6d0). Sizing from par().nNoise --
  // never from table contents (dry-run safety)
  for (auto &name : outputNames()) {
    if (par().nNoise == 1) {
      envCreate(PropagatorField, name, 1, envGetGrid(PropagatorField));
      envGet(PropagatorField, name) = Zero();
    } else {
      envCreate(std::vector<PropagatorField>, name, 1, par().nNoise,
                envGetGrid(PropagatorField));
      auto &props = envGet(std::vector<PropagatorField>, name);
      for (auto &p : props) {
        p = Zero();
      }
    }
  }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TLMAMesonFieldPropMILC<FImpl, Pack>::execute(void) {
  // content checks run here, not in setup(): the scheduler's memory-
  // profiling pass dry-runs every module's setup() BEFORE anything
  // executes, when the loader tables are still nt empty 0x0 matrices
  auto gammas = parseGammas();
  auto mfs = mesonFieldList(gammas);
  auto &mat = envGet(FMat, par().action);
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
    // the table must COVER the requested noise windows: window n reads
    // columns (noiseIndex + n)*Nc .. (noiseIndex + n)*Nc + Nc-1, so the
    // noise index needs (noiseIndex + nNoise)*Nc columns. Superset
    // tables are legal (e.g. an all-timeslices time-diluted table:
    // window 0 is its t=0 noise)
    if (static_cast<unsigned int>(mf[0].cols()) <
        (par().noiseIndex + par().nNoise) * FImpl::Dimension) {
      HADRONS_ERROR(Size,
                    "meson field '" + mfs[g] + "' has " +
                        std::to_string(mf[0].cols()) +
                        " columns, too few for noise windows [noiseIndex=" +
                        std::to_string(par().noiseIndex) +
                        ", noiseIndex+nNoise=" +
                        std::to_string(par().noiseIndex + par().nNoise) +
                        ") (needs " +
                        std::to_string((par().noiseIndex + par().nNoise) *
                                       FImpl::Dimension) +
                        " columns: (noiseIndex+nNoise)*3)");
    }
    // gamma/file pairing cross-check against the loader's published
    // metadata, on the applyG5-conjugated VALUE (the producing
    // StagA2AMesonField writes gamma_spin/gamma_taste under its own
    // conjugated gammas, MesonField.hpp:532): the raw label never
    // appears in the file. A miswired pairing would otherwise just
    // mislabel the output names
    auto &md = envGet(MContraction::MesonFieldMILCMetadata,
                      mfs[g] + "_metadata");
    if ((md.gamma_spin != gammas[g].second.first) ||
        (md.gamma_taste != gammas[g].second.second)) {
      // distinguish a never-filled side object (non-HDF5 build or
      // legacy file) from a genuine gamma/file miswire
      const std::string fileSt =
          StagGamma::GetName(md.gamma_spin, md.gamma_taste);
      if (fileSt.find("undef") != std::string::npos) {
        HADRONS_ERROR(Argument,
                      "meson field '" + mfs[g] + "' carries undefined "
                      "spin-taste metadata (non-HDF5 build or legacy "
                      "file): cannot cross-check gamma '" +
                          gammas[g].first + "'");
      }
      HADRONS_ERROR(Argument,
                    "gamma '" + gammas[g].first + "' (file key '" +
                        StagGamma::GetName(gammas[g].second) +
                        "' under applyG5 " +
                        (par().spinTaste.applyG5 ? "true" : "false") +
                        ") is configured for meson field '" + mfs[g] +
                        "', but the file holds spin-taste '" + fileSt +
                        "'");
    }
  }

  // optional pairing/normalization self-check (unchanged from the former
  // solver family): the file pair sum of the first eigenpair, summed
  // over ALL timeslices, equals P * <e|eta_j,E> with the live
  // checkerboarded eigenvector; their ratio exposes the production
  // constants and catches pairing/order/normalization mistakes. The
  // live reference is gamma-independent, so the checkable table is the
  // one whose file CONTENT is the identity pairing: scan the
  // applyG5-conjugated VALUES (with applyG5=true the identity table is
  // reached through its conjugated file key); other gammas get a skip
  // notice
  if (!par().noise.empty()) {
    auto &noise = envGet(std::vector<FermionField>, par().noise);
    if ((par().noiseIndex + par().nNoise) * FImpl::Dimension >
        noise.size()) {
      HADRONS_ERROR(Size, "noise windows [noiseIndex, noiseIndex+nNoise) "
                          "out of range for noise object '" +
                              par().noise + "'");
    }
    int gIdentity = -1;
    for (unsigned int g = 0; g < gammas.size(); ++g) {
      if ((gammas[g].second.first == StagGamma::StagAlgebra::G1) &&
          (gammas[g].second.second == StagGamma::StagAlgebra::G1)) {
        gIdentity = g;
        break;
      }
    }
    if (gIdentity < 0) {
      LOG(Message) << "Self-check skipped: no (G1,G1) gamma VALUE in "
                      "the list (the live reference <e|eta> is "
                      "gamma-independent)"
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

      // one self-check per noise window (D5): window n checks
      // noise[(noiseIndex + n)*3] against the same table column pair,
      // keeping the object-index == column convention (window-unit
      // noiseIndex: global noise g = noiseIndex + n)
      for (unsigned int n = 0; n < par().nNoise; ++n) {
        const unsigned int j =
            (par().noiseIndex + n) * FImpl::Dimension;

        rbNoise = Zero();
        rbNoise.Checkerboard() = cb;
        pickCheckerboard(cb, rbNoise, noise[j]);

        ComplexD ipFull =
            TensorRemove(innerProduct(epack.evec[eigStart], rbNoise));
        ComplexD sumFile = 0.;
        for (int t = 0; t < nt; ++t) {
          sumFile += ComplexD(mf[t](2 * eigStart, j)) +
                     ComplexD(mf[t](2 * eigStart + 1, j));
        }
        // |ip| through .real()/.imag(): std::abs has no overload for the
        // ComplexD (thrust::complex) of GPU builds
        if (std::hypot(ipFull.real(), ipFull.imag()) > 1.e-12) {
          ComplexD pLive = sumFile / ipFull;
          LOG(Message) << "Self-check (gamma '" << gammas[gIdentity].first
                       << "', file key '"
                       << StagGamma::GetName(gammas[gIdentity].second)
                       << "', noise window " << n << "): file-derived "
                       << "production constant P = " << pLive
                       << " (configured pairScale = " << pairScale
                       << ")" << std::endl;
          // modulus via .real()/.imag(): std::abs has no ComplexD
          // overload on GPU builds (see EigenPackCheck)
          const ComplexD dP = pLive - static_cast<RealD>(pairScale);
          if (std::hypot(dP.real(), dP.imag()) >
              0.05 * std::abs(pairScale)) {
            LOG(Warning) << "Meson-field pair normalization mismatch "
                            "(noise window " << n << "): derived P = "
                         << pLive << " but pairScale = " << pairScale
                         << " -- check the |e+o>/|e-o> pair ordering and "
                            "the production normalization"
                         << std::endl;
          }
        } else {
          LOG(Warning) << "Self-check skipped (noise window " << n
                       << "): live inner product <e_eigStart|eta_" << j
                       << "> vanishes" << std::endl;
        }
      }
    }
  }

  // reconstruction: one output per (t, gamma) name -- a scalar
  // PropagatorField when nNoise == 1, nNoise propagators otherwise.
  // Color slot c of noise n is reconstructed from the color-diluted
  // source's table column j = (noiseIndex + n)*3 + c (adjacent columns
  // = adjacent colors, TimeDilutedSpinColorDiagonal color-fast layout)
  // and assembled FermToProp-style -- the GaugeProp solveField
  // color-loop pattern applied to meson-field columns, once per noise
  // window
  unsigned int eigStart = par().eigStart;
  int nEigs = par().nEigs;
  if (nEigs < 1) {
    nEigs = epack.evec.size();
  }
  bool negFirst = (par().negFirst == "true");
  RealD pairScale = std::sqrt(2.0);
  if (!par().pairScale.empty()) {
    pairScale = strToVec<RealD>(par().pairScale)[0];
  }

  envGetTmp(FermionField, rbTemp0);
  envGetTmp(FermionField, rbTemp1);
  envGetTmp(FermionField, rbTemp2);
  envGetTmp(FermionField, rbTempNeg0);
  envGetTmp(FermionField, rbTempNeg1);
  envGetTmp(FermionField, rbTempNeg2);
  envGetTmp(FermionField, rbFermNeg);
  envGetTmp(FermionField, sol);

  // SUM/DIFF accumulator columns (entry c reconstructs table column
  // nBase + c of the current noise window); array sugar over the named
  // env temps for the per-column post-processing loop below
  FermionField *rbTempC[FImpl::Dimension] = {&rbTemp0, &rbTemp1, &rbTemp2};
  FermionField *rbTempNegC[FImpl::Dimension] = {&rbTempNeg0, &rbTempNeg1,
                                                &rbTempNeg2};

  int cb = epack.evec[0].Checkerboard();
  RealD norm = 1. / ::sqrt(norm2(epack.evec[0]));

  auto ts = sliceTimes();
  for (unsigned int g = 0; g < gammas.size(); ++g) {
    auto &mf =
        envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, mfs[g]);
    for (auto &t : ts) {
      const A2AMatrix<HADRONS_A2AM_IO_TYPE> &mft = mf[t];
      // scalar collapse: nNoise == 1 publishes a bare PropagatorField
      // (legacy consumers); otherwise a noise-major vector consumed
      // element-wise by MContraction::Meson (source[i] <-> sink[i])
      std::vector<PropagatorField> *propVec = nullptr;
      PropagatorField *propScalar = nullptr;
      if (par().nNoise == 1) {
        propScalar = &envGet(PropagatorField, outputName(g, t));
      } else {
        propVec =
            &envGet(std::vector<PropagatorField>, outputName(g, t));
      }
      // per-noise reconstruction (per-noise eigenpass loops): noise
      // window n is global noise noiseIndex + n, occupying columns
      // (noiseIndex+n)*3 .. (noiseIndex+n)*3 + 2; its three color
      // columns reconstruct that noise's propagator in isolation -- the
      // accumulators, kernel passes and Meooe post-processing below
      // run once per noise
      for (unsigned int n = 0; n < par().nNoise; ++n) {
        PropagatorField &prop =
            (par().nNoise == 1) ? *propScalar : (*propVec)[n];
        const unsigned int nBase =
            (par().noiseIndex + n) * FImpl::Dimension;
        // load-bearing zero-init: the color-only FermToProp below
        // writes ONLY color slot c, so every other slot must start
        // zeroed
        prop = Zero();

        // zero the six per-column accumulators (SUM channels rbTempC,
        // DIFF channels rbTempNegC): one fused pass per eigenvector below
        // updates every color column, so all accumulators start clean
        // before the single eigenpair loop
        for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
          *rbTempC[c] = Zero();
          rbTempC[c]->Checkerboard() = cb;
          *rbTempNegC[c] = Zero();
          rbTempNegC[c]->Checkerboard() = cb;
        }

        // accumulate the two parity channels from the file row pairs; the
        // subtraction channel is accumulated as DIFF/lam so that the live
        // Meooe below completes the 1/lam_D^2 weighting of LowModeProj.
        // ALL color columns are built from the SAME eigenvectors, so one
        // fused element-wise pass per eigenpair updates all six
        // accumulators: the per-element expressions replicate the former
        // per-column axpy calls verbatim (bit-identical arithmetic) while
        // streaming e once for all three columns -- a third of the kernel
        // invocations and eigenvector traffic of the per-column version.
        // The hoisted aliased Read+Write view pairs mirror Grid's own
        // in-place axpy idiom; the extra scope closes every view before
        // the Meooe/setCheckerboard sequence below: those open CPU-mode
        // views, and the Grid memory manager asserts against mixing
        // concurrent view families on one buffer
        {
          autoView(t0W_v, rbTemp0, AcceleratorWrite);
          autoView(t0R_v, rbTemp0, AcceleratorRead);
          autoView(t1W_v, rbTemp1, AcceleratorWrite);
          autoView(t1R_v, rbTemp1, AcceleratorRead);
          autoView(t2W_v, rbTemp2, AcceleratorWrite);
          autoView(t2R_v, rbTemp2, AcceleratorRead);
          autoView(n0W_v, rbTempNeg0, AcceleratorWrite);
          autoView(n0R_v, rbTempNeg0, AcceleratorRead);
          autoView(n1W_v, rbTempNeg1, AcceleratorWrite);
          autoView(n1R_v, rbTempNeg1, AcceleratorRead);
          autoView(n2W_v, rbTempNeg2, AcceleratorWrite);
          autoView(n2R_v, rbTempNeg2, AcceleratorRead);
          // per-eigenpair table coefficients (expressions identical to the
          // former inline computation -- value- and order-identical)
          auto coeffs = [&](const int k, ComplexD sumC[], ComplexD negC[]) {
            const RealD lam_D = epack.eval[k].imag();
            for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
              const unsigned int j = nBase + c;
              ComplexD sum =
                  ComplexD(mft(2 * k, j)) + ComplexD(mft(2 * k + 1, j));
              ComplexD diff =
                  ComplexD(mft(2 * k, j)) - ComplexD(mft(2 * k + 1, j));
              if (negFirst) {
                diff = -diff;
              }
              sumC[c] = sum;
              negC[c] = diff / lam_D;
            }
          };

          // batches of four eigenpairs per kernel launch: thread dispatch
          // dominates these small kernels, and batching leaves the
          // per-element accumulation order exactly as before (k strictly
          // descending; batch statements in k order), so results stay
          // bit-identical. The trailing nEigs % 4 eigenpairs fall through
          // to the single-eigenpair kernel below
          int kHi = static_cast<int>(eigStart) + nEigs - 1;
          for (; kHi - 3 >= int(eigStart); kHi -= 4) {
            const int kA = kHi, kB = kHi - 1, kC = kHi - 2, kD = kHi - 3;
            const FermionField &eA = epack.evec[kA];
            const FermionField &eB = epack.evec[kB];
            const FermionField &eC = epack.evec[kC];
            const FermionField &eD = epack.evec[kD];

            ComplexD sumC[FImpl::Dimension], negC[FImpl::Dimension];
            coeffs(kA, sumC, negC);
            const ComplexD sA0 = sumC[0], sA1 = sumC[1], sA2 = sumC[2];
            const ComplexD dA0 = negC[0], dA1 = negC[1], dA2 = negC[2];
            coeffs(kB, sumC, negC);
            const ComplexD sB0 = sumC[0], sB1 = sumC[1], sB2 = sumC[2];
            const ComplexD dB0 = negC[0], dB1 = negC[1], dB2 = negC[2];
            coeffs(kC, sumC, negC);
            const ComplexD sC0 = sumC[0], sC1 = sumC[1], sC2 = sumC[2];
            const ComplexD dC0 = negC[0], dC1 = negC[1], dC2 = negC[2];
            coeffs(kD, sumC, negC);
            const ComplexD sD0 = sumC[0], sD1 = sumC[1], sD2 = sumC[2];
            const ComplexD dD0 = negC[0], dD1 = negC[1], dD2 = negC[2];

            autoView(eAR_v, eA, AcceleratorRead);
            autoView(eBR_v, eB, AcceleratorRead);
            autoView(eCR_v, eC, AcceleratorRead);
            autoView(eDR_v, eD, AcceleratorRead);
            accelerator_for(ss, eAR_v.size(),
                            FermionField::vector_type::Nsimd(), {
              auto evA = coalescedRead(eAR_v[ss]);
              auto evB = coalescedRead(eBR_v[ss]);
              auto evC = coalescedRead(eCR_v[ss]);
              auto evD = coalescedRead(eDR_v[ss]);
              // statements in k order (A,B,C,D): each accumulator's update
              // sequence per element matches the unbatched kernel exactly
              // (the R and W views alias one buffer, so later statements
              // read the value written by earlier ones)
              coalescedWrite(t0W_v[ss],
                             sA0 * evA + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sB0 * evB + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sC0 * evC + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sD0 * evD + coalescedRead(t0R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sA1 * evA + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sB1 * evB + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sC1 * evC + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sD1 * evD + coalescedRead(t1R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sA2 * evA + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sB2 * evB + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sC2 * evC + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sD2 * evD + coalescedRead(t2R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dA0 * evA + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dB0 * evB + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dC0 * evC + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dD0 * evD + coalescedRead(n0R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dA1 * evA + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dB1 * evB + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dC1 * evC + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dD1 * evD + coalescedRead(n1R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dA2 * evA + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dB2 * evB + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dC2 * evC + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dD2 * evD + coalescedRead(n2R_v[ss]));
            });
          }
          // trailing eigenpairs (nEigs % 4): the original single-eigenpair
          // kernel
          for (; kHi >= int(eigStart); --kHi) {
            const FermionField &e = epack.evec[kHi];

            ComplexD sumC[FImpl::Dimension], negC[FImpl::Dimension];
            coeffs(kHi, sumC, negC);
            const ComplexD s0 = sumC[0], s1 = sumC[1], s2 = sumC[2];
            const ComplexD d0 = negC[0], d1 = negC[1], d2 = negC[2];

            autoView(eR_v, e, AcceleratorRead);
            accelerator_for(ss, eR_v.size(),
                            FermionField::vector_type::Nsimd(), {
              auto ev = coalescedRead(eR_v[ss]);
              coalescedWrite(t0W_v[ss], s0 * ev + coalescedRead(t0R_v[ss]));
              coalescedWrite(t1W_v[ss], s1 * ev + coalescedRead(t1R_v[ss]));
              coalescedWrite(t2W_v[ss], s2 * ev + coalescedRead(t2R_v[ss]));
              coalescedWrite(n0W_v[ss], d0 * ev + coalescedRead(n0R_v[ss]));
              coalescedWrite(n1W_v[ss], d1 * ev + coalescedRead(n1R_v[ss]));
              coalescedWrite(n2W_v[ss], d2 * ev + coalescedRead(n2R_v[ss]));
            });
          }
        }

        // ferm_c = (norm/pairScale) *
        //          [ SUM-channel - i * Meooe(DIFF/lam-channel) ]
        // for each color column c, from the c-th accumulators above
        for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
          rbFermNeg = Zero();
          rbFermNeg.Checkerboard() = (cb == Even) ? Odd : Even;

          mat.Meooe(*rbTempNegC[c], rbFermNeg);

          // fused per-column assembly: ONE pass over the full-grid sol
          // replaces the former FOUR element-wise passes (the -i axpy
          // on rbFermNeg, the two setCheckerboard copies, and the sol
          // scale). Per element the op sequence is IDENTICAL to the
          // unfused sequence -- cb-parity sites: sum * scale; other
          // -parity sites: ((0,-1) * meooe) * scale -- with the same
          // scalar types and operand order the lattice-level ops
          // lowered to (Lattice *= lowers to (*this)*r, Lattice_base.h;
          // the complex scalar product lowers to the same tensor
          // operator* the eigenpass coefficients use), so results are
          // bit-identical. The full-grid -> rb-grid site mapping
          // replicates Grid's acceleratorSetCheckerboard
          // (Lattice_transfer.h): coordinate from _rdimensions, parity
          // from _checker_dim_mask, rb index from _ostride with the
          // checker dim halved. Kernel launches per color column drop
          // from five to two (assembly + FermToProp; Meooe unchanged)
          {
            const GridBase *halfGrid = rbFermNeg.Grid();
            const Coordinate rdimFull = sol.Grid()->_rdimensions;
            const Coordinate rdimHalf = halfGrid->_rdimensions;
            const Coordinate cbMask = halfGrid->_checker_dim_mask;
            const Coordinate ostride = halfGrid->_ostride;
            const int ndim = halfGrid->_ndimension;
            const RealD scale = norm / pairScale;
            const ComplexD negI(0., -1.);
            const int cbSum = cb;

            autoView(solW, sol, AcceleratorWrite);
            // named reference first: autoView(n, *ptr[c], m) expands to
            // *ptr[c].View(m) -- '.' binds tighter than '*' (ledger)
            FermionField &sumC = *rbTempC[c];
            autoView(sumR, sumC, AcceleratorRead);
            autoView(negR, rbFermNeg, AcceleratorRead);
            accelerator_for(ss, sol.Grid()->oSites(),
                            FermionField::vector_type::Nsimd(), {
              Coordinate coor;
              int linear = 0;
              Lexicographic::CoorFromIndex(coor, ss, rdimFull);
              for (int d = 0; d < ndim; ++d) {
                if (cbMask[d]) {
                  linear += coor[d];
                }
              }
              int ssh = 0;
              for (int d = 0; d < ndim; ++d) {
                if (d == 0) {
                  ssh += ostride[d] * ((coor[d] / 2) % rdimHalf[d]);
                } else {
                  ssh += ostride[d] * (coor[d] % rdimHalf[d]);
                }
              }
              if ((linear & 0x1) == cbSum) {
                coalescedWrite(solW[ss],
                               coalescedRead(sumR[ssh]) * scale);
              } else {
                coalescedWrite(solW[ss],
                               (negI * coalescedRead(negR[ssh])) * scale);
              }
            });
          }

          FermToProp<FImpl>(prop, sol, c);
        }

      }
      LOG(Message) << "Reconstructed '" << outputName(g, t) << "' ("
                   << par().nNoise << " noise window(s)) from columns "
                   << par().noiseIndex * FImpl::Dimension << ".."
                   << ((par().noiseIndex + par().nNoise) *
                       FImpl::Dimension - 1)
                   << " of '" << mfs[g] << "'" << std::endl;
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MFermion_LMAMesonFieldProp_hpp_
