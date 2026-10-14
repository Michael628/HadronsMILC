/*
 * ASlashOneLinkMesonField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MContraction_A2AASlashOneLinkMesonField_hpp_
#define HadronsMILC_MContraction_A2AASlashOneLinkMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <A2AVectors.hpp>
#include <EigenPack.hpp>
#include <A2AMatrix.hpp>
#include <A2AStencilKernel.hpp>
#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

// One-link J.A flavour of the stochastic A2A meson-field module: the
// A2AWorkerSpinTasteStencil J.A one-link worker. Each stochastic sample
// provides Nd RAW photon components A_k = PeekIndex(AField, k) -- NO
// applyCoeffsAndPhase fold, NO per-direction sum (the task dresses the
// gauge links itself, W~_k(x) = A_k(x) . U_k(x), and bakes the standard
// eta_k(x) . (1/2) phase: feeding folded/summed fields would double the
// coefficient). One output channel per component (J_kA_k); the per-sample
// loop constructs one worker per sample (the ctor takes exactly one
// sample's Nd fields). Output mu axis 4*nMom, g-major (mu = k*nMom + p).
// CB low modes via cbPairsLeft/cbPairsRight pair sources (the canonical
// ASlash module surface).
class ASlashOneLinkMesonFieldMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ASlashOneLinkMesonFieldMILCPar, int, block,
                                  std::string, lowModes, std::string, left,
                                  std::string, cbPairsLeft, std::string,
                                  cbPairsRight, std::string, right,
                                  std::string, output, std::string, gauge,
                                  std::string, EmFunc, int,
                                  nEmFields, std::string, EmSeedString,
                                  std::vector<std::string>, mom);
  ASlashOneLinkMesonFieldMILCPar() : nEmFields(0) {}
};
// gauge: name of the gauge-field module (REQUIRED: the one-link operator
//        transports through the photon-dressed link W~_k = A_k . U_k)
// (NO gammas parameter: the directional spec gammas are synthesized
//  inside Grid via jaOneLinkSpecs -- (gmu[k], G1) per direction)

class ASlashOneLinkMesonFieldMILCMetadata : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ASlashOneLinkMesonFieldMILCMetadata,
                                  std::vector<RealF>, momentum, std::string,
                                  amu_seed, int, amu_index);

  ASlashOneLinkMesonFieldMILCMetadata() : momentum{}, amu_index(-1), amu_seed("N/A") {}
};
// per file: momentum = the file's momentum vector; amu_index = k (the
// J_kA_k direction, 0..3); amu_seed = "<EmSeedString>_<j>" (the
// sample-qualified seed string, so name<->metadata stay 1:1)

template <typename T, typename FImpl>
class ASlashOneLinkMesonFieldKernelMILC : public A2AStencilKernelMILC<T, FImpl> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  ASlashOneLinkMesonFieldKernelMILC(GridBase *grid,
                                    const std::vector<ComplexField> &mom,
                                    const std::vector<ComplexField> &Amu,
                                    LatticeGaugeField *U,
                                    int orthogDir)
      : A2AStencilKernelMILC<T, FImpl>(grid) {
    GridCartesian *cart = dynamic_cast<GridCartesian *>(grid);
    if (cart == nullptr) {
      HADRONS_ERROR(Implementation,
                    "ASlashOneLinkMesonField requires a Cartesian grid");
    }
    // J.A one-link: Amu holds exactly Nd RAW photon components (the task
    // dresses W~_k = A_k . U_k and bakes eta_k . (1/2) itself -- never
    // feed applyCoeffsAndPhase-folded fields).
    this->setWorker(cart, mom, Amu, U, orthogDir);
  }

  virtual ~ASlashOneLinkMesonFieldKernelMILC(void) = default;
};

template <typename FImpl, typename Pack>
class TASlashOneLinkMesonFieldMILC : public Module<ASlashOneLinkMesonFieldMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );
  typedef A2AMatrixBlockComputationMILC<
      Complex, FermionField, ASlashOneLinkMesonFieldMILCMetadata, HADRONS_A2AM_IO_TYPE>
      Computation;
  typedef ASlashOneLinkMesonFieldKernelMILC<Complex, FImpl> Kernel;
  typedef PhotonR::GaugeLinkField EmComp;
  typedef PhotonR::GaugeField EmField;

public:
  // constructor
  TASlashOneLinkMesonFieldMILC(const std::string name);
  // destructor
  virtual ~TASlashOneLinkMesonFieldMILC(void) {};

  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);

private:
  std::string _momphName, _EmName;
  std::vector<std::vector<Real>> _mom;
};

MODULE_REGISTER(StagA2AASlashOneLinkMesonField, ARG(TASlashOneLinkMesonFieldMILC<STAGIMPL, MassShiftEigenPack<STAGIMPL>>), MContraction);

/******************************************************************************
 *                  TASlashOneLinkMesonFieldMILC implementation *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TASlashOneLinkMesonFieldMILC<FImpl, Pack>::TASlashOneLinkMesonFieldMILC(
    const std::string name)
    : Module<ASlashOneLinkMesonFieldMILCPar>(name), _momphName(name + "_momph"),
      _EmName(name + "_em") {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TASlashOneLinkMesonFieldMILC<FImpl, Pack>::getInput(void) {
  std::vector<std::string> in = {};
  if (!par().left.empty())
    in.push_back(par().left);
  if (!par().right.empty())
    in.push_back(par().right);

  if (!par().lowModes.empty()) {
    if (!par().cbPairsLeft.empty()) {
      in.push_back(par().cbPairsLeft);
      in.push_back(par().cbPairsRight);
    }
    in.push_back(par().lowModes);
  }

  if (!par().EmFunc.empty()) {
    in.push_back(par().EmFunc);
  }

  if (par().gauge.empty()) {
    HADRONS_ERROR(Argument, "ASlashOneLinkMesonField requires the 'gauge' "
                            "module name");
  }
  in.push_back(par().gauge);

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TASlashOneLinkMesonFieldMILC<FImpl, Pack>::getOutput(void) {
  std::vector<std::string> out = {};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TASlashOneLinkMesonFieldMILC<FImpl, Pack>::setup(void) {
  if (par().EmFunc.empty()) {
    HADRONS_ERROR(Argument, "ASlashOneLinkMesonField requires the 'EmFunc' "
                            "stochastic photon module name");
  }
  if (par().gauge.empty()) {
    HADRONS_ERROR(Argument, "ASlashOneLinkMesonField requires the 'gauge' "
                            "module name (the one-link operator dresses the "
                            "gauge links: W~_k = A_k . U_k)");
  }

  _mom.clear();

  for (auto &pstr : par().mom) {
    auto p = strToVec<Real>(pstr);

    if (p.size() != env().getNd() - 1) {
      HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size()) +
                              " components instead of " +
                              std::to_string(env().getNd() - 1));
    }
    _mom.push_back(p);
  }
  int nmom = _mom.size();
  // All-zero single-momentum collapse to zero phase fields: scoped to the
  // single-entry case (a phase field is exp(0)=1; skipping it would
  // silently truncate the flat mu axis). NOTE: scoped from the legacy
  // module's verbatim form, whose outside-the-guard `allzero` collapsed
  // multi-momentum lists to zero phase fields (invisible on the legacy
  // path, fatal on the stencil path -- MesonField.hpp pattern).
  if (par().mom.size() == 1) {
    bool allzero = true;
    for (auto p : _mom[0]) {
      if (p != 0)
        allzero = false;
    }
    if (allzero)
      nmom = 0;
  }

  if (!par().cbPairsLeft.empty() != !par().cbPairsRight.empty()) {
    HADRONS_ERROR(Argument,
                  "ASlashOneLinkMesonField: 'cbPairsLeft' and 'cbPairsRight' "
                  "must be set together (two EigenPackCBPairs module "
                  "instances)");
  }
  if (!par().cbPairsLeft.empty() && par().cbPairsLeft == par().cbPairsRight) {
    HADRONS_ERROR(Argument, "ASlashOneLinkMesonField: 'cbPairsLeft' and "
                            "'cbPairsRight' must name distinct "
                            "EigenPackCBPairs instances");
  }
  if (!par().cbPairsLeft.empty() && par().lowModes.empty()) {
    HADRONS_ERROR(Argument, "ASlashOneLinkMesonField: 'cbPairsLeft'/"
                            "'cbPairsRight' require 'lowModes'");
  }
  if (!par().lowModes.empty() && par().cbPairsLeft.empty()) {
    // red-black eigenvectors without CB pair sources would take the
    // full-grid path and read each half-volume evec out of bounds
    auto &lowModes = envGet(Pack, par().lowModes);
    if (!lowModes.evec.empty() && lowModes.evec[0].Grid()->_isCheckerBoarded) {
      HADRONS_ERROR(Argument,
                    "ASlashOneLinkMesonField: 'lowModes' ('" + par().lowModes +
                        "') holds red-black checkerboard eigenvectors; set "
                        "'cbPairsLeft'/'cbPairsRight' (two "
                        "EigenPackCBPairs instances)");
    }
  }

  envCache(std::vector<ComplexField>, _momphName, 1, nmom,
           envGetGrid(ComplexField));

  // Nd raw photon component fields, refilled per stochastic sample
  envCache(std::vector<ComplexField>, _EmName, 1, env().getNd(),
           envGetGrid(ComplexField));
  envTmpLat(ComplexField, "coor");
  // Computation dims: ext = 4 components (J_kA_k), str = momenta (the
  // g-major flat-mu contract: linear slot e*nstr + s == k*nMom + p). One
  // Computation object serves every sample (execute re-entrant per call;
  // per-sample IO lambdas distinguish the output files).
  envTmp(Computation, "computation", 1, envGetGrid(FermionField),
         env().getNd() - 1, env().getNd(), _mom.size(), par().block, this);
  envTmp(std::vector<FermionField>, "dummy", 1, 0, envGetGrid(FermionField));
  envTmpLat(PhotonR::GaugeField, "AField");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TASlashOneLinkMesonFieldMILC<FImpl, Pack>::execute(void) {
  bool hasLowModes = (!par().lowModes.empty());
  bool isCheckerBoarded = (!par().cbPairsLeft.empty());

  std::vector<FermionField> *left, *right;

  envGetTmp(std::vector<FermionField>, dummy);
  if (!par().left.empty()) {
    left = &(envGet(std::vector<FermionField>, par().left));
  } else {
    left = &dummy;
  }
  if (!par().right.empty()) {
    right = &(envGet(std::vector<FermionField>, par().right));
  } else {
    right = &dummy;
  }

  int nt = env().getDim().back();
  int N_i = left->size();
  int N_j = right->size();

  if (hasLowModes) {
    auto &lowModes = envGet(Pack, par().lowModes);
    if (N_j != 0 && N_i == 0) {
      N_i += (isCheckerBoarded ? 2 : 1) * lowModes.evec.size();
    } else if (N_i != 0 && N_j == 0) {
      N_j += (isCheckerBoarded ? 2 : 1) * lowModes.evec.size();
    } else {
      N_i += (isCheckerBoarded ? 2 : 1) * lowModes.evec.size();
      N_j += (isCheckerBoarded ? 2 : 1) * lowModes.evec.size();
    }
  }
  int block = par().block;

  LOG(Message) << "Computing all-to-all one-link J.A meson fields" << std::endl;
  if (hasLowModes)
    LOG(Message) << "Low Modes: '" << par().lowModes << "'" << std::endl;
  if (!par().left.empty())
    LOG(Message) << "Left: '" << par().left << "'" << std::endl;
  if (!par().right.empty())
    LOG(Message) << "Right: '" << par().right << "'" << std::endl;
  LOG(Message) << "Momenta:" << std::endl;
  for (auto &p : _mom) {
    LOG(Message) << "  " << p << std::endl;
  }
  LOG(Message) << "Samples: " << par().nEmFields
               << " x " << env().getNd() << " J_kA_k channels" << std::endl;
  LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j
               << " (filesize "
               << sizeString(nt * N_i * N_j * sizeof(HADRONS_A2AM_IO_TYPE))
               << "/momentum/channel)" << std::endl;

  auto &ph = envGet(std::vector<ComplexField>, _momphName);
  auto &Amu = envGet(std::vector<ComplexField>, _EmName);

  startTimer("Momentum phases");
  for (unsigned int j = 0; j < ph.size(); ++j) {
    Complex i(0.0, 1.0);

    envGetTmp(ComplexField, coor);
    ph[j] = Zero();
    for (unsigned int mu = 0; mu < _mom[j].size(); mu++) {
      LatticeCoordinate(coor, mu);
      ph[j] = ph[j] + (_mom[j][mu] / env().getDim(mu)) * coor;
    }
    ph[j] = exp((Real)(2 * M_PI) * i * ph[j]);
  }
  stopTimer("Momentum phases");

  auto &photon = envGet(PhotonR, par().EmFunc);
  auto &w = envGet(EmComp, "_" + par().EmFunc + "_weight");
  auto &Ufield = envGet(LatticeGaugeField, par().gauge);

  envGetTmp(Computation, computation);

  int orthogDir = env().getNd() - 1;

  // Per-sample loop: one worker + one Computation pass per stochastic
  // sample (the J.A one-link ctor takes exactly one sample's Nd raw
  // fields). The rng is seeded ONCE before the loop: sample j draws the
  // (j+1)-th StochasticField, the same sequence the local A_slash module
  // consumes for its sample j (cross-module reproducibility).
  auto &rng = rng4d();
  if (!par().EmSeedString.empty())
    rng.SeedUniqueString(par().EmSeedString);

  envGetTmp(EmField, AField);

  for (int j = 0; j < par().nEmFields; ++j) {
    startTimer("Stochastic photon sample");
    photon.StochasticField(AField, rng, w);
    // RAW components: NO applyCoeffsAndPhase fold, NO per-direction sum
    // (the task dresses W~_k = A_k . U_k and bakes eta_k . (1/2) itself
    // -- folding here would double the coefficient).
    for (unsigned int k = 0; k < env().getNd(); ++k) {
      Amu[k] = PeekIndex<LorentzIndex>(AField, k);
    }
    stopTimer("Stochastic photon sample");

    // per-sample IO lambdas: (k, m) = (component, momentum); filename
    // EmSeedString_j_k_momcomponents; amu_seed sample-qualified.
    auto AmuIOnameFn = [this, j](const unsigned int k,
                                 const unsigned int m) {
      std::stringstream ss;

      ss << par().EmSeedString << "_" << j << "_" << int(k) << "_";

      for (unsigned int mu = 0; mu < _mom[m].size(); ++mu) {
        ss << _mom[m][mu] << ((mu == _mom[m].size() - 1) ? "" : "_");
      }

      return ss.str();
    };

    auto AmuFilenameFn = [this, &AmuIOnameFn](const unsigned int k,
                                              const unsigned int m) {
      return par().output + "." + std::to_string(vm().getTrajectory()) + "/" +
             AmuIOnameFn(k, m) + ".h5";
    };

    auto AmuMetadataFn = [this, j](const unsigned int k,
                                   const unsigned int m) {
      ASlashOneLinkMesonFieldMILCMetadata md;

      for (auto pmu : _mom[m]) {
        md.momentum.push_back(pmu);
      }

      md.amu_seed = par().EmSeedString.empty()
                        ? "<blank> see output file for seed string!"
                        : par().EmSeedString + "_" + std::to_string(j);
      md.amu_index = k;

      return md;
    };

    Kernel kernel(envGetGrid(FermionField), ph, Amu, &Ufield, orthogDir);

    if (hasLowModes) {
      auto &lowModes = envGet(Pack, par().lowModes);
      if (isCheckerBoarded) {
        auto &pairLeft = envGet(A2ALowModePairSourceMILC<FermionField>,
                                par().cbPairsLeft);
        auto &pairRight = envGet(A2ALowModePairSourceMILC<FermionField>,
                                 par().cbPairsRight);
        computation.execute(*left, *right, kernel, AmuIOnameFn, AmuFilenameFn,
                            AmuMetadataFn, &lowModes.evec, lowModes.eval,
                            nullptr, &pairLeft, &pairRight);
      } else {
        computation.execute(*left, *right, kernel, AmuIOnameFn, AmuFilenameFn,
                            AmuMetadataFn, &lowModes.evec, lowModes.eval);
      }
    } else {
      computation.execute(*left, *right, kernel, AmuIOnameFn, AmuFilenameFn,
                          AmuMetadataFn);
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MContraction_A2AASlashOneLinkMesonField_hpp_
