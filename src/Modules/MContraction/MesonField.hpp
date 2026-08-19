/*
 * MesonFieldMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
 * Author: paboyle <paboyle@ph.ed.ac.uk>
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
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef HadronsMILC_MContraction_A2AMesonField_hpp_
#define HadronsMILC_MContraction_A2AMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <A2AVectors.hpp>
#include <EigenPack.hpp>
#include <A2AMatrix.hpp>
#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

// Stencil flavour of the A2A meson-field module: the
// A2AWorkerSpinTasteStencil full-grid worker. v1: the stencil path
// implements ContractType::Full only (checkerboarded low modes and momentum
// projection are rejected at setup) and the kernel base is still the 4-field
// A2AKernelMILC shared with the legacy module until the A2AMatrix fork. The
// pre-stencil implementation lives on as StagA2AMesonFieldLegacy in
// MesonFieldLegacy.hpp.
class MesonFieldMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFieldMILCPar, int, block, std::string,
                                  lowModes, std::string, left, std::string,
                                  action, std::string, right, std::string,
                                  output, SpinTasteParams, spinTaste,
                                  std::vector<std::string>, mom);
  MesonFieldMILCPar() {}
};

class MesonFieldMILCMetadata : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFieldMILCMetadata, std::vector<RealF>,
                                  momentum, StagGamma::StagAlgebra, gamma_spin,
                                  StagGamma::StagAlgebra, gamma_taste);

  MesonFieldMILCMetadata()
      : momentum{}, gamma_spin(StagGamma::StagAlgebra::undef),
        gamma_taste(StagGamma::StagAlgebra::undef) {}
};

template <typename T, typename FImpl>
class MesonFieldKernelMILC
    : public A2AKernelMILC<T, typename FImpl::FermionField> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  MesonFieldKernelMILC(GridBase *grid) {
    _vol = 1.;
    for (auto &d : grid->GlobalDimensions()) {
      _vol *= d;
    }
  }
  virtual ~MesonFieldKernelMILC(void) {};

  virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left_e,
                          const FermionField *left_o,
                          const FermionField *right_e,
                          const FermionField *right_o) {
    // Full-array entry on the full-grid arrays; CB blocks are rejected at
    // module setup (v1: stencil implements ContractType::Full only).
    MesonFunctionStencil<FImpl>(m, left_e, right_e);
  }

  virtual double flops(const unsigned int blockSizei,
                       const unsigned int blockSizej, int cbDiv = 1) {
    return _vol / cbDiv * _stencilWorker->getFlops() * blockSizei *
           blockSizej;
  }

  virtual double bytes(const unsigned int blockSizei,
                       const unsigned int blockSizej) {
    return -1.0;
  }

  virtual double kernelTime() { return _stencilWorker->_t_kernel; }
  virtual double globalSumTime() { return _stencilWorker->_t_gsum; }
  void setWorkerStencil(GridCartesian *grid,
                        const std::vector<ComplexField> &mom,
                        const std::vector<StagGamma::SpinTastePair> &gammas,
                        int orthogDir, LatticeGaugeField *U) {
    _stencilWorker = std::make_unique<A2AWorkerSpinTasteStencil<FImpl>>(
        grid, mom, gammas, U, orthogDir);
  }

private:
  template <typename TFImpl, typename... Args>
  IfNotStag<TFImpl, void> MesonFunctionStencil(Args &&...) {
    assert(0);
  }

  template <typename TFImpl>
  IfStag<TFImpl, void> MesonFunctionStencil(A2AMatrixSet<T> &m,
                                            const FermionField *left_e,
                                            const FermionField *right_e) {
    // Full-array entry on the full-grid arrays; CB blocks are rejected at
    // module setup.
    _stencilWorker->StagMesonField(m, left_e, right_e, (int)m.dimension(3),
                                   (int)m.dimension(4));
  }

private:
  double _vol;
  std::unique_ptr<A2AWorkerSpinTasteStencil<FImpl>> _stencilWorker;
};

template <typename FImpl, typename Pack>
class TMesonFieldMILC : public Module<MesonFieldMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );
  typedef A2AMatrixBlockComputationMILC<
      Complex, FermionField, MesonFieldMILCMetadata, HADRONS_A2AM_IO_TYPE>
      Computation;
  typedef MesonFieldKernelMILC<Complex, FImpl> Kernel;

public:
  // constructor
  TMesonFieldMILC(const std::string name);
  // destructor
  virtual ~TMesonFieldMILC(void) {};

  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);

private:
  std::string _momphName;
  std::vector<StagGamma::SpinTastePair> _gammas;
  std::vector<std::vector<Real>> _mom;
};

MODULE_REGISTER(StagA2AMesonField,
                ARG(TMesonFieldMILC<STAGIMPL, MassShiftEigenPack<STAGIMPL>>),
                MContraction);

/******************************************************************************
 *                  TMesonFieldMILC implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TMesonFieldMILC<FImpl, Pack>::TMesonFieldMILC(const std::string name)
    : Module<MesonFieldMILCPar>(name), _momphName(name + "_momph") {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TMesonFieldMILC<FImpl, Pack>::getInput(void) {
  std::vector<std::string> in = {};
  if (!par().left.empty())
    in.push_back(par().left);
  if (!par().right.empty())
    in.push_back(par().right);

  if (!par().lowModes.empty()) {
    if (!par().action.empty())
      in.push_back(par().action);
    in.push_back(par().lowModes);
  }

  if (!par().spinTaste.gauge.empty()) {
    in.push_back(par().spinTaste.gauge);
  }

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TMesonFieldMILC<FImpl, Pack>::getOutput(void) {
  std::vector<std::string> out = {};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TMesonFieldMILC<FImpl, Pack>::setup(void) {
  _gammas = StagGamma::ParseSpinTasteString(par().spinTaste.gammas,
                                            par().spinTaste.applyG5);

  if (_gammas.empty()) {
    LOG(Warning) << "MesonField: empty spin-taste gamma list; no meson "
                    "fields will be computed"
                 << std::endl;
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
  bool allzero = true;
  if (par().mom.size() == 1) {
    for (auto p : _mom[0]) {
      if (p != 0)
        allzero = false;
    }
  }
  if (allzero)
    nmom = 0;

  if (!par().action.empty()) {
    HADRONS_ERROR(Implementation,
                  "MesonField: checkerboarded low modes are not supported "
                  "yet (stencil implements ContractType::Full only); unset "
                  "'action' or use StagA2AMesonFieldLegacy");
  }
  bool anyMomentum = false;
  for (auto &p : _mom)
    for (auto pmu : p)
      if (pmu != 0.0)
        anyMomentum = true;
  if (anyMomentum) {
    HADRONS_ERROR(Argument,
                  "MesonField: momentum projection is not implemented on "
                  "the stencil path; use zero momentum or "
                  "StagA2AMesonFieldLegacy");
  }
  if (_mom.size() > 1) {
    HADRONS_ERROR(Argument, "MesonField: at most one (zero) momentum is "
                            "supported; got " +
                                std::to_string(_mom.size()));
  }
  if (par().spinTaste.gauge.empty()) {
    HADRONS_ERROR(Argument, "MesonField requires 'spinTaste.gauge'");
  }

  envCache(std::vector<ComplexField>, _momphName, 1, nmom,
           envGetGrid(ComplexField));

  envTmpLat(ComplexField, "coor");

  envTmp(Computation, "computationStencil", 1, envGetGrid(FermionField),
         env().getNd() - 1, _mom.size(), _gammas.size(), par().block, this);

  envTmp(std::vector<FermionField>, "dummy", 1, 0, envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TMesonFieldMILC<FImpl, Pack>::execute(void) {
  bool hasLowModes = (!par().lowModes.empty());
  bool isCheckerBoarded = (!par().action.empty());

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

  LOG(Message) << "Computing all-to-all meson fields" << std::endl;
  if (hasLowModes)
    LOG(Message) << "Low Modes: '" << par().lowModes << "'" << std::endl;

  if (!(par().left.empty() && par().right.empty())) {
    if (!par().left.empty())
      LOG(Message) << "Left: '" << par().left << "'" << std::endl;
    if (!par().right.empty())
      LOG(Message) << "Right: '" << par().right << "'" << std::endl;
  }

  LOG(Message) << "Momenta:" << std::endl;

  for (auto &p : _mom) {
    LOG(Message) << "  " << p << std::endl;
  }

  LOG(Message) << "Spin bilinears:" << std::endl;

  for (auto &g : _gammas) {
    LOG(Message) << "  " << StagGamma::GetName(g) << std::endl;
  }

  LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j
               << " (filesize "
               << sizeString(nt * N_i * N_j * sizeof(HADRONS_A2AM_IO_TYPE))
               << "/momentum/bilinear)" << std::endl;

  auto &ph = envGet(std::vector<ComplexField>, _momphName);

  startTimer("Momentum phases");
  for (unsigned int j = 0; j < ph.size(); ++j) {
    Complex i(0.0, 1.0);
    std::vector<Real> p;

    envGetTmp(ComplexField, coor);
    ph[j] = Zero();
    for (unsigned int mu = 0; mu < _mom[j].size(); mu++) {
      LatticeCoordinate(coor, mu);
      ph[j] = ph[j] + (_mom[j][mu] / env().getDim(mu)) * coor;
    }
    ph[j] = exp((Real)(2 * M_PI) * i * ph[j]);
  }
  stopTimer("Momentum phases");

  auto gammaIOnameFn = [this](const unsigned int m, const unsigned int g) {
    std::stringstream ss;

    ss << StagGamma::GetName(_gammas[g]) << "_";

    for (unsigned int mu = 0; mu < _mom[m].size(); ++mu) {
      ss << _mom[m][mu] << ((mu == _mom[m].size() - 1) ? "" : "_");
    }

    return ss.str();
  };

  auto gammaFilenameFn = [this, &gammaIOnameFn](const unsigned int m,
                                                const unsigned int g) {
    return par().output + "." + std::to_string(vm().getTrajectory()) + "/" +
           gammaIOnameFn(m, g) + ".h5";
  };

  auto gammaMetadataFn = [this](const unsigned int m, const unsigned int g) {
    MesonFieldMILCMetadata md;

    for (auto pmu : _mom[m]) {
      md.momentum.push_back(pmu);
    }

    md.gamma_spin = _gammas[g].first;
    md.gamma_taste = _gammas[g].second;

    return md;
  };

  envGetTmp(Computation, computationStencil);

  Kernel kernel(envGetGrid(FermionField));

  GaugeField *U = nullptr;
  if (!par().spinTaste.gauge.empty()) {
    U = env().template getObject<GaugeField>(par().spinTaste.gauge);
  }

  int orthogDir = env().getNd() - 1;

  // One stencil worker serves every gamma: the task handles mixed
  // popcounts natively (per-gamma endpoint tables), amortizing the
  // gauge-chain setup across the whole run.
  if (_gammas.size() > 0) {
    GridCartesian *grid =
        dynamic_cast<GridCartesian *>(envGetGrid(FermionField));
    if (grid == nullptr) {
      HADRONS_ERROR(Implementation, "MesonField requires a Cartesian grid");
    }
    kernel.setWorkerStencil(grid, ph, _gammas, orthogDir, U);
    if (hasLowModes) {
      auto &lowModes = envGet(Pack, par().lowModes);
      computationStencil.execute(*left, *right, kernel, gammaIOnameFn,
                                 gammaFilenameFn, gammaMetadataFn,
                                 &lowModes.evec, lowModes.eval);
    } else {
      computationStencil.execute(*left, *right, kernel, gammaIOnameFn,
                                 gammaFilenameFn, gammaMetadataFn);
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonFieldMILC_hpp_
