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
// A2AWorkerSpinTasteStencil full-grid worker. Supports full-grid and
// checkerboarded low modes (the latter via two distinct EigenPackCBPairs
// module instances named by cbPairsLeft/cbPairsRight; the kernel packs the
// CB arrays into full-grid objects, contracts with the worker's CB
// ContractType modes, and converts the raw parity partials to the legacy
// interleaved layout). Momentum projection is still rejected at setup (the
// stencil worker does not implement it). The pre-stencil implementation
// lives on as StagA2AMesonFieldLegacy in MesonFieldLegacy.hpp.
class MesonFieldMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFieldMILCPar, int, block, std::string,
                                  lowModes, std::string, left, std::string,
                                  cbPairsLeft, std::string, cbPairsRight,
                                  std::string, right, std::string, output,
                                  SpinTasteParams, spinTaste,
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
    if (left_o == nullptr && right_o == nullptr) {
      // Full-array entry on the full-grid arrays (high modes).
      MesonFunctionStencil<FImpl>(m, left_e, right_e);
    } else {
      MesonFunctionStencilCB<FImpl>(m, left_e, left_o, right_e, right_o);
    }
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
    _sigma.clear();
    StagGamma spinTaste;
    for (auto &g : gammas) {
      spinTaste.setSpinTaste(g);
      int pc = StagGamma::popcountShift(spinTaste._spin, spinTaste._taste);
      _sigma.push_back((pc & 1) ? -1.0 : 1.0);
    }
  }

private:
  template <typename TFImpl, typename... Args>
  IfNotStag<TFImpl, void> MesonFunctionStencil(Args &&...) {
    HADRONS_ERROR(Implementation, "MesonField stencil kernel requires a "
                                  "staggered fermion implementation");
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

  template <typename TFImpl, typename... Args>
  IfNotStag<TFImpl, void> MesonFunctionStencilCB(Args &&...) {
    HADRONS_ERROR(Implementation, "MesonField stencil kernel requires a "
                                  "staggered fermion implementation");
  }

  template <typename TFImpl>
  IfStag<TFImpl, void> MesonFunctionStencilCB(A2AMatrixSet<T> &m,
                                              const FermionField *left_e,
                                              const FermionField *left_o,
                                              const FermionField *right_e,
                                              const FermionField *right_o) {
    int N_ii = m.dimension(3), N_jj = m.dimension(4);
    bool leftCB = (left_o != nullptr), rightCB = (right_o != nullptr);
    if (leftCB && (N_ii % 2 != 0)) {
      HADRONS_ERROR(Size, "MesonField kernel: checkerboarded left block "
                          "size must be even, got " + std::to_string(N_ii));
    }
    if (rightCB && (N_jj % 2 != 0)) {
      HADRONS_ERROR(Size, "MesonField kernel: checkerboarded right block "
                          "size must be even, got " + std::to_string(N_jj));
    }
    int sizeL = leftCB ? N_ii / 2 : N_ii;
    int sizeR = rightCB ? N_jj / 2 : N_jj;

    ContractType ct = (leftCB && rightCB) ? ContractType::BothHalf
                      : leftCB            ? ContractType::LeftHalf
                                          : ContractType::RightHalf;

    const FermionField *lhs, *rhs;
    if (leftCB) {
      growPack(_packL, sizeL, left_e, left_o);
      lhs = _packL.data();
    } else {
      lhs = left_e;
    }
    if (rightCB) {
      growPack(_packR, sizeR, right_e, right_o);
      rhs = _packR.data();
    } else {
      rhs = right_e;
    }

    // The worker emits raw parity partials in one uniform block layout
    // (2*sizeL, sizeR): rows [0,sizeL) = M0 (even source sites), rows
    // [sizeL,2*sizeL) = M1 (odd). That layout equals mBlock only for
    // LeftHalf (RightHalf has equal element count but different dims,
    // BothHalf half the elements) -- always contract into a scratch at
    // the worker-native shape, then convert to the legacy interleaved
    // (M, M-dagger) slot layout the framework bake and HDF5 consumers
    // expect.
    _scratch.resize(m.dimension(0) * m.dimension(1) * m.dimension(2) *
                    (2 * sizeL) * sizeR);
    A2AMatrixSet<T> mScratch(_scratch.data(), m.dimension(0), m.dimension(1),
                             m.dimension(2), 2 * sizeL, sizeR);

    // Contents change at fixed buffer addresses: force the worker to
    // re-pack (its address cache would otherwise contract stale data).
    _stencilWorker->resetCache();
    _stencilWorker->StagMesonField(mScratch, lhs, rhs, sizeL, sizeR, ct);

    reconstructLegacy(m, mScratch, sizeL, sizeR, leftCB, rightCB);
  }

  // Pack n CB pairs into cached full-grid objects: E copy on even
  // sites, O copy on odd (the setCheckerboard convention the worker's CB
  // modes expect; Test_a2a_stencil.cc packing precedent).
  void growPack(std::vector<FermionField> &pack, const int n,
                const FermionField *even, const FermionField *odd) {
    if ((int)pack.size() < n) {
      pack.resize(n, _stencilWorker->_grid);
    }
    for (int k = 0; k < n; ++k) {
      pack[k] = Zero();
      setCheckerboard(pack[k], even[k]);
      setCheckerboard(pack[k], odd[k]);
    }
  }

  // Convert the raw (2*sizeL, sizeR) parity partials M0/M1 into the legacy
  // interleaved slot layout -- the exact algebraic image of the legacy
  // simdSumHalf/simdSumMixed tables, validated value-for-value against the
  // legacy 4-arg worker in Grid's Test_a2a_stencil.cc reconstructLegacySlot:
  // sigma = -1 for odd popcount(spin^taste), entering RightHalf's rc=1 slot
  // and BothHalf's rc=1 column slots only (LeftHalf is sigma-free).
  void reconstructLegacy(A2AMatrixSet<T> &m, const A2AMatrixSet<T> &mScratch,
                         const int sizeL, const int sizeR, const bool leftCB,
                         const bool rightCB) {
    int next = m.dimension(0), nstr = m.dimension(1), nt = m.dimension(2);
    for (int e = 0; e < next; ++e)
      for (int s = 0; s < nstr; ++s) {
        RealD sig = _sigma[s];
        for (int t = 0; t < nt; ++t) {
          for (int l = 0; l < sizeL; ++l)
            for (int r = 0; r < sizeR; ++r) {
              T M0 = mScratch(e, s, t, l, r);
              T M1 = mScratch(e, s, t, sizeL + l, r);
              if (leftCB && rightCB) {
                m(e, s, t, 2 * l, 2 * r) = M0 + M1;
                m(e, s, t, 2 * l, 2 * r + 1) = sig * (M0 - M1);
                m(e, s, t, 2 * l + 1, 2 * r) = M0 - M1;
                m(e, s, t, 2 * l + 1, 2 * r + 1) = sig * (M0 + M1);
              } else if (leftCB) {
                m(e, s, t, 2 * l, r) = M0 + M1;
                m(e, s, t, 2 * l + 1, r) = M0 - M1;
              } else {
                m(e, s, t, l, 2 * r) = M0 + M1;
                m(e, s, t, l, 2 * r + 1) = sig * (M0 - M1);
              }
            }
        }
      }
  }

private:
  double _vol;
  std::unique_ptr<A2AWorkerSpinTasteStencil<FImpl>> _stencilWorker;
  std::vector<FermionField> _packL, _packR; // full-grid packing buffers (grow-only)
  std::vector<T> _scratch;                  // raw (2*sizeL, sizeR) worker output
  std::vector<RealD> _sigma;                // per-gamma sign (+1/-1)
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
    if (!par().cbPairsLeft.empty()) {
      in.push_back(par().cbPairsLeft);
      in.push_back(par().cbPairsRight);
    }
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

  if (!par().cbPairsLeft.empty() != !par().cbPairsRight.empty()) {
    HADRONS_ERROR(Argument,
                  "MesonField: 'cbPairsLeft' and 'cbPairsRight' must be "
                  "set together (two EigenPackCBPairs module instances)");
  }
  if (!par().cbPairsLeft.empty() &&
      par().cbPairsLeft == par().cbPairsRight) {
    HADRONS_ERROR(Argument, "MesonField: 'cbPairsLeft' and 'cbPairsRight' "
                            "must name distinct EigenPackCBPairs instances");
  }
  if (!par().cbPairsLeft.empty() && par().lowModes.empty()) {
    HADRONS_ERROR(Argument, "MesonField: 'cbPairsLeft'/'cbPairsRight' "
                            "require 'lowModes'");
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
      if (isCheckerBoarded) {
        auto &pairLeft = envGet(A2ALowModePairSourceMILC<FermionField>,
                                par().cbPairsLeft);
        auto &pairRight = envGet(A2ALowModePairSourceMILC<FermionField>,
                                 par().cbPairsRight);
        computationStencil.execute(*left, *right, kernel, gammaIOnameFn,
                                   gammaFilenameFn, gammaMetadataFn,
                                   &lowModes.evec, lowModes.eval, nullptr,
                                   &pairLeft, &pairRight);
      } else {
        computationStencil.execute(*left, *right, kernel, gammaIOnameFn,
                                   gammaFilenameFn, gammaMetadataFn,
                                   &lowModes.evec, lowModes.eval);
      }
    } else {
      computationStencil.execute(*left, *right, kernel, gammaIOnameFn,
                                 gammaFilenameFn, gammaMetadataFn);
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonFieldMILC_hpp_
