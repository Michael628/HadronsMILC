/*
 * A2AStencilKernel.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_A2AStencilKernel_hpp_
#define HadronsMILC_A2AStencilKernel_hpp_

#include <Hadrons/Global.hpp>
#include <A2AMatrix.hpp>
#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *        Shared stencil A2A kernel base (A2AWorkerSpinTasteStencil)          *
 ******************************************************************************/
// Full/CB dispatch, growPack, and the per-side content-stamp pack cache
// shared by every kernel that drives A2AWorkerSpinTasteStencil (the
// MesonField and ASlash module kernels). The worker is constructed via
// one of the three setWorker overloads, mirroring the stencil worker's
// ctors one-to-one (parameter order identical to the worker's):
//   setWorker(grid, mom, gammas, U, orthogDir)  -- momentum-aware
//   setWorker(grid, mom, Amu, orthogDir)        -- J.A local
//   setWorker(grid, mom, Amu, U, orthogDir)     -- J.A one-link
template <typename T, typename FImpl>
class A2AStencilKernelMILC
    : public A2AKernelMILC<T, typename FImpl::FermionField> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  A2AStencilKernelMILC(GridBase *grid) {
    _vol = 1.;
    for (auto &d : grid->GlobalDimensions()) {
      _vol *= d;
    }
  }
  virtual ~A2AStencilKernelMILC(void) {};

  virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left_e,
                          const FermionField *left_o,
                          const FermionField *right_e,
                          const FermionField *right_o) {
    if (left_o == nullptr && right_o == nullptr) {
      // Full-array entry on the full-grid arrays (high modes).
      contractStencil<FImpl>(m, left_e, right_e);
    } else {
      contractStencilCB<FImpl>(m, left_e, left_o, right_e, right_o,
                               A2ASideStamp{}, A2ASideStamp{});
    }
  }

  // Stamped dispatcher override -- Full keeps its pointer-based fast
  // path (raw staging pointers advance per block; the worker's own
  // address cache is a correct validity test there); CB blocks take the
  // content-stamped path.
  virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left_e,
                          const FermionField *left_o,
                          const FermionField *right_e,
                          const FermionField *right_o,
                          const A2ASideStamp &lstamp,
                          const A2ASideStamp &rstamp) {
    if (left_o == nullptr && right_o == nullptr) {
      contractStencil<FImpl>(m, left_e, right_e);
    } else {
      contractStencilCB<FImpl>(m, left_e, left_o, right_e, right_o,
                               lstamp, rstamp);
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

  // Worker construction: the three A2AWorkerSpinTasteStencil ctors.
  // Momentum-aware: mom holds one exp(i p.x) phase field per momentum;
  // the output mu axis becomes nSpec*nMom, g-major (mu = spec*nMom + p).
  void setWorker(GridCartesian *grid,
                 const std::vector<ComplexField> &mom,
                 const std::vector<StagGamma> &gammas,
                 LatticeGaugeField *U, int orthogDir) {
    _stencilWorker = std::make_unique<A2AWorkerSpinTasteStencil<FImpl>>(
        grid, mom, gammas, U, orthogDir);
    // (the per-gamma sigma sign table moved into the task: built from
    // the gammas the task already holds, applied at the fused CB store
    // sites in assembleMat)
  }
  // J.A local: the caller's Amu fields ARE the per-channel site phase
  // (any coefficient convention already folded caller-side -- the task
  // adds nothing).
  void setWorker(GridCartesian *grid,
                 const std::vector<ComplexField> &mom,
                 const std::vector<ComplexField> &Amu, int orthogDir) {
    _stencilWorker = std::make_unique<A2AWorkerSpinTasteStencil<FImpl>>(
        grid, mom, Amu, orthogDir);
  }
  // J.A one-link: exactly Nd (4) raw photon fields dressing the gauge
  // links (the task builds W~_k = A_k . U_k and bakes the standard
  // eta_k . (1/2) phase itself -- feed RAW components, never
  // applyCoeffsAndPhase-folded ones).
  void setWorker(GridCartesian *grid,
                 const std::vector<ComplexField> &mom,
                 const std::vector<ComplexField> &Amu,
                 LatticeGaugeField *U, int orthogDir) {
    _stencilWorker = std::make_unique<A2AWorkerSpinTasteStencil<FImpl>>(
        grid, mom, Amu, U, orthogDir);
  }

  // Pack n CB pairs into cached full-grid objects: E copy on even
  // sites, O copy on odd (the setCheckerboard convention the worker's CB
  // modes expect; Test_a2a_stencil.cc packing precedent). ONE fused
  // accelerator pass per pair replaces the former Zero + two
  // setCheckerboard passes (three kernels): every site is written (both
  // parities covered), so no zero-init is needed, and the copies are
  // value-identical to the unfused sequence. The full-grid -> rb-grid
  // site mapping replicates Grid's acceleratorSetCheckerboard
  // (Lattice_transfer.h): coordinate from _rdimensions, parity from
  // _checker_dim_mask, rb index from _ostride with the checker dim
  // halved. Public because nvcc forbids extended __host__ __device__
  // lambdas (accelerator_for) inside private/protected member functions.
  void growPack(std::vector<FermionField> &pack, const int n,
                const FermionField *even, const FermionField *odd) {
    if ((int)pack.size() < n) {
      pack.resize(n, _stencilWorker->_grid);
    }
    for (int k = 0; k < n; ++k) {
      autoView(packW, pack[k], AcceleratorWrite);
      autoView(evenR, even[k], AcceleratorRead);
      autoView(oddR, odd[k], AcceleratorRead);
      const GridBase *rbGrid = even[k].Grid();
      const Coordinate rdimFull = pack[k].Grid()->_rdimensions;
      const Coordinate rdimHalf = rbGrid->_rdimensions;
      const Coordinate cbMask = rbGrid->_checker_dim_mask;
      const Coordinate ostride = rbGrid->_ostride;
      const int ndim = rbGrid->_ndimension;
      const int cbEven = even[k].Checkerboard();
      accelerator_for(ss, pack[k].Grid()->oSites(),
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
        if ((linear & 0x1) == cbEven) {
          coalescedWrite(packW[ss], coalescedRead(evenR[ssh]));
        } else {
          coalescedWrite(packW[ss], coalescedRead(oddR[ssh]));
        }
      });
    }
  }

protected:
  template <typename TFImpl, typename... Args>
  IfNotStag<TFImpl, void> contractStencil(Args &&...) {
    HADRONS_ERROR(Implementation, "stencil A2A kernels require a "
                                  "staggered fermion implementation");
  }

  template <typename TFImpl>
  IfStag<TFImpl, void> contractStencil(A2AMatrixSet<T> &m,
                                       const FermionField *left_e,
                                       const FermionField *right_e) {
    // Full-array entry on the full-grid arrays; CB blocks are rejected at
    // module setup.
    _stencilWorker->StagMesonField(m, left_e, right_e, (int)m.dimension(3),
                                   (int)m.dimension(4));
  }

  template <typename TFImpl, typename... Args>
  IfNotStag<TFImpl, void> contractStencilCB(Args &&...) {
    HADRONS_ERROR(Implementation, "stencil A2A kernels require a "
                                  "staggered fermion implementation");
  }

  template <typename TFImpl>
  IfStag<TFImpl, void> contractStencilCB(A2AMatrixSet<T> &m,
                                         const FermionField *left_e,
                                         const FermionField *left_o,
                                         const FermionField *right_e,
                                         const FermionField *right_o,
                                         const A2ASideStamp &lstamp,
                                         const A2ASideStamp &rstamp) {
    int N_ii = m.dimension(3), N_jj = m.dimension(4);
    bool leftCB = (left_o != nullptr), rightCB = (right_o != nullptr);
    if (leftCB && (N_ii % 2 != 0)) {
      HADRONS_ERROR(Size, "A2A stencil kernel: checkerboarded left block "
                          "size must be even, got " + std::to_string(N_ii));
    }
    if (rightCB && (N_jj % 2 != 0)) {
      HADRONS_ERROR(Size, "A2A stencil kernel: checkerboarded right block "
                          "size must be even, got " + std::to_string(N_jj));
    }
    int sizeL = leftCB ? N_ii / 2 : N_ii;
    int sizeR = rightCB ? N_jj / 2 : N_jj;

    ContractType ct = (leftCB && rightCB) ? ContractType::BothHalf
                      : leftCB            ? ContractType::LeftHalf
                                          : ContractType::RightHalf;

    const FermionField *lhs, *rhs;
    if (leftCB) {
      if (lstamp != _lstampL) {
        growPack(_packL, sizeL, left_e, left_o);
        _stencilWorker->resetLeft();
        _lstampL = lstamp;
      }
      lhs = _packL.data();
    } else {
      lhs = left_e;
    }
    if (rightCB) {
      if (rstamp != _lstampR) {
        growPack(_packR, sizeR, right_e, right_o);
        _stencilWorker->resetRight();
        _lstampR = rstamp;
      }
      rhs = _packR.data();
    } else {
      rhs = right_e;
    }

    // The task's CB modes write the legacy interleaved (M, M-dagger)
    // slot layout into m (raw (2L, R) parity partials are stored on
    // device and expanded to the interleaved slots inside the worker's
    // finishMatOutput). m's dims are already final -- (N_ii, N_jj)
    // with N_ii = 2*sizeL / N_jj = 2*sizeR on the CB-flagged sides.
    //
    // Staleness contract: a side's pack is reused iff its stamp is
    // unchanged. Stamps advance exactly when that side's staged
    // contents change -- EigenPackCBPairSource::get() refills its
    // grow-only even_/odd_ buffers at stable addresses every call and
    // the SwapFn path refills _lowBuf_x[0] in place per block, so RAW
    // POINTER COMPARISON IS NOT A VALIDITY TEST on this path (the
    // pointer-keyed first draft produced stale packs and a bitwise
    // FAIL, 2026-10-12). The i==j alias hands the right side's buffers
    // over as left inputs under a fresh left stamp (mode-i content):
    // one redundant left pack per sweep, same as pre-stamp behavior.
    // growPack growth that relocates the pack still re-fires the
    // worker's own address guard (an extra set, never stale).
    // Cross-sweep note: i restarts at 0 each j-sweep, so left stamps
    // REPEAT across sweeps and left packs are legitimately reused
    // (staging is deterministic per coordinate); setLeft ~NBlock_i and
    // setRight ~NBlock_j per execute, down from one of each per (i,j)
    // block.
    _stencilWorker->StagMesonField(m, lhs, rhs, sizeL, sizeR, ct);
  }

protected:
  double _vol;
  std::unique_ptr<A2AWorkerSpinTasteStencil<FImpl>> _stencilWorker;
  std::vector<FermionField> _packL, _packR; // full-grid packing buffers (grow-only)
  A2ASideStamp _lstampL{}, _lstampR{}; // last-packed side stamps; the
                                       // block=-1 default never matches
};

END_HADRONS_NAMESPACE

#endif // HadronsMILC_A2AStencilKernel_hpp_
