/*
 * A2AVectorsMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
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
#ifndef A2AVectors_hpp_
#define A2AVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Class to generate V & W all-to-all vectors                 *
 ******************************************************************************/
template <typename FImpl> class A2AVectorsMILC {
public:
  FERM_TYPE_ALIASES(FImpl, );
  SOLVER_TYPE_ALIASES(FImpl, );
  HADRONS_DEFINE_SCHUR_OP(SchurOp, FImpl);

public:
  A2AVectorsMILC(FMat &action, Solver &solver);
  virtual ~A2AVectorsMILC(void) = default;

  void removeLowModeProjection(std::vector<FermionField> &wout,
                               const std::vector<FermionField> &evecs,
                               const std::vector<ComplexD> evals);
  inline void makeLowModeCBeooe(FermionField &out, const FermionField &evec,
                                const Complex eval);
  void makeLowModePairs(typename std::vector<FermionField>::iterator vecOut,
                        const typename std::vector<FermionField>::iterator evec,
                        const Complex eval);
  void makeLowModePairs(typename std::vector<FermionField>::iterator vout,
                        typename std::vector<FermionField>::iterator wout,
                        const typename std::vector<FermionField>::iterator evec,
                        const Complex eval);
  void makeLowModePairs(typename std::vector<FermionField>::iterator vecOut,
                        typename std::vector<Complex>::iterator evalOut,
                        const typename std::vector<FermionField>::iterator evec,
                        const Complex eval);
  void
  makeLowModePairs5D(typename std::vector<FermionField>::iterator vout,
                     typename std::vector<FermionField>::iterator vout5,
                     typename std::vector<FermionField>::iterator wout,
                     typename std::vector<FermionField>::iterator wout5,
                     const typename std::vector<FermionField>::iterator evec,
                     const Complex eval);

  void makeHighModeV(FermionField &vout, const FermionField &noise);
  void makeHighModeV(FermionField &vout, const FermionField &noise,
                     const FermionField &guess);

protected:
  FMat &action_;
  Solver &solver_;
  GridBase *frbGrid_, *gGrid_, *fGrid_;
  FermionField sol_rb_, temp_, temp5_;

  SchurOp<FMat, FermionField> op_;
};

/******************************************************************************
 *            A2AVectorsMILC template implementation                  *
 ******************************************************************************/
template <typename FImpl>
A2AVectorsMILC<FImpl>::A2AVectorsMILC(FMat &action, Solver &solver)
    : action_(action), solver_(solver), fGrid_(action_.FermionGrid()),
      frbGrid_(action_.FermionRedBlackGrid()), gGrid_(action_.GaugeGrid()),
      sol_rb_(frbGrid_), temp_(frbGrid_), temp5_(fGrid_), op_(action_) {}

template <typename FImpl>
inline void A2AVectorsMILC<FImpl>::makeLowModeCBeooe(FermionField &out,
                                                     const FermionField &evec,
                                                     const Complex eval) {
  int cb = evec.Checkerboard();
  int cbNeg = (cb == Even) ? Odd : Even;

  temp_ = Zero();
  out.Checkerboard() = cbNeg;
  temp_.Checkerboard() = cbNeg;
  action_.Meooe(evec, temp_);
  out = (1.0 / eval) * temp_;
}

template <typename FImpl>
void A2AVectorsMILC<FImpl>::removeLowModeProjection(
    std::vector<FermionField> &wout, const std::vector<FermionField> &evecs,
    const std::vector<ComplexD> evals) {
  // Alternate form
  int cb = evecs[0].Checkerboard();
  int cbNeg = (cb == Even) ? Odd : Even;

  FermionField rbw(frbGrid_), rbwNeg(frbGrid_);

  rbw = Zero();
  rbwNeg = Zero();

  // Normalize vectors so that checkerboard has magnitude 1/sqrt(2)
  RealD norm = 1 / ::sqrt(norm2(evecs[0]));

  for (auto &w : wout) {
    rbw.Checkerboard() = cb;
    rbwNeg.Checkerboard() = cbNeg;

    pickCheckerboard(cb, rbw, w);
    pickCheckerboard(cbNeg, rbwNeg, w);

    // Add up W vector projection onto provided evec checkerboard
    temp_ = Zero();
    temp_.Checkerboard() = cb;
    for (int i = 0; i < evecs.size(); i++) {
      const FermionField &e = evecs[i];
      axpy(temp_, TensorRemove(innerProduct(e, rbw)), e, temp_);
    }

    // Subtract projected component from original. (factor of 2 compensates for
    // normalization of checkerboard to 1/2)
    axpy(rbw, -norm, temp_, rbw);
    setCheckerboard(w, rbw);

    action_.MeooeDag(rbwNeg, rbw); // Move cbNeg component of W to cb

    // Add up cbNeg checkerboard of W vector projection
    temp_ = Zero();
    temp_.Checkerboard() = cb;
    for (int i = evecs.size() - 1; i >= 0; i--) {
      RealD eval_Dinv =
          1.0 / pow(evals[i].imag(),
                    2); // using Meooe twice brings two factors of 1/eval_D
      const FermionField &e = evecs[i];
      axpy(temp_, eval_Dinv * TensorRemove(innerProduct(e, rbw)), e, temp_);
    }
    rbw.Checkerboard() = cbNeg;
    action_.Meooe(temp_, rbw); // Move projection back to cbNeg checkerboard
    axpy(rbwNeg, -norm, rbw,
         rbwNeg); // Subtract projected component from original.
    setCheckerboard(w, rbwNeg);
  }

  // int cb = evecs[0].Checkerboard();
  // int cbNeg = (cb==Even) ? Odd : Even;

  // FermionField Mevec(fGrid_), Mdagevec(fGrid_), evecNeg(frbGrid_);

  // evecNeg.Checkerboard() = cbNeg;

  // for (int i=0;i<evecs.size();i++) {
  //     ComplexD eval_D = ComplexD(0.0,evals[i].imag());
  //     makeLowModeCBeooe(evecNeg,evecs[i],eval_D);

  //     setCheckerboard(Mevec,evecNeg);
  //     setCheckerboard(Mevec,evecs[i]);

  //     if (cb == Even) {
  //         temp_.Checkerboard() = cbNeg;
  //         temp_ = -evecNeg;
  //         setCheckerboard(Mdagevec,temp_);
  //         setCheckerboard(Mdagevec,evecs[i]);
  //     } else {
  //         temp_.Checkerboard() = cb;
  //         temp_ = -evecs[i];
  //         setCheckerboard(Mdagevec,temp_);
  //         setCheckerboard(Mdagevec,evecNeg);
  //     }

  //     std::cout << GridLogMessage <<  norm2(Mevec) << norm2(Mdagevec) <<
  //     std::endl; for (auto &w:wout) {
  //         auto ip = innerProduct(Mevec,w)/norm2(Mevec);
  //         w = w - ip*Mevec;
  //         ip = innerProduct(Mdagevec,w)/norm2(Mdagevec);
  //         w = w - ip*Mdagevec;
  //     }
  // }
}

template <typename FImpl>
void A2AVectorsMILC<FImpl>::makeLowModePairs(
    typename std::vector<FermionField>::iterator vecOut,
    const typename std::vector<FermionField>::iterator evec,
    const Complex eval) {
  double cbEven = (*evec).Checkerboard() == Even;
  int cbParity = cbEven ? Even : Odd;
  int cbParityNeg = !cbEven ? Even : Odd;

  // Expects eigenvalues of M
  ComplexD eval_D = ComplexD(0.0, eval.imag());

  makeLowModeCBeooe(sol_rb_, evec, eval_D);

  setCheckerboard(*vecOut, sol_rb_);
  setCheckerboard(*vecOut, *evec);

  if (cbEven) {
    pickCheckerboard(cbParityNeg, temp_, *(vecOut + 1));
    temp_ = -1.0 * sol_rb_;

    setCheckerboard(*(vecOut + 1), temp_);
    setCheckerboard(*(vecOut + 1), *evec);
  } else {
    pickCheckerboard(cbParity, temp_, *(vecOut + 1));
    temp_ = -1.0 * (*evec);

    setCheckerboard(*(vecOut + 1), temp_);
    setCheckerboard(*(vecOut + 1), sol_rb_);
  }
}

template <typename FImpl>
void A2AVectorsMILC<FImpl>::makeLowModePairs(
    typename std::vector<FermionField>::iterator vecOut,
    typename std::vector<Complex>::iterator evalOut,
    const typename std::vector<FermionField>::iterator evec,
    const Complex eval) {
  makeLowModePairs(vecOut, evec, eval);

  *evalOut = 1.0 / eval;
  *(evalOut + 1) = 1.0 / conjugate(eval);
}

template <typename FImpl>
void A2AVectorsMILC<FImpl>::makeLowModePairs(
    typename std::vector<FermionField>::iterator vout,
    typename std::vector<FermionField>::iterator wout,
    const typename std::vector<FermionField>::iterator evec,
    const Complex eval) {
  std::vector<ComplexD> evals(2);

  makeLowModePairs(wout, evals.begin(), evec, eval);

  *vout = evals[0] * (*wout);
  *(vout + 1) = evals[1] * (*(wout + 1));
}

template <typename FImpl>
void A2AVectorsMILC<FImpl>::makeLowModePairs5D(
    typename std::vector<FermionField>::iterator vout,
    typename std::vector<FermionField>::iterator vout5,
    typename std::vector<FermionField>::iterator wout,
    typename std::vector<FermionField>::iterator wout5,
    const typename std::vector<FermionField>::iterator evec,
    const Complex eval) {
  makeLowModePairs(vout5, wout5, evec, eval);
  action_.ExportPhysicalFermionSolution(*vout5, *vout);
  action_.ExportPhysicalFermionSolution(*(vout5 + 1), *(vout + 1));
  action_.DminusDag(temp5_, *wout5);
  action_.ExportPhysicalFermionSolution(temp5_, *wout);
  action_.DminusDag(temp5_, *(wout5 + 1));
  action_.ExportPhysicalFermionSolution(temp5_, *(wout + 1));
}

template <typename FImpl>
void A2AVectorsMILC<FImpl>::makeHighModeV(FermionField &vout,
                                          const FermionField &noise) {
  solver_(vout, noise);
}

template <typename FImpl>
void A2AVectorsMILC<FImpl>::makeHighModeV(FermionField &vout,
                                          const FermionField &noise,
                                          const FermionField &guess) {
  solver_(vout, noise, guess);
}

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_
