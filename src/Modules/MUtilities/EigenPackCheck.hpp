/*
 * EigenPackCheck.hpp, part of HadronsMILC
 *
 * Copyright (C) 2026
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
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef HadronsMILC_MUtilities_EigenPackCheck_hpp_
#define HadronsMILC_MUtilities_EigenPackCheck_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <EigenPack.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *        Consistency check of a checkerboarded staggered eigenpack           *
 ******************************************************************************
 * Checks the pack the way StagLMA, StagLMAFull and A2AVectorsMILC use it. The
 * pack holds one checkerboard e_k of each eigenvector (the Lanczos output) and
 * the Dirac eigenvalue 2m + i lambda_k. The other checkerboard is rebuilt as
 * w_k = Meooe e_k / lambda_k, which is the other half of an eigenvector of
 * D = 2m + Dhop only if e_k is an eigenvector of Meooe^dag Meooe with
 * eigenvalue exactly lambda_k^2. Everything is measured against the action's
 * Meooe directly, so no operator convention is assumed. Per mode k:
 *
 *   lambdaRayleigh      sqrt(||Meooe e||^2 / ||e||^2), the lambda the vector
 *                       actually carries
 *   normMismatch        ||w||^2 / ||e||^2 - 1 = (lambdaRayleigh / lambda)^2 - 1;
 *                       the two checkerboards of the rebuilt vector have equal
 *                       norm only if this is 0
 *   pairOverlap         <v+|v-> / ||v+||^2 for the +/-lambda partners built as
 *                       in StagLMAFull; they are orthogonal only if this is 0
 *   resStored           ||Meooe^dag Meooe e - lambda^2 e|| / (lambda^2 ||e||)
 *   resRayleigh         the same with lambdaRayleigh^2: the eigenvector's own
 *                       quality, with any eigenvalue error taken out
 *   resFull             ||Dhop v+ - i lambda v+|| / (lambda ||v+||). The rebuilt
 *                       checkerboard satisfies its half of the equation exactly,
 *                       so all of this sits on e's checkerboard
 *   evalShiftOverMass2  (lambda^2 - lambdaRayleigh^2) / (2m)^2. Near 1 if the
 *                       Lanczos operator carried the mass (a Schur operator,
 *                       eigenvalue (2m)^2 + lambda^2) but the pack was built as
 *                       if it were massless; near 0 otherwise
 *   antiHermiticity     ||Meooe w + Meooe^dag w|| / ||Meooe^dag w||, checking
 *                       Dhop^dag = -Dhop, which the reconstruction assumes
 *   gramStored          max normalised |<e_j|e_k>| over the gramBand preceding
 *                       modes
 *   gramRebuilt         the same for the rebuilt checkerboard w
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EigenPackCheckMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(EigenPackCheckMILCPar,
                                  std::string, action,
                                  std::string, eigenPack,
                                  int, nEigs,
                                  unsigned int, gramBand,
                                  std::string, output);
};

class EigenPackCheckMILCResult : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(EigenPackCheckMILCResult,
                                  int, checkerboard,
                                  double, massDiag,
                                  std::vector<double>, lambdaStored,
                                  std::vector<double>, lambdaRayleigh,
                                  std::vector<double>, normStored,
                                  std::vector<double>, normMismatch,
                                  std::vector<double>, pairOverlap,
                                  std::vector<double>, resStored,
                                  std::vector<double>, resRayleigh,
                                  std::vector<double>, resFull,
                                  std::vector<double>, evalShiftOverMass2,
                                  std::vector<double>, antiHermiticity,
                                  std::vector<double>, gramStored,
                                  std::vector<double>, gramRebuilt);
};

// median and max of |x| over the finite entries, with the lambda at the max
inline void eigenPackCheckSummary(const std::string &name,
                                  const std::vector<double> &x,
                                  const std::vector<double> &lambda) {
  std::vector<double> a;
  int kMax = -1;
  double vMax = 0.;

  a.reserve(x.size());
  for (unsigned int k = 0; k < x.size(); ++k) {
    if (!std::isfinite(x[k])) {
      continue;
    }
    a.push_back(std::abs(x[k]));
    if (kMax < 0 || std::abs(x[k]) > vMax) {
      vMax = std::abs(x[k]);
      kMax = k;
    }
  }

  std::ostringstream line;
  line << std::setw(20) << std::left << name;
  if (a.empty()) {
    line << "no finite values";
  } else {
    std::nth_element(a.begin(), a.begin() + a.size() / 2, a.end());
    line << std::scientific << std::setprecision(3)
         << "median |x| " << a[a.size() / 2] << "   max |x| " << vMax
         << " at mode " << kMax << " (lambda " << lambda[kMax] << ")";
  }
  LOG(Message) << line.str() << std::endl;
}

template <typename FImpl, typename Pack>
class TEigenPackCheckMILC : public Module<EigenPackCheckMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TEigenPackCheckMILC(const std::string name);
  // destructor
  virtual ~TEigenPackCheckMILC(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(
    EigenPackCheckMILC,
    ARG(TEigenPackCheckMILC<STAGIMPL, MassShiftEigenPack<STAGIMPL>>),
    MUtilities);

/******************************************************************************
 *                    TEigenPackCheckMILC implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TEigenPackCheckMILC<FImpl, Pack>::TEigenPackCheckMILC(const std::string name)
    : Module<EigenPackCheckMILCPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TEigenPackCheckMILC<FImpl, Pack>::getInput(void) {
  std::vector<std::string> in = {par().action, par().eigenPack};

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TEigenPackCheckMILC<FImpl, Pack>::getOutput(void) {
  std::vector<std::string> out = {};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TEigenPackCheckMILC<FImpl, Pack>::setup(void) {
  if (env().getObjectLs(par().action) > 1) {
    HADRONS_ERROR(Argument, "Ls > 1 not implemented");
  }

  const unsigned int band = par().gramBand;

  envTmp(FermionField, "w", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "y", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "z", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "d", 1, envGetRbGrid(FermionField));
  envTmp(std::vector<FermionField>, "ringE", 1, band,
         envGetRbGrid(FermionField));
  envTmp(std::vector<FermionField>, "ringW", 1, band,
         envGetRbGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TEigenPackCheckMILC<FImpl, Pack>::execute(void) {
  auto &mat = envGet(FMat, par().action);
  auto &epack = envGet(Pack, par().eigenPack);

  envGetTmp(FermionField, w);
  envGetTmp(FermionField, y);
  envGetTmp(FermionField, z);
  envGetTmp(FermionField, d);
  envGetTmp(std::vector<FermionField>, ringE);
  envGetTmp(std::vector<FermionField>, ringW);

  const int nPack = epack.evec.size();
  const int nEigs =
      (par().nEigs < 1 || par().nEigs > nPack) ? nPack : par().nEigs;
  const unsigned int band = par().gramBand;
  const int cb = epack.evec[0].Checkerboard();
  const int cbNeg = (cb == Even) ? Odd : Even;
  const double nan = std::numeric_limits<double>::quiet_NaN();

  LOG(Message) << "Checking " << nEigs << " of " << nPack
               << " eigenvectors of '" << par().eigenPack << "' (stored on the "
               << ((cb == Even) ? "even" : "odd")
               << " checkerboard) against action '" << par().action << "'"
               << std::endl;

  EigenPackCheckMILCResult res;
  std::vector<double> *perMode[] = {
      &res.lambdaStored,       &res.lambdaRayleigh, &res.normStored,
      &res.normMismatch,       &res.pairOverlap,    &res.resStored,
      &res.resRayleigh,        &res.resFull,        &res.evalShiftOverMass2,
      &res.antiHermiticity,    &res.gramStored,     &res.gramRebuilt};
  for (auto v : perMode) {
    v->assign(nEigs, nan);
  }
  res.checkerboard = cb;
  res.massDiag = epack.eval[0].real();

  std::vector<double> ringNe(band, 0.), ringNw(band, 0.);
  // |<a|b>| / sqrt(na nb), through .real()/.imag(): std::abs has no overload
  // for the ComplexD of GPU builds
  auto cosine = [](const FermionField &a, const FermionField &b, double na,
                   double nb) {
    const ComplexD ip = TensorRemove(innerProduct(a, b));
    return std::hypot(ip.real(), ip.imag()) / std::sqrt(na * nb);
  };

  for (int k = 0; k < nEigs; ++k) {
    const FermionField &e = epack.evec[k];
    const double mu = epack.eval[k].real();
    const double lam = epack.eval[k].imag();
    const double lam2 = lam * lam;

    w.Checkerboard() = cbNeg;
    y.Checkerboard() = cb;
    z.Checkerboard() = cb;
    mat.Meooe(e, w);    // the other checkerboard, before the 1/lambda
    mat.MeooeDag(w, y); // Meooe^dag Meooe e
    mat.Meooe(w, z);    // Dhop twice, for the full residual

    const double ne = norm2(e);
    const double nw = norm2(w);
    const double rho = nw / ne; // Rayleigh quotient of Meooe^dag Meooe

    res.lambdaStored[k] = lam;
    res.lambdaRayleigh[k] = std::sqrt(rho);
    res.normStored[k] = ne;

    if (lam > 0.) {
      const double nwScaled = nw / lam2; // ||w / lambda||^2

      res.normMismatch[k] = rho / lam2 - 1.;
      res.pairOverlap[k] = (nwScaled - ne) / (nwScaled + ne);
      d = y - lam2 * e;
      res.resStored[k] = std::sqrt(norm2(d) / ne) / lam2;
      d = z + lam2 * e;
      res.resFull[k] = std::sqrt(norm2(d) / (ne + nwScaled)) / lam;
      if (mu != 0.) {
        res.evalShiftOverMass2[k] = (lam2 - rho) / (mu * mu);
      }
    }
    d = y - rho * e;
    res.resRayleigh[k] = std::sqrt(norm2(d) / ne) / rho;
    d = z + y;
    res.antiHermiticity[k] = std::sqrt(norm2(d) / norm2(y));

    // orthogonality to the preceding gramBand modes, on each checkerboard
    if (band > 0) {
      double gE = 0., gW = 0.;
      const int nPrev = std::min<int>(band, k);

      for (int b = 0; b < nPrev; ++b) {
        const unsigned int s = (k - 1 - b) % band;

        gE = std::max(gE, cosine(ringE[s], e, ringNe[s], ne));
        gW = std::max(gW, cosine(ringW[s], w, ringNw[s], nw));
      }
      if (k > 0) {
        res.gramStored[k] = gE;
        res.gramRebuilt[k] = gW;
      }

      const unsigned int s = k % band;
      ringE[s] = e;
      ringW[s] = w;
      ringNe[s] = ne;
      ringNw[s] = nw;
    }
  }

  // StagLMAFull normalises every mode with mode 0's norm
  std::vector<double> normRelMode0(nEigs);
  for (int k = 0; k < nEigs; ++k) {
    normRelMode0[k] = res.normStored[k] / res.normStored[0] - 1.;
  }

  LOG(Message) << "Diagonal 2m = " << res.massDiag << ", lambda from "
               << res.lambdaStored.front() << " to " << res.lambdaStored.back()
               << std::endl;
  eigenPackCheckSummary("normMismatch", res.normMismatch, res.lambdaStored);
  eigenPackCheckSummary("pairOverlap", res.pairOverlap, res.lambdaStored);
  eigenPackCheckSummary("resStored", res.resStored, res.lambdaStored);
  eigenPackCheckSummary("resRayleigh", res.resRayleigh, res.lambdaStored);
  eigenPackCheckSummary("resFull", res.resFull, res.lambdaStored);
  eigenPackCheckSummary("evalShiftOverMass2", res.evalShiftOverMass2,
                        res.lambdaStored);
  eigenPackCheckSummary("antiHermiticity", res.antiHermiticity,
                        res.lambdaStored);
  eigenPackCheckSummary("normRelMode0", normRelMode0, res.lambdaStored);
  if (band > 0) {
    eigenPackCheckSummary("gramStored", res.gramStored, res.lambdaStored);
    eigenPackCheckSummary("gramRebuilt", res.gramRebuilt, res.lambdaStored);
  }

  saveResult(par().output, "eigenPackCheck", res);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MUtilities_EigenPackCheck_hpp_
