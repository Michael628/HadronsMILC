/*
 * EigenPackFullPairs.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MUtilities_EigenPackFullPairs_hpp_
#define HadronsMILC_MUtilities_EigenPackFullPairs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *        Full-volume |e+o>/|e-o> pair pack from a checkerboarded pack       *
 ******************************************************************************/
/*  Builds the production-layout eigenvector representation used by
    MIO::LoadMesonField / MFermion::StagLMAMesonFieldProp files: for each
    checkerboarded eigenvector e_E^k (parity cb, eigenvalue 2m + i*lam_k),
    construct the full-volume pair

        |e+o>_k = (e_E^k +  Meooe(e_E^k)/(i*lam_k)) / sqrt(2)
        |e-o>_k = (e_E^k -  Meooe(e_E^k)/(i*lam_k)) / sqrt(2)

    (the same e_O relation as MesonField's swapEvecCheckerFn). Outputs:
      <name>_vec  std::vector<FermionField>, 2*nEvec full-grid fields
                  (rows 2k/2k+1 of eigenpair k) -- usable directly as a
                  MesonField 'left' probe set
      <name>      MassShiftEigenPack referencing <name>_vec -- pass as
                  MesonField 'lowModes' with an EMPTY action to write
                  full-volume |e+o>/|e-o> pair files (norm2_full = 1 per
                  pair for unit-checkerboard evecs, so the file's baked norm
                  is 1 and StagLMAMesonFieldProp's default pairScale sqrt(2)
                  applies)
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EigenPackFullPairsPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(EigenPackFullPairsPar,
                                  std::string, eigenPack,
                                  std::string, action);
};

template <typename FImpl, typename Pack>
class TEigenPackFullPairs : public Module<EigenPackFullPairsPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TEigenPackFullPairs(const std::string name);
  // destructor
  virtual ~TEigenPackFullPairs(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  virtual DependencyMap getObjectDependencies(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(EigenPackFullPairs,
                    ARG(TEigenPackFullPairs<STAGIMPL,
                                            MassShiftEigenPack<STAGIMPL>>),
                    MUtilities);

/******************************************************************************
 *                  TEigenPackFullPairs implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TEigenPackFullPairs<FImpl, Pack>::TEigenPackFullPairs(const std::string name)
    : Module<EigenPackFullPairsPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TEigenPackFullPairs<FImpl, Pack>::getInput(void) {
  std::vector<std::string> in{par().eigenPack, par().action};

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TEigenPackFullPairs<FImpl, Pack>::getOutput(void) {
  std::vector<std::string> out{getName(), getName() + "_vec"};

  return out;
}

template <typename FImpl, typename Pack>
DependencyMap TEigenPackFullPairs<FImpl, Pack>::getObjectDependencies(void) {
  DependencyMap dep;

  dep.insert({par().eigenPack, getName()});
  dep.insert({par().action, getName()});

  return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TEigenPackFullPairs<FImpl, Pack>::setup(void) {
  auto &epack = envGet(Pack, par().eigenPack);
  int Ls = env().getObjectLs(par().eigenPack);
  unsigned int nEvec = epack.evec.size();

  envCreate(std::vector<FermionField>, getName() + "_vec", Ls, 2 * nEvec,
            envGetGrid(FermionField));

  // rebuild the raw M†M eigenvalues (lambda^2) expected by the adaptor
  // constructor, one per row (each pair duplicates its eigenpair's value)
  std::vector<RealD> evalMm(2 * nEvec, 0.);
  for (unsigned int k = 0; k < nEvec; ++k) {
    RealD lam2 = epack.eval[k].imag() * epack.eval[k].imag();

    evalMm[2 * k] = lam2;
    evalMm[2 * k + 1] = lam2;
  }
  envCreate(Pack, getName(), Ls,
            envGet(std::vector<FermionField>, getName() + "_vec"), evalMm,
            epack.mass);

  envTmp(FermionField, "tempOdd", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "tempEven", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "tempScaled", 1, envGetRbGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TEigenPackFullPairs<FImpl, Pack>::execute(void) {
  auto &epack = envGet(Pack, par().eigenPack);
  auto &pairs = envGet(std::vector<FermionField>, getName() + "_vec");
  auto &mat = envGet(FMat, par().action);
  RealD invSqrt2 = 1. / std::sqrt(2.);

  envGetTmp(FermionField, tempOdd);
  envGetTmp(FermionField, tempEven);
  envGetTmp(FermionField, tempScaled);

  LOG(Message) << "Building " << 2 * epack.evec.size()
               << " full-volume |e+o>/|e-o> pairs from '"
               << par().eigenPack << "'" << std::endl;

  for (unsigned int k = 0; k < epack.evec.size(); ++k) {
    const FermionField &e = epack.evec[k];
    int cb = e.Checkerboard();
    int cbNeg = (cb == Even) ? Odd : Even;
    ComplexD eval_D = ComplexD(0., epack.eval[k].imag());

    // e_O = Meooe(e_E) / (i*lam) -- the swapEvecCheckerFn relation
    tempOdd = Zero();
    tempOdd.Checkerboard() = cbNeg;
    mat.Meooe(e, tempOdd);
    tempOdd = (1.0 / eval_D) * tempOdd;

    tempEven = invSqrt2 * e;
    tempScaled = invSqrt2 * tempOdd;
    // expression-template assignments do NOT propagate the checkerboard
    // attribute (Grid Lattice_arith/ET), and setCheckerboard(full, half)
    // dispatches on half's attribute -- label the halves explicitly or the
    // |e-o> content lands on the wrong full-grid parity
    tempEven.Checkerboard() = cb;
    tempScaled.Checkerboard() = cbNeg;

    // |e+o> = (e_E + e_O)/sqrt(2)
    pairs[2 * k] = Zero();
    setCheckerboard(pairs[2 * k], tempEven);
    setCheckerboard(pairs[2 * k], tempScaled);

    // |e-o> = (e_E - e_O)/sqrt(2)
    tempScaled = -tempScaled;
    pairs[2 * k + 1] = Zero();
    setCheckerboard(pairs[2 * k + 1], tempEven);
    setCheckerboard(pairs[2 * k + 1], tempScaled);
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MUtilities_EigenPackFullPairs_hpp_
