/*
 * ModifyEigenPackMILC.hpp, part of Hadrons
 * (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef HadronsMILC_MUtilities_ModifyEigenPack_hpp_
#define HadronsMILC_MUtilities_ModifyEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Load eigen vectors/values package                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class ModifyEigenPackMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ModifyEigenPackMILCPar, std::string,
                                  eigenPack, bool, evenEigen, double, mass);
};

template <typename FImpl, typename Pack>
class TModifyEigenPackMILC : public Module<ModifyEigenPackMILCPar> {
public:
  typedef typename Pack::Field Field;
  typedef BaseEigenPack<Field> BasePack;

  FERM_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TModifyEigenPackMILC(const std::string name);
  // destructor
  virtual ~TModifyEigenPackMILC(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  virtual DependencyMap getObjectDependencies(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(
    ModifyEigenPackMILC,
    ARG(TModifyEigenPackMILC<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),
    MUtilities);

/******************************************************************************
 *                    TModifyEigenPackMILC implementation *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TModifyEigenPackMILC<FImpl, Pack>::TModifyEigenPackMILC(const std::string name)
    : Module<ModifyEigenPackMILCPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TModifyEigenPackMILC<FImpl, Pack>::getInput(void) {
  std::vector<std::string> in = {par().eigenPack};

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TModifyEigenPackMILC<FImpl, Pack>::getOutput(void) {
  std::vector<std::string> out = {getName()};

  return out;
}

template <typename FImpl, typename Pack>
DependencyMap TModifyEigenPackMILC<FImpl, Pack>::getObjectDependencies(void) {
  DependencyMap dep;

  dep.insert({par().eigenPack, getName()});

  return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TModifyEigenPackMILC<FImpl, Pack>::setup(void) {
  auto Ls = env().getObjectLs(par().eigenPack);

  auto &epack = envGet(BasePack, par().eigenPack);

  envCreate(MassShiftEigenPack<FImpl>, getName(), Ls, epack.evec, epack.eval,
            par().mass);
  // envTmp(FermionField, "tempRb", 1, envGetRbGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TModifyEigenPackMILC<FImpl, Pack>::execute(void) {
  int Ls = env().getObjectLs(par().eigenPack);
  auto &epack = envGet(BasePack, par().eigenPack);

  bool evenEigen = par().evenEigen;
  double mass = par().mass;

  LOG(Message) << "Creating modified eigenpack with eigenvectors from "
               << par().eigenPack
               << " and Dirac matrix eigenvalues mass + i*lambda_D"
               << std::endl;

  int cb = (evenEigen ? Even : Odd);
  /*    if (!par().checkerSwapAction.empty()) {
          LOG(Message) << "Swapping checkerboard from " << (evenEigen?"Odd to
     ":"Even to ") << (evenEigen?"Even":"Odd")<< std::endl;

          evec.resize(epack.evec.size(),envGetRbGrid(Field, Ls));

          auto &action = envGet(FMat, par().checkerSwapAction);
          int cbNeg = (!evenEigen ? Even : Odd);

          envGetTmp(FermionField,tempRb);
          tempRb = Zero();
          tempRb.Checkerboard() = cb;

          for (int i = 0; i < epack.evec.size();i++) {

              epack.evec[i].Checkerboard() = cbNeg;
              action.Meooe(epack.evec[i],tempRb);

              if (mass == 0.0) {
                  evec[i] = (1.0/evalM[i])*tempRb;
              } else {
                  evec[i] = (1.0/(evalM[i]-evalM[i].real()))*tempRb;
              }
          }
    */

  // Bad design. Possibly changing checkerboard of another module's data.
  // Should probably add this functionality to LoadEigenPack module.
  for (auto &e : epack.evec) {
    e.Checkerboard() = cb;
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_ModifyEigenPack_hpp_
