/*
 * SpinTasteMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2021
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
 * Author: Michael Lynch <michaellynch628@gmail.com>
 * Author: Carleton DeTar <detar@physics.utah.edu>
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

#ifndef HadronsMILC_MFermion_SpinTaste_hpp_
#define HadronsMILC_MFermion_SpinTaste_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <spin/StagGamma.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 SpinTasteMILC *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

template <typename FImpl>
class TSpinTasteMILC : public Module<SpinTasteParams> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TSpinTasteMILC(const std::string name);
  // destructor
  virtual ~TSpinTasteMILC(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  virtual DependencyMap getObjectDependencies(void);

protected:
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);

private:
  bool hasPhase_{false};
  std::string phName_;
};

MODULE_REGISTER_TMP(SpinTaste, TSpinTasteMILC<STAGIMPL>, MFermion);

/******************************************************************************
 *                          TSpinTasteMILC implementation *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSpinTasteMILC<FImpl>::TSpinTasteMILC(const std::string name)
    : Module<SpinTasteParams>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSpinTasteMILC<FImpl>::getInput(void) {
  std::vector<std::string> in;

  if (!par().gauge.empty()) {
    in.push_back(par().gauge);
  }

  return in;
}

template <typename FImpl>
std::vector<std::string> TSpinTasteMILC<FImpl>::getOutput(void) {
  std::vector<std::string> out = {getName()};

  return out;
}

template <typename FImpl>
DependencyMap TSpinTasteMILC<FImpl>::getObjectDependencies(void) {
  DependencyMap dep;

  if (!par().gauge.empty()) {
    dep.insert({par().gauge, getName()});
  }

  return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl> void TSpinTasteMILC<FImpl>::setup(void) {
  auto gammaList = strToVec<StagGamma::SpinTastePair>(par().gammas);
  envCreate(std::vector<StagGamma>, getName(), 1, gammaList.size(),
            StagGamma());

  auto &spinTaste = envGet(std::vector<StagGamma>, getName());

  for (int i = 0; i < gammaList.size(); i++) {
    spinTaste[i].setSpinTaste(gammaList[i]);
    if (!par().gauge.empty()) {
      auto &gauge = envGet(LatticeGaugeField, par().gauge);
      spinTaste[i].setGaugeField(gauge);
    }
  }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl> void TSpinTasteMILC<FImpl>::execute(void) {}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_SpinTaste_hpp_
