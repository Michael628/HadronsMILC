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
#include <GridMilc/GridMilc.h>

#include <map>
#include <string>
#include <vector>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     TGammaMap: label-keyed per-gamma maps                   *
 ******************************************************************************/
// Cross-module gamma association: the key is each gamma's RAW label
// (StagGamma::getLabelName -- naming follows the label, physics follows
// the object), the value the per-gamma product (propagator fields,
// per-timeslice maps, ...). Entries are constructed from the gammas
// module's vector INSIDE the envCreate window: the environment's memory
// profiler sizes objects by the allocation delta around construction
// (Hadrons/Environment.hpp createObject), so post-hoc fills would be
// invisible to the scheduler's peak-memory objective. Duplicate labels
// dedupe by first occurrence (map::emplace no-op on existing keys),
// matching the first-wins semantics of the former std::map::insert keys.
template <typename T>
class TGammaMap : public std::map<std::string, T> {
public:
  TGammaMap(void) = default;

  // One entry per gamma label; every value copy-constructed from the
  // same ctor arguments (a GridBase* for lattices, a size and a
  // GridBase* for per-source vectors).
  template <typename... Args>
  TGammaMap(const std::vector<StagGamma> &gammas, Args &&...args)
      : std::map<std::string, T>() {
    for (auto &g : gammas) {
      this->emplace(g.getLabelName(), args...);
    }
  }
};

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
  // The ONLY place spin-taste operators are built: MakeSpinTasteOps parses,
  // constructs, binds the gauge and applies the eps fold per object, in
  // string order (sequential strToVec -- directional lists keep index
  // semantics). applyG5=true yields immutable eps-folded objects: consumers
  // must bind them through const& (setSpin/setTaste/setSpinTaste re-derive
  // _negated from the stored P pair and silently destroy the fold).
  LatticeGaugeField *U = nullptr;
  if (!par().gauge.empty()) {
    U = &envGet(LatticeGaugeField, par().gauge);
  }
  auto ops = StagGamma::MakeSpinTasteOps(par().gammas, par().applyG5, U);

  if (ops.empty()) {
    LOG(Warning) << "SpinTaste module '" << getName()
                 << "': empty gamma list; publishing an empty vector"
                 << std::endl;
  } else {
    LOG(Message) << "Publishing " << ops.size()
                 << " spin-taste operator(s) (applyG5 "
                 << (par().applyG5 ? "true" : "false") << ", gauge '"
                 << par().gauge << "')" << std::endl;
  }

  envCreate(std::vector<StagGamma>, getName(), 1, ops);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl> void TSpinTasteMILC<FImpl>::execute(void) {}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_SpinTaste_hpp_
