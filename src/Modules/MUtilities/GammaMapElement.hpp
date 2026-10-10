/*
 * GammaMapElement.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MUtilities_GammaMapElement_hpp_
#define HadronsMILC_MUtilities_GammaMapElement_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <GridMilc/GridMilc.h>
#include <Modules/MFermion/SpinTaste.hpp> // TGammaMap + gammas module

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             TGammaMap element -> bare environment object bridge            *
 ******************************************************************************/
/*  Bridges one label-keyed TGammaMap element onto a bare environment
    object under this module's own name, for consumers that take plain
    fields rather than maps (the MUtilities::PropToFermions 'prop'
    parameter, StagA2AMesonField-style right-hand sides, ...). Copies
    ONE element per module instance -- probe/test scale by design;
    map-native consumers should read the map directly.

    map         name of the TGammaMap<T> environment object (e.g. a
                StagGaugeProp output, or a StagLMAMesonFieldProp
                per-timeslice map "<name>_t<t>")
    label       the element key (a gamma's RAW label,
                StagGamma::getLabelName)

    Module lifecycle: the input edge on 'map' orders this module after
    the map's producer (Hadrons runs each module's setup()+execute()
    back-to-back in program order, Module.hpp operator()), so by this
    module's setup the entry already holds the producer's values in a
    real run; setup allocates the bare output by copy-constructing the
    entry inside the envCreate window (the copy is zeroed only under
    the scheduler's dry-run memory profile, which recursively profiles
    producers first -- the allocation stays profiler-visible either
    way), and execute re-assigns the element as belt-and-braces against
    scheduling variation. Consumers of this module name it in their own
    getInput, ordering them after the copy.
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class GammaMapElementPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(GammaMapElementPar, std::string, map,
                                  std::string, label);
};

template <typename FImpl>
class TGammaMapElement : public Module<GammaMapElementPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TGammaMapElement(const std::string name);
  // destructor
  virtual ~TGammaMapElement(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  template <typename T> void setupHelper(void);
  virtual void setup(void);
  // execution
  template <typename T> void executeHelper(void);
  virtual void execute(void);

private:
  // find an element with a loud error listing the available keys
  // (fails at setup, before any consumer executes)
  template <typename T>
  typename TGammaMap<T>::iterator findElement(TGammaMap<T> &map) const;
};

MODULE_REGISTER_TMP(GammaMapElement, TGammaMapElement<STAGIMPL>, MUtilities);

/******************************************************************************
 *                      TGammaMapElement implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGammaMapElement<FImpl>::TGammaMapElement(const std::string name)
    : Module<GammaMapElementPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGammaMapElement<FImpl>::getInput(void) {
  std::vector<std::string> in = {par().map};

  return in;
}

template <typename FImpl>
std::vector<std::string> TGammaMapElement<FImpl>::getOutput(void) {
  std::vector<std::string> out = {getName()};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename T>
typename TGammaMap<T>::iterator
TGammaMapElement<FImpl>::findElement(TGammaMap<T> &map) const {
  auto it = map.find(par().label);
  if (it == map.end()) {
    std::string available;
    for (auto &p : map) {
      available += (available.empty() ? "" : ", ") + p.first;
    }
    HADRONS_ERROR(Argument, "gamma map '" + par().map + "' has no entry '" +
                                par().label + "' (available: " + available +
                                ")");
  }

  return it;
}

template <typename FImpl>
template <typename T> void TGammaMapElement<FImpl>::setupHelper(void) {
  // validate the key at setup, before any consumer executes (setup
  // runs after the producer's via the input edge on 'map')
  auto &map = envGet(TGammaMap<T>, par().map);
  auto it = findElement(map);

  // allocation-only output: copy-construct from the entry inside the
  // envCreate window -- the memory profiler sizes objects by the
  // allocation delta around construction; the VALUE assignment happens
  // at execute as belt-and-braces
  envCreate(T, getName(), 1, it->second);
}

template <typename FImpl> void TGammaMapElement<FImpl>::setup(void) {
  if (par().map.empty() || par().label.empty()) {
    HADRONS_ERROR(Argument,
                  "GammaMapElement requires both 'map' and 'label'");
  }
  if (envHasType(TGammaMap<PropagatorField>, par().map)) {
    setupHelper<PropagatorField>();
  } else if (envHasType(TGammaMap<std::vector<PropagatorField>>, par().map)) {
    setupHelper<std::vector<PropagatorField>>();
  } else if (envHasType(TGammaMap<FermionField>, par().map)) {
    setupHelper<FermionField>();
  } else if (envHasType(TGammaMap<std::vector<FermionField>>, par().map)) {
    setupHelper<std::vector<FermionField>>();
  } else {
    HADRONS_ERROR(Argument, "'map' parameter '" + par().map +
                                "' is not a gamma map (expected "
                                "TGammaMap<PropagatorField>, "
                                "TGammaMap<std::vector<PropagatorField>>, "
                                "TGammaMap<FermionField> or "
                                "TGammaMap<std::vector<FermionField>>)");
  }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename T> void TGammaMapElement<FImpl>::executeHelper(void) {
  auto &map = envGet(TGammaMap<T>, par().map);
  auto it = findElement(map);

  envGet(T, getName()) = it->second;
  LOG(Message) << "Copied element '" << par().label << "' of gamma map '"
               << par().map << "'" << std::endl;
}

template <typename FImpl> void TGammaMapElement<FImpl>::execute(void) {
  if (envHasType(TGammaMap<PropagatorField>, par().map)) {
    executeHelper<PropagatorField>();
  } else if (envHasType(TGammaMap<std::vector<PropagatorField>>, par().map)) {
    executeHelper<std::vector<PropagatorField>>();
  } else if (envHasType(TGammaMap<FermionField>, par().map)) {
    executeHelper<FermionField>();
  } else if (envHasType(TGammaMap<std::vector<FermionField>>, par().map)) {
    executeHelper<std::vector<FermionField>>();
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MUtilities_GammaMapElement_hpp_
