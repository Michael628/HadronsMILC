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
// Cross-module gamma association: the key is each gamma's EFFECTIVE label
// (a SpinTaste module's custom `labels` override when given, else
// StagGamma::getLabelName() -- naming follows the label, physics follows
// the object), the value the per-gamma product (propagator fields,
// per-timeslice maps, ...). Entries are constructed from the gammas
// module's vector or label map INSIDE the envCreate window: the
// environment's memory profiler sizes objects by the allocation delta
// around construction (Hadrons/Environment.hpp createObject), so post-hoc
// fills would be invisible to the scheduler's peak-memory objective.
// Duplicate labels dedupe by first occurrence (map::emplace no-op on
// existing keys), matching the first-wins semantics of the former
// std::map::insert keys.
template <typename T>
class TGammaMap : public std::map<std::string, T> {
public:
  TGammaMap(void) = default;

  // One entry per gamma label; every value copy-constructed from the
  // same ctor arguments (a GridBase* for lattices, a size and a
  // GridBase* for per-source vectors). Keys by StagGamma::getLabelName()
  // -- used only where no label map is available (should not occur for
  // producer modules downstream of a SpinTaste module's `_map` output;
  // kept for direct StagGamma::MakeSpinTasteOps() callers with no
  // custom-label concept).
  template <typename... Args>
  TGammaMap(const std::vector<StagGamma> &gammas, Args &&...args)
      : std::map<std::string, T>() {
    for (auto &g : gammas) {
      this->emplace(g.getLabelName(), args...);
    }
  }

  // One entry per (label, StagGamma) pair already resolved by a
  // producer -- reuses the label map's OWN keys instead of re-deriving
  // via getLabelName(), so a SpinTaste module's custom `labels` override
  // flows through unchanged. This is the constructor GaugeProp and
  // LMAMesonFieldProp use against a SpinTaste `_map` object.
  template <typename... Args>
  TGammaMap(const TGammaMap<StagGamma> &labelMap, Args &&...args)
      : std::map<std::string, T>() {
    for (auto &p : labelMap) {
      this->emplace(p.first, args...);
    }
  }
};

/******************************************************************************
 *                                 SpinTasteMILC *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

// Hadrons-local Par: nests Grid's SpinTasteParams (gammas/gauge/applyG5,
// StagGamma.h:23-27) unchanged, plus a sibling `labels` override -- keeps
// the custom-label mechanism entirely on the Hadrons side, no StagGamma/
// Grid changes. Modeled on ImplicitlyRestartedLanczos.hpp:18-27's
// LanczosParams-nested-in-ImplicitlyRestartedLanczosMILCPar precedent.
//
// labels: optional whitespace-separated list, positionally parallel to
//         the parsed spinTaste.gammas list; "" = every gamma defaults to
//         StagGamma::getLabelName(). When non-empty, size must match the
//         parsed gamma count (fatal otherwise). The EFFECTIVE label
//         (custom or default) is what keys this module's `_map` output
//         and feeds the duplicate-label check -- naming follows the
//         effective label, physics still follows the StagGamma object.
class SpinTasteMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(SpinTasteMILCPar, SpinTasteParams,
                                  spinTaste, std::string, labels);
};

template <typename FImpl>
class TSpinTasteMILC : public Module<SpinTasteMILCPar> {
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
  // the effective label for ops[i]: par().labels[i] when the (validated,
  // count-matched) override list is non-empty, else ops[i].getLabelName().
  // Single derivation point for setup()'s opsMap keying AND the
  // duplicate-label check, so the two can never disagree.
  std::vector<std::string> effectiveLabels(
      const std::vector<StagGamma> &ops) const;
};

MODULE_REGISTER_TMP(SpinTaste, TSpinTasteMILC<STAGIMPL>, MFermion);

/******************************************************************************
 *                          TSpinTasteMILC implementation *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSpinTasteMILC<FImpl>::TSpinTasteMILC(const std::string name)
    : Module<SpinTasteMILCPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSpinTasteMILC<FImpl>::getInput(void) {
  std::vector<std::string> in;

  if (!par().spinTaste.gauge.empty()) {
    in.push_back(par().spinTaste.gauge);
  }

  return in;
}

template <typename FImpl>
std::vector<std::string> TSpinTasteMILC<FImpl>::getOutput(void) {
  std::vector<std::string> out = {getName(), getName() + "_map"};

  return out;
}

template <typename FImpl>
DependencyMap TSpinTasteMILC<FImpl>::getObjectDependencies(void) {
  DependencyMap dep;

  if (!par().spinTaste.gauge.empty()) {
    dep.insert({par().spinTaste.gauge, getName()});
  }

  return dep;
}

// effective-label derivation //////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSpinTasteMILC<FImpl>::effectiveLabels(
    const std::vector<StagGamma> &ops) const {
  auto overrides = strToVec<std::string>(par().labels);
  std::vector<std::string> labels;
  labels.reserve(ops.size());

  if (overrides.empty()) {
    for (auto &g : ops) {
      labels.push_back(g.getLabelName());
    }
    return labels;
  }

  if (overrides.size() != ops.size()) {
    HADRONS_ERROR(Argument,
                  "SpinTaste module '" + getName() + "': 'labels' has " +
                      std::to_string(overrides.size()) +
                      " entries but 'spinTaste.gammas' parsed " +
                      std::to_string(ops.size()) +
                      " gamma(s) -- labels must be a positionally parallel "
                      "list, one entry per gamma");
  }

  return overrides;
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
  if (!par().spinTaste.gauge.empty()) {
    U = &envGet(LatticeGaugeField, par().spinTaste.gauge);
  }
  auto ops = StagGamma::MakeSpinTasteOps(par().spinTaste.gammas,
                                        par().spinTaste.applyG5, U);

  if (ops.empty()) {
    LOG(Warning) << "SpinTaste module '" << getName()
                 << "': empty gamma list; publishing an empty vector"
                 << std::endl;
  } else {
    LOG(Message) << "Publishing " << ops.size()
                 << " spin-taste operator(s) (applyG5 "
                 << (par().spinTaste.applyG5 ? "true" : "false") << ", gauge '"
                 << par().spinTaste.gauge << "')" << std::endl;
  }

  // Single derivation point: the effective label (custom override or
  // default getLabelName()) feeds BOTH the duplicate check below and the
  // opsMap keying -- they can never disagree.
  auto labels = effectiveLabels(ops);

  // Centralized duplicate-label validation on the EFFECTIVE label (not
  // getLabelName() directly): a duplicate label is caught here before any
  // downstream module can act on it. Intentional label sharing ACROSS
  // separate SpinTaste module instances is unaffected -- only within-list
  // duplicates are fatal here.
  for (unsigned int i = 0; i < labels.size(); ++i) {
    for (unsigned int j = i + 1; j < labels.size(); ++j) {
      if (labels[i] == labels[j]) {
        HADRONS_ERROR(Argument,
                      "duplicate gamma label '" + labels[i] +
                          "' in SpinTaste module '" + getName() +
                          "' (gamma-map entries would collide)");
      }
    }
  }

  envCreate(std::vector<StagGamma>, getName(), 1, ops);

  // Companion output: one TGammaMap<StagGamma> entry per gamma, keyed by
  // its EFFECTIVE label (custom override or default getLabelName()).
  TGammaMap<StagGamma> opsMap;
  for (unsigned int i = 0; i < ops.size(); ++i) {
    opsMap.emplace(labels[i], ops[i]);
  }
  envCreate(TGammaMap<StagGamma>, getName() + "_map", 1, opsMap);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl> void TSpinTasteMILC<FImpl>::execute(void) {}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_SpinTaste_hpp_
