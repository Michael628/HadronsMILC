/*
 * GaugePropMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Guido Cossu <guido.cossu@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Nils Asmussen <n.asmussen@soton.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: pretidav <david.preti@csic.es>
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

#ifndef HadronsMILC_MFermion_GaugeProp_hpp_
#define HadronsMILC_MFermion_GaugeProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <GridMilc/GridMilc.h>

#include <Modules/MFermion/SpinTaste.hpp> // TGammaMap + the gammas module

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                GaugePropMILC   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class GaugePropMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(GaugePropMILCPar, std::string, source,
                                  std::string, gammas, std::string, solver,
                                  std::string, guess, std::string, sourceLabel);
};
// gammas:   REQUIRED name of an MFermion::SpinTaste module; this module
//           consumes the SpinTaste module's `_map` companion output
//           (par().gammas + "_map", a TGammaMap<StagGamma>) -- never the
//           bare vector. Empty is a setup-time error (StagGaugePropLegacy
//           preserves the former empty-gammas identity behavior).
// guess:    name of a TGammaMap object (one entry per gamma label, e.g.
//           an LMAMesonFieldProp per-timeslice map); per-gamma guess is
//           guessMap.at(gammaLabel)
// sourceLabel: when source is a TGammaMap, the element key used for EVERY
//           gamma; "" = each gamma reads its own label's element

template <typename FImpl>
class TGaugePropMILC : public Module<GaugePropMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );
  SOLVER_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TGaugePropMILC(const std::string name);
  // destructor
  virtual ~TGaugePropMILC(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);

protected:
  // setup
  template <typename TField>
  void setupHelper(const TGammaMap<StagGamma> &gammas);
  virtual void setup(void);
  // execution
  template <typename TField>
  EnableIf<is_lattice<TField>, void>
  executeHelper(TField &sol, const TField &src, const StagGamma &gamma,
                const TField *guess = nullptr);
  template <typename TField>
  EnableIf<is_lattice<TField>, void>
  executeHelper(std::vector<TField> &sol, const std::vector<TField> &src,
                const StagGamma &gamma, const std::vector<TField> *guess = nullptr);
  // one overload per source container shape; each feeds gammaLoop with a
  // per-gamma source fetcher
  template <typename TField> void executeHelper(const TField &src);
  template <typename TField> void executeHelper(const std::vector<TField> &src);
  template <typename TField> void executeHelper(const TGammaMap<TField> &srcMap);
  template <typename TField>
  void executeHelper(const TGammaMap<std::vector<TField>> &srcMap);
  // the per-gamma loop: V is the per-gamma value type (TField, or
  // std::vector<TField> for per-source-element solves)
  template <typename V, typename SrcFetch> void gammaLoop(SrcFetch &&srcFor);
  virtual void execute(void);

private:
  void solveField(FermionField &prop, const FermionField &src,
                  const FermionField *guess = nullptr);
  void solveField(PropagatorField &prop, const PropagatorField &src,
                  const PropagatorField *guess = nullptr);

private:
  bool _hasGuess;
};

MODULE_REGISTER_TMP(StagGaugeProp, TGaugePropMILC<STAGIMPL>, MFermion);

/******************************************************************************
 *                      TGaugePropMILC implementation   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaugePropMILC<FImpl>::TGaugePropMILC(const std::string name)
    : Module<GaugePropMILCPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGaugePropMILC<FImpl>::getInput(void) {
  std::vector<std::string> in = {par().source, par().solver,
                                 par().gammas + "_map"};

  if (!par().guess.empty()) {
    in.push_back(par().guess);
  }

  return in;
}

template <typename FImpl>
std::vector<std::string> TGaugePropMILC<FImpl>::getOutput(void) {
  std::vector<std::string> out = {getName()};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename TField>
void TGaugePropMILC<FImpl>::setupHelper(const TGammaMap<StagGamma> &gammas) {
  envTmpLat(TField, "field");

  if (envHasType(TField, par().source) ||
      envHasType(TGammaMap<TField>, par().source)) {
    envCreate(TGammaMap<TField>, getName(), 1, gammas, envGetGrid(TField));
    for (auto &p : envGet(TGammaMap<TField>, getName())) {
      p.second = Zero();
    }
  } else if (envHasType(TGammaMap<std::vector<TField>>, par().source)) {
    auto &src = envGet(TGammaMap<std::vector<TField>>, par().source);
    envCreate(TGammaMap<std::vector<TField>>, getName(), 1, gammas,
              std::vector<TField>(src.begin()->second.size(),
                                  envGetGrid(TField)));
    for (auto &p : envGet(TGammaMap<std::vector<TField>>, getName())) {
      for (auto &s : p.second) {
        s = Zero();
      }
    }
  } else {
    auto &src = envGet(std::vector<TField>, par().source);
    envCreate(TGammaMap<std::vector<TField>>, getName(), 1, gammas,
              std::vector<TField>(src.size(), envGetGrid(TField)));
    for (auto &p : envGet(TGammaMap<std::vector<TField>>, getName())) {
      for (auto &s : p.second) {
        s = Zero();
      }
    }
  }
}

template <typename FImpl> void TGaugePropMILC<FImpl>::setup(void) {
  _hasGuess = !par().guess.empty();

  if (par().gammas.empty()) {
    HADRONS_ERROR(Argument,
                  "GaugeProp requires the 'gammas' SpinTaste module name "
                  "(use StagGaugePropLegacy for the former empty-gammas "
                  "identity behavior)");
  }
  const auto &gammas = envGet(TGammaMap<StagGamma>, par().gammas + "_map");
  if (gammas.empty()) {
    HADRONS_ERROR(Argument,
                  "SpinTaste module '" + par().gammas + "' published an "
                  "empty gamma map -- GaugeProp requires at least one "
                  "gamma");
  }

  if (envHasType(PropagatorField, par().source) ||
      envHasType(std::vector<PropagatorField>, par().source) ||
      envHasType(TGammaMap<PropagatorField>, par().source) ||
      envHasType(TGammaMap<std::vector<PropagatorField>>, par().source)) {

    // Additional temp storage for propagator field calculations
    envTmpLat(FermionField, "fermIn");
    envTmpLat(FermionField, "fermOut");
    envTmpLat(FermionField, "fermGuess");
    envGetTmp(FermionField, fermIn);
    envGetTmp(FermionField, fermOut);
    envGetTmp(FermionField, fermGuess);
    fermIn = Zero();
    fermOut = Zero();
    fermGuess = Zero();

    setupHelper<PropagatorField>(gammas);

  } else if (envHasType(FermionField, par().source) ||
             envHasType(std::vector<FermionField>, par().source) ||
             envHasType(TGammaMap<FermionField>, par().source) ||
             envHasType(TGammaMap<std::vector<FermionField>>, par().source)) {
    setupHelper<FermionField>(gammas);
  } else {
    HADRONS_ERROR(Logic,
                  "Type of source '" + par().source + "' not recognized.");
  }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugePropMILC<FImpl>::solveField(FermionField &sol,
                                       const FermionField &src,
                                       const FermionField *guess) {
  auto &solver = envGet(Solver, par().solver);

  if (guess != nullptr) {
    solver(sol, src, *guess);
  } else {
    solver(sol, src);
  }
}

template <typename FImpl>
void TGaugePropMILC<FImpl>::solveField(PropagatorField &sol,
                                       const PropagatorField &src,
                                       const PropagatorField *guess) {
  auto &solver = envGet(Solver, par().solver);

  envGetTmp(FermionField, fermIn);
  envGetTmp(FermionField, fermOut);
  envGetTmp(FermionField, fermGuess);

  for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
    PropToFerm<FImpl>(fermIn, src, c);
    if (guess != nullptr) {
      PropToFerm<FImpl>(fermGuess, *guess, c);
      solver(fermOut, fermIn, fermGuess);
    } else {
      solver(fermOut, fermIn);
    }
    FermToProp<FImpl>(sol, fermOut, c);
  }
}

template <typename FImpl>
template <typename TField>
EnableIf<is_lattice<TField>, void>
TGaugePropMILC<FImpl>::executeHelper(TField &sol, const TField &src,
                                     const StagGamma &gamma,
                                     const TField *guess) {
  envGetTmp(TField, field);

  gamma(field, src);

  if (_hasGuess) {
    solveField(sol, field, guess);
  } else {
    solveField(sol, field);
  }
}

template <typename FImpl>
template <typename TField>
EnableIf<is_lattice<TField>, void> TGaugePropMILC<FImpl>::executeHelper(
    std::vector<TField> &sol, const std::vector<TField> &src,
    const StagGamma &gamma, const std::vector<TField> *guess) {
  envGetTmp(TField, field);

  for (unsigned int i = 0; i < src.size(); i++) {

    gamma(field, src[i]);

    if (_hasGuess) {
      const TField *guessTemp = &(guess->at(i));
      solveField(sol[i], field, guessTemp);
    } else {
      solveField(sol[i], field);
    }
  }
}

template <typename FImpl>
template <typename V, typename SrcFetch>
void TGaugePropMILC<FImpl>::gammaLoop(SrcFetch &&srcFor) {
  const auto &gammas = envGet(TGammaMap<StagGamma>, par().gammas + "_map");
  auto &solMap = envGet(TGammaMap<V>, getName());
  const TGammaMap<V> *guessMap = nullptr;

  if (!par().guess.empty()) {
    if (!envHasType(TGammaMap<V>, par().guess)) {
      HADRONS_ERROR(Argument,
                    "guess parameter '" + par().guess + "' must be a gamma "
                    "map with the same value structure as the source '" +
                    par().source + "'");
    }
    guessMap = &envGet(TGammaMap<V>, par().guess);
  }

  // Iterates the map's OWN (label, StagGamma) pairs directly -- never
  // re-derives the label from the gamma object, so a SpinTaste module's
  // custom `labels` override flows through unchanged.
  for (const auto &p : gammas) {
    const std::string &label = p.first;
    const StagGamma &gamma = p.second;
    LOG(Message) << "Solve for '" << par().source << "' with spin-taste: '"
                 << label << "'" << std::endl;

    auto solIt = solMap.find(label);
    if (solIt == solMap.end()) {
      HADRONS_ERROR(Argument, "internal: output map '" + getName() +
                                  "' lacks an entry for gamma '" + label +
                                  "'");
    }
    const V *guess = nullptr;
    if (guessMap != nullptr) {
      auto guessIt = guessMap->find(label);
      if (guessIt == guessMap->end()) {
        HADRONS_ERROR(Argument, "guess map '" + par().guess +
                                    "' has no entry for gamma '" + label +
                                    "'");
      }
      guess = &guessIt->second;
    }

    executeHelper(solIt->second, srcFor(gamma, label), gamma, guess);
  }
}

template <typename FImpl>
template <typename TField>
void TGaugePropMILC<FImpl>::executeHelper(const TField &src) {
  gammaLoop<TField>(
      [&src](const StagGamma &, const std::string &) -> const TField & {
        return src;
      });
}

template <typename FImpl>
template <typename TField>
void TGaugePropMILC<FImpl>::executeHelper(const std::vector<TField> &src) {
  gammaLoop<std::vector<TField>>(
      [&src](const StagGamma &,
             const std::string &) -> const std::vector<TField> & {
        return src;
      });
}

template <typename FImpl>
template <typename TField>
void TGaugePropMILC<FImpl>::executeHelper(const TGammaMap<TField> &srcMap) {
  gammaLoop<TField>(
      [this, &srcMap](const StagGamma &,
                      const std::string &label) -> const TField & {
        const std::string key =
            par().sourceLabel.empty() ? label : par().sourceLabel;
        auto it = srcMap.find(key);
        if (it == srcMap.end()) {
          HADRONS_ERROR(Argument, "source map '" + par().source +
                                      "' has no entry '" + key +
                                      "' for gamma '" + label + "'");
        }
        return it->second;
      });
}

template <typename FImpl>
template <typename TField>
void TGaugePropMILC<FImpl>::executeHelper(
    const TGammaMap<std::vector<TField>> &srcMap) {
  gammaLoop<std::vector<TField>>(
      [this, &srcMap](const StagGamma &,
                      const std::string &label) -> const std::vector<TField> & {
        const std::string key =
            par().sourceLabel.empty() ? label : par().sourceLabel;
        auto it = srcMap.find(key);
        if (it == srcMap.end()) {
          HADRONS_ERROR(Argument, "source map '" + par().source +
                                      "' has no entry '" + key +
                                      "' for gamma '" + label + "'");
        }
        return it->second;
      });
}

template <typename FImpl> void TGaugePropMILC<FImpl>::execute(void) {
  LOG(Message) << "Computing quark propagator '" << getName() << "'"
               << std::endl;

  if (envHasType(PropagatorField, par().source)) {
    executeHelper(envGet(PropagatorField, par().source));
  } else if (envHasType(std::vector<PropagatorField>, par().source)) {
    executeHelper(envGet(std::vector<PropagatorField>, par().source));
  } else if (envHasType(TGammaMap<PropagatorField>, par().source)) {
    executeHelper(envGet(TGammaMap<PropagatorField>, par().source));
  } else if (envHasType(TGammaMap<std::vector<PropagatorField>>,
                        par().source)) {
    executeHelper(
        envGet(TGammaMap<std::vector<PropagatorField>>, par().source));
  } else if (envHasType(FermionField, par().source)) {
    executeHelper(envGet(FermionField, par().source));
  } else if (envHasType(std::vector<FermionField>, par().source)) {
    executeHelper(envGet(std::vector<FermionField>, par().source));
  } else if (envHasType(TGammaMap<FermionField>, par().source)) {
    executeHelper(envGet(TGammaMap<FermionField>, par().source));
  } else if (envHasType(TGammaMap<std::vector<FermionField>>,
                        par().source)) {
    executeHelper(envGet(TGammaMap<std::vector<FermionField>>, par().source));
  } else {
    HADRONS_ERROR(Logic,
                  "Type of source '" + par().source + "' not recognized.");
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_GaugeProp_hpp_
