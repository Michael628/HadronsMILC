/*
 * MesonMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
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

#ifndef HadronsMILC_MContraction_Meson_hpp_
#define HadronsMILC_MContraction_Meson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <GridMilc/GridMilc.h>

#include <Modules/MFermion/SpinTaste.hpp> // TGammaMap + gammas module

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                MesonMILC                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class MesonMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonMILCPar, std::string, source,
                                  std::string, sink, std::string, sinkGammas,
                                  std::string, sinkFunc, std::string,
                                  sourceShift, std::string, output);
};
// sinkGammas: REQUIRED name of an MFermion::SpinTaste module; this module
//             consumes the SpinTaste module's `_map` companion output
//             (par().sinkGammas + "_map", a TGammaMap<StagGamma>) as the
//             LOOP DRIVER and single source of truth for key matching.
//             Empty is a setup-time error (StagMesonLegacy preserves the
//             former sourceGammas/sinkSpinTaste workflow).
// source/sink: each independently may be a bare TField/std::vector<TField>
//             (used UNCHANGED for every gamma) or a TGammaMap<TField>/
//             TGammaMap<std::vector<TField>> (per-gamma lookup by label,
//             every sinkGammas key validated present at setup)

template <typename FImpl> class TMesonMILC : public Module<MesonMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );
  BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
  SINK_TYPE_ALIASES(Scalar);

  class Result : Serializable {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Result, std::string, sourceGamma,
                                    std::string, sinkGamma,
                                    std::vector<Complex>, corr,
                                    std::vector<std::vector<Complex>>, srcCorrs,
                                    std::vector<Integer>, timeShifts, Real,
                                    scaling);
  };

public:
  // constructor
  TMesonMILC(const std::string name);
  // destructor
  virtual ~TMesonMILC(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);

protected:
  template <typename TField>
  EnableIf<is_lattice<TField>, void>
  contract(Result &result, const TField &source, const TField &sink,
           const StagGamma &gamma);
  template <typename TField>
  EnableIf<is_lattice<TField>, void>
  contract(Result &result, const std::vector<TField> &source,
           const std::vector<TField> &sink, const StagGamma &gamma);

  // Generalizes GaugeProp's gammaLoop/SrcFetch-lambda pattern
  // (GaugeProp.hpp) to Meson's symmetric source+sink dispatch: bare V ->
  // fixed object for every gamma; TGammaMap<V> -> per-label lookup
  // (validated present by checkKeys at setup).
  template <typename V>
  std::function<const V &(const std::string &)> fetchFor(const std::string &name);
  // Setup-time key-set validation (GammaMapElement.hpp precedent): if
  // `name` names a TGammaMap<V>, every sinkGammas key must be present
  // (fatal otherwise, dynamic available-keys diagnostic); bare objects
  // need no check (used unchanged for every gamma).
  template <typename V>
  void checkKeys(const std::string &name,
                const TGammaMap<StagGamma> &sinkGammas) const;
  template <typename V>
  void executeHelper(std::vector<Result> &results,
                     const TGammaMap<StagGamma> &sinkGammas);

  inline void buildProp(PropagatorField &result, const FermionField &source,
                        const FermionField &sink) {
    result = outerProduct(sink, source);
  }
  inline void buildProp(PropagatorField &result, const PropagatorField &source,
                        const PropagatorField &sink) {
    result = sink * adj(source);
  }

  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);

private:
  std::string _sinkSuffix = "";
  Integer _Nt;
};

MODULE_REGISTER_TMP(StagMeson, TMesonMILC<STAGIMPL>, MContraction);

/******************************************************************************
 *                       TMesonMILC implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMesonMILC<FImpl>::TMesonMILC(const std::string name)
    : Module<MesonMILCPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMesonMILC<FImpl>::getInput(void) {
  std::vector<std::string> in = {par().sinkFunc, par().source, par().sink,
                                 par().sinkGammas + "_map"};

  if (!par().sourceShift.empty()) {
    in.push_back(par().sourceShift);
  }

  return in;
}

template <typename FImpl>
std::vector<std::string> TMesonMILC<FImpl>::getOutput(void) {
  std::vector<std::string> out = {};

  return out;
}

// setup-time key-set validation //////////////////////////////////////////////
template <typename FImpl>
template <typename V>
void TMesonMILC<FImpl>::checkKeys(const std::string &name,
                                 const TGammaMap<StagGamma> &sinkGammas) const {
  if (!envHasType(TGammaMap<V>, name)) {
    return; // bare object: used unchanged for every gamma, no key-set check
  }
  const auto &map = envGet(TGammaMap<V>, name);
  for (const auto &p : sinkGammas) {
    if (map.find(p.first) == map.end()) {
      std::string available;
      for (auto &q : map) {
        available += (available.empty() ? "" : ", ") + q.first;
      }
      HADRONS_ERROR(Argument, "map '" + name + "' has no entry for gamma '" +
                                  p.first + "' -- sinkGammas module '" +
                                  par().sinkGammas + "' and '" + name +
                                  "'s gamma list disagree (available: " +
                                  available + ")");
    }
  }
}

// per-gamma fetch dispatch ////////////////////////////////////////////////////
template <typename FImpl>
template <typename V>
std::function<const V &(const std::string &)>
TMesonMILC<FImpl>::fetchFor(const std::string &name) {
  if (envHasType(TGammaMap<V>, name)) {
    auto &map = envGet(TGammaMap<V>, name);
    return [&map](const std::string &label) -> const V & {
      return map.at(label);
    };
  }
  auto &obj = envGet(V, name);
  return [&obj](const std::string &) -> const V & { return obj; };
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl> void TMesonMILC<FImpl>::setup(void) {
  envTmpLat(PropagatorField, "prop");
  _Nt = env().getDim(Tp);

  if (par().sinkGammas.empty()) {
    HADRONS_ERROR(Argument,
                  "Meson requires the 'sinkGammas' SpinTaste module name "
                  "(use StagMesonLegacy for the former sourceGammas/"
                  "sinkSpinTaste workflow)");
  }
  const auto &sinkGammas =
      envGet(TGammaMap<StagGamma>, par().sinkGammas + "_map");
  if (sinkGammas.empty()) {
    HADRONS_ERROR(Argument,
                  "SpinTaste module '" + par().sinkGammas + "' published "
                  "an empty gamma map -- Meson requires at least one "
                  "gamma");
  }

  if (envHasType(PropagatorField, par().sink + _sinkSuffix) ||
      envHasType(TGammaMap<PropagatorField>, par().sink + _sinkSuffix)) {
    envTmpLat(PropagatorField, "field");
    if (!envHasType(PropagatorField, par().source) &&
        !envHasType(TGammaMap<PropagatorField>, par().source)) {
      HADRONS_ERROR(Argument, "source parameter '" + par().source +
                                  "' must have the same field-and-container "
                                  "structure (PropagatorField, bare or "
                                  "TGammaMap thereof) as sink parameter '" +
                                  par().sink + _sinkSuffix + "'");
    }
    checkKeys<PropagatorField>(par().source, sinkGammas);
    checkKeys<PropagatorField>(par().sink + _sinkSuffix, sinkGammas);
  } else if (envHasType(std::vector<PropagatorField>,
                        par().sink + _sinkSuffix) ||
             envHasType(TGammaMap<std::vector<PropagatorField>>,
                        par().sink + _sinkSuffix)) {
    envTmpLat(PropagatorField, "field");
    if (!envHasType(std::vector<PropagatorField>, par().source) &&
        !envHasType(TGammaMap<std::vector<PropagatorField>>, par().source)) {
      HADRONS_ERROR(Argument, "source parameter '" + par().source +
                                  "' must have the same field-and-container "
                                  "structure (std::vector<PropagatorField>, "
                                  "bare or TGammaMap thereof) as sink "
                                  "parameter '" + par().sink + _sinkSuffix +
                                  "'");
    }
    checkKeys<std::vector<PropagatorField>>(par().source, sinkGammas);
    checkKeys<std::vector<PropagatorField>>(par().sink + _sinkSuffix,
                                            sinkGammas);
  } else if (envHasType(FermionField, par().sink + _sinkSuffix) ||
             envHasType(TGammaMap<FermionField>, par().sink + _sinkSuffix)) {
    envTmpLat(FermionField, "field");
    if (!envHasType(FermionField, par().source) &&
        !envHasType(TGammaMap<FermionField>, par().source)) {
      HADRONS_ERROR(Argument, "source parameter '" + par().source +
                                  "' must have the same field-and-container "
                                  "structure (FermionField, bare or "
                                  "TGammaMap thereof) as sink parameter '" +
                                  par().sink + _sinkSuffix + "'");
    }
    checkKeys<FermionField>(par().source, sinkGammas);
    checkKeys<FermionField>(par().sink + _sinkSuffix, sinkGammas);
  } else if (envHasType(std::vector<FermionField>,
                        par().sink + _sinkSuffix) ||
             envHasType(TGammaMap<std::vector<FermionField>>,
                        par().sink + _sinkSuffix)) {
    envTmpLat(FermionField, "field");
    if (!envHasType(std::vector<FermionField>, par().source) &&
        !envHasType(TGammaMap<std::vector<FermionField>>, par().source)) {
      HADRONS_ERROR(Argument, "source parameter '" + par().source +
                                  "' must have the same field-and-container "
                                  "structure (std::vector<FermionField>, "
                                  "bare or TGammaMap thereof) as sink "
                                  "parameter '" + par().sink + _sinkSuffix +
                                  "'");
    }
    checkKeys<std::vector<FermionField>>(par().source, sinkGammas);
    checkKeys<std::vector<FermionField>>(par().sink + _sinkSuffix,
                                         sinkGammas);
  } else {
    HADRONS_ERROR(Argument, "Sink parameter '" + par().sink +
                                "' must be a PropagatorField, FermionField, "
                                "a TGammaMap thereof, or a std::vector of "
                                "these fields.");
  }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename TField>
EnableIf<is_lattice<TField>, void>
TMesonMILC<FImpl>::contract(Result &result, const TField &source,
                            const TField &sink, const StagGamma &gamma) {

  int offset;
  std::vector<TComplex> buf;

  SinkFnScalar &sinkFunc = envGet(SinkFnScalar, par().sinkFunc);

  envGetTmp(PropagatorField, prop);
  envGetTmp(TField, field);

  gamma(field, sink);

  buildProp(prop, field, source);

  buf = sinkFunc(trace(prop));

  offset = 0;
  if (result.timeShifts.size() > 0) {
    LOG(Message) << "Shifting correlator to (t0 = " << result.timeShifts[0]
                 << ")" << std::endl;
    offset = result.timeShifts[0];
  }

  for (unsigned int t = 0; t < _Nt; ++t) {
    auto ct = TensorRemove(buf[offset]);
    result.srcCorrs[0][t] = ct; // Save corr for individual source

    offset = mod(offset + 1, _Nt);
  }
}

template <typename FImpl>
template <typename TField>
EnableIf<is_lattice<TField>, void>
TMesonMILC<FImpl>::contract(Result &result, const std::vector<TField> &source,
                            const std::vector<TField> &sink, const StagGamma &gamma) {

  int offset;
  std::vector<TComplex> buf;

  SinkFnScalar &sinkFunc = envGet(SinkFnScalar, par().sinkFunc);

  result.srcCorrs.resize(sink.size(), std::vector<Complex>(_Nt, 0.0));
  result.scaling = sink.size();

  envGetTmp(PropagatorField, prop);
  envGetTmp(TField, field);

  for (int i = 0; i < result.srcCorrs.size(); i++) {

    LOG(Message) << "Contracting element i = " << i << "." << std::endl;

    gamma(field, sink[i]);
    buildProp(prop, field, source[i]);

    buf = sinkFunc(trace(prop));

    offset = 0;
    if (result.timeShifts.size() > 0) {
      LOG(Message) << "Shifting correlator " << i
                   << " to (t0 = " << result.timeShifts[i] << ")" << std::endl;
      offset = result.timeShifts[i];
    }

    for (unsigned int t = 0; t < _Nt; ++t) {
      auto ct = TensorRemove(buf[offset]);
      result.srcCorrs[i][t] = ct; // Save corr for individual source

      offset = mod(offset + 1, _Nt);
    }
  }
}

template <typename FImpl>
template <typename V>
void TMesonMILC<FImpl>::executeHelper(std::vector<Result> &results,
                                     const TGammaMap<StagGamma> &sinkGammas) {
  auto sourceFetch = fetchFor<V>(par().source);
  auto sinkFetch = fetchFor<V>(par().sink + _sinkSuffix);

  int i = 0;
  for (const auto &p : sinkGammas) {
    const std::string &label = p.first;
    const StagGamma &gamma = p.second;

    results[i].sourceGamma = label;
    results[i].sinkGamma = label;
    LOG(Message) << "Contracting with gamma: " << label << std::endl;

    contract(results[i], sourceFetch(label), sinkFetch(label), gamma);

    for (int j = 0; j < results[i].srcCorrs.size(); j++) {
      for (int t = 0; t < _Nt; t++) {
        results[i].corr[t] += (results[i].srcCorrs[j])[t] / results[i].scaling;
      }
    }

    i++;
  }
}

template <typename FImpl> void TMesonMILC<FImpl>::execute(void) {
  LOG(Message) << "Computing meson contractions '" << getName() << "' using"
               << " quarks '" << par().source << "' and '" << par().sink
               << "'" << std::endl;

  const auto &sinkGammas =
      envGet(TGammaMap<StagGamma>, par().sinkGammas + "_map");

  std::vector<Result> results(sinkGammas.size());
  int i = 0;
  for (const auto &p : sinkGammas) {
    results[i].srcCorrs.resize(1, std::vector<Complex>(_Nt, 0.0));
    results[i].corr.resize(_Nt, 0.0);
    results[i].scaling = 1.0;
    if (!par().sourceShift.empty()) {
      results[i].timeShifts = envGet(std::vector<Integer>, par().sourceShift);
    }
    i++;
  }

  if (envHasType(PropagatorField, par().sink + _sinkSuffix) ||
      envHasType(TGammaMap<PropagatorField>, par().sink + _sinkSuffix)) {
    executeHelper<PropagatorField>(results, sinkGammas);
  } else if (envHasType(std::vector<PropagatorField>,
                        par().sink + _sinkSuffix) ||
             envHasType(TGammaMap<std::vector<PropagatorField>>,
                        par().sink + _sinkSuffix)) {
    executeHelper<std::vector<PropagatorField>>(results, sinkGammas);
  } else if (envHasType(FermionField, par().sink + _sinkSuffix) ||
             envHasType(TGammaMap<FermionField>, par().sink + _sinkSuffix)) {
    executeHelper<FermionField>(results, sinkGammas);
  } else if (envHasType(std::vector<FermionField>,
                        par().sink + _sinkSuffix) ||
             envHasType(TGammaMap<std::vector<FermionField>>,
                        par().sink + _sinkSuffix)) {
    executeHelper<std::vector<FermionField>>(results, sinkGammas);
  }
  // no `else`: setup() already validated one of the four branches above
  // matches, per GammaMapElement's "full enumeration in setup(), silent
  // re-dispatch in execute()" precedent -- this closes the pre-existing
  // gap where execute()'s ladder had no else AND no map branches.

  saveResult(par().output, "meson", results);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Meson_hpp_
