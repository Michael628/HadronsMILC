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
#include <spin/StagGamma.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                MesonMILC                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class MesonMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonMILCPar, std::string, source,
                                  std::string, sink, std::string, sourceGammas,
                                  SpinTasteParams, sinkSpinTaste, std::string,
                                  sinkFunc, std::string, sourceShift,
                                  std::string, output);
};

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
           StagGamma &gamma);
  template <typename TField>
  EnableIf<is_lattice<TField>, void>
  contract(Result &result, const std::vector<TField> &source,
           const std::vector<TField> &sink, StagGamma &gamma);
  template <typename TField>
  void executeHelper(std::vector<Result> &results, const TField &sink);

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
  void parseGammas(void);

  std::string _sinkSuffix = "";
  std::vector<std::string> _sourceGammas;
  std::map<std::string, StagGamma::SpinTastePair> _mapSinkGammas;
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

template <typename FImpl> void TMesonMILC<FImpl>::parseGammas(void) {

  _sourceGammas.clear();
  _mapSinkGammas.clear();

  if (!par().sourceGammas.empty()) {
    for (auto gamma : strToVec<StagGamma::SpinTastePair>(par().sourceGammas)) {
      _sourceGammas.push_back(StagGamma::GetName(gamma));
    }
  }
  if (!par().sinkSpinTaste.gammas.empty()) {
    auto gamma_vals = StagGamma::ParseSpinTasteString(
        par().sinkSpinTaste.gammas, par().sinkSpinTaste.applyG5);
    auto gamma_keys =
        StagGamma::ParseSpinTasteString(par().sinkSpinTaste.gammas);
    for (int i = 0; i < gamma_vals.size(); ++i) {
      _mapSinkGammas.insert({StagGamma::GetName(gamma_keys[i]), gamma_vals[i]});
    }
  }
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMesonMILC<FImpl>::getInput(void) {
  parseGammas();
  std::vector<std::string> in = {par().sinkFunc};

  if (!par().sourceShift.empty()) {
    in.push_back(par().sourceShift);
  }

  if (!par().sinkSpinTaste.gauge.empty()) {
    in.push_back(par().sinkSpinTaste.gauge);
  }

  for (auto iter = _sourceGammas.begin(); iter != _sourceGammas.end(); ++iter) {
    in.push_back(par().source + *iter);
  }

  if (_sourceGammas.empty())
    in.push_back(par().source);

  // std::string identityName =
  // StagGamma::GetName(StagGamma::StagAlgebra::G5,StagGamma::StagAlgebra::G5);

  // if (env().hasObject(par().sink + identityName)) {
  // _sinkSuffix = identityName;
  // in.push_back(par().sink+identityName);
  // } else {
  in.push_back(par().sink);
  // }

  return in;
}

template <typename FImpl>
std::vector<std::string> TMesonMILC<FImpl>::getOutput(void) {
  std::vector<std::string> out = {};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl> void TMesonMILC<FImpl>::setup(void) {

  envTmpLat(PropagatorField, "prop");

  _Nt = env().getDim(Tp);

  if (_sourceGammas.size() > 0 &&
      _sourceGammas.size() != _mapSinkGammas.size()) {
    HADRONS_ERROR(Argument,
                  "Parameter 'sourceGammas' must be empty or have the same "
                  "number of operators as 'sinkSpinTaste.gammas'.");
  }
  if (envHasType(PropagatorField, par().sink + _sinkSuffix) ||
      envHasType(std::vector<PropagatorField>, par().sink + _sinkSuffix)) {
    envTmpLat(PropagatorField, "field");
  } else if (envHasType(FermionField, par().sink + _sinkSuffix) ||
             envHasType(std::vector<FermionField>, par().sink + _sinkSuffix)) {
    envTmpLat(FermionField, "field");
  } else {
    HADRONS_ERROR(Argument, "Sink parameter '" + par().sink +
                                "' must be a PropagatorField, FermionField, or "
                                "a std::vector of these fields.");
  }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename TField>
EnableIf<is_lattice<TField>, void>
TMesonMILC<FImpl>::contract(Result &result, const TField &source,
                            const TField &sink, StagGamma &gamma) {

  int offset;
  std::vector<TComplex> buf;

  SinkFnScalar &sinkFunc = envGet(SinkFnScalar, par().sinkFunc);

  envGetTmp(PropagatorField, prop);
  envGetTmp(TField, field);

  gamma(field, source);

  buildProp(prop, field, sink);

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
                            const std::vector<TField> &sink, StagGamma &gamma) {

  int offset;
  std::vector<TComplex> buf;

  SinkFnScalar &sinkFunc = envGet(SinkFnScalar, par().sinkFunc);

  result.srcCorrs.resize(sink.size(), std::vector<Complex>(_Nt, 0.0));
  result.scaling = sink.size();

  envGetTmp(PropagatorField, prop);
  envGetTmp(TField, field);

  for (int i = 0; i < result.srcCorrs.size(); i++) {

    LOG(Message) << "Contracting element i = " << i << "." << std::endl;

    gamma(field, source[i]);
    buildProp(prop, field, sink[i]);

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
template <typename TField>
void TMesonMILC<FImpl>::executeHelper(std::vector<Result> &results,
                                      const TField &sink) {

  std::string srcName;
  StagGamma spinTaste;

  if (!par().sinkSpinTaste.gauge.empty()) {
    auto &U = envGet(GaugeField, par().sinkSpinTaste.gauge);
    spinTaste.setGaugeField(U);
  }

  int i = 0;
  for (auto iter = _mapSinkGammas.begin(); iter != _mapSinkGammas.end();
       ++iter) {
    spinTaste.setSpinTaste(iter->second);

    LOG(Message) << "Contracting with gamma: " << iter->first << std::endl;

    srcName = par().source;

    if (!par().sourceGammas.empty()) {
      srcName += _sourceGammas[i];
      results[i].sourceGamma = _sourceGammas[i];
      LOG(Message) << "Using source gamma: '" << results[i].sourceGamma << "'."
                   << std::endl;
    } else {
      results[i].sourceGamma = "N/A";
    }

    auto &source = envGet(TField, srcName);

    contract(results[i], source, sink, spinTaste);

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
               << " quarks '" << par().source << "' and '" << par().sink << "'"
               << std::endl;

  std::vector<Result> results;

  results.resize(_mapSinkGammas.size());

  int i = 0;
  for (auto iter = _mapSinkGammas.begin(); iter != _mapSinkGammas.end();
       ++iter) {
    results[i].sinkGamma = iter->first;
    results[i].srcCorrs.resize(1, std::vector<Complex>(_Nt, 0.0));
    results[i].corr.resize(_Nt, 0.0);
    results[i].scaling = 1.0;
    if (!par().sourceShift.empty()) {
      results[i].timeShifts = envGet(std::vector<Integer>, par().sourceShift);
    }

    i++;
  }

  if (envHasType(PropagatorField, par().sink + _sinkSuffix)) {
    auto &sink = envGet(PropagatorField, par().sink + _sinkSuffix);
    executeHelper(results, sink);

  } else if (envHasType(std::vector<PropagatorField>,
                        par().sink + _sinkSuffix)) {
    auto &sink = envGet(std::vector<PropagatorField>, par().sink + _sinkSuffix);
    executeHelper(results, sink);

  } else if (envHasType(FermionField, par().sink + _sinkSuffix)) {
    auto &sink = envGet(FermionField, par().sink + _sinkSuffix);
    executeHelper(results, sink);

  } else if (envHasType(std::vector<FermionField>, par().sink + _sinkSuffix)) {
    auto &sink = envGet(std::vector<FermionField>, par().sink + _sinkSuffix);
    executeHelper(results, sink);
  }

  saveResult(par().output, "meson", results);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Meson_hpp_
