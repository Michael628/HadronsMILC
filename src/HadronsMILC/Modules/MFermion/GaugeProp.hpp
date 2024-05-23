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
#include "../../../spin/StagGamma.h"

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                GaugePropMILC                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class GaugePropMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugePropMILCPar,
                                    std::string, source,
                                    SpinTasteParams, spinTaste,
                                    std::string, solver,
                                    std::string, guess);
};

template <typename FImpl>
class TGaugePropMILC: public Module<GaugePropMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);

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
    void setupHelper(void);
    virtual void setup(void);
    // execution
    template <typename TField>
    EnableIf<is_lattice<TField>,void> executeHelper(TField &sol, const TField &src, StagGamma& gamma, const TField *guess = nullptr);
    template <typename TField>
    EnableIf<is_lattice<TField>,void> executeHelper(std::vector<TField> &sol, const std::vector<TField> &src, StagGamma& gamma, const std::vector<TField> *guess = nullptr);
    template <typename TField>
    void executeHelper(const TField &src);
    virtual void execute(void);
private:
    void parseGammas(void);
    void solveField(FermionField &prop, const FermionField &src, const FermionField* guess = nullptr);
    void solveField(PropagatorField &prop, const PropagatorField &src, const PropagatorField* guess = nullptr);
private:
    bool _hasGuess;
    std::map<std::string,StagGamma::SpinTastePair> _mapGammas;
};

MODULE_REGISTER_TMP(StagGaugeProp, TGaugePropMILC<STAGIMPL>, MFermion);

/******************************************************************************
 *                      TGaugePropMILC implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaugePropMILC<FImpl>::TGaugePropMILC(const std::string name)
: Module<GaugePropMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
void TGaugePropMILC<FImpl>::parseGammas(void) {
    _mapGammas.clear();
    if(! par().spinTaste.gammas.empty()) {
        auto gamma_vals = StagGamma::ParseSpinTasteString(par().spinTaste.gammas,par().spinTaste.applyG5);
        auto gamma_keys = StagGamma::ParseSpinTasteString(par().spinTaste.gammas);
        for (int i = 0; i < gamma_vals.size(); ++i)
        {
            _mapGammas.insert({StagGamma::GetName(gamma_keys[i]),gamma_vals[i]});
        }
    }
}
template <typename FImpl>
std::vector<std::string> TGaugePropMILC<FImpl>::getInput(void)
{

    parseGammas();

    std::vector<std::string> in = {par().source, par().solver};

    if (!par().spinTaste.gauge.empty()) {
        in.push_back(par().spinTaste.gauge);
    }

    if (!par().guess.empty() && !_mapGammas.empty()) {
        for (auto iter = _mapGammas.begin(); iter != _mapGammas.end(); ++iter)
        {
            in.push_back(par().guess+iter->first);
        }
    } else if (!par().guess.empty()) {
        in.push_back(par().guess);
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TGaugePropMILC<FImpl>::getOutput(void)
{
    parseGammas();
    std::vector<std::string> out;

    if (!_mapGammas.empty()) {
        for (auto iter = _mapGammas.begin(); iter != _mapGammas.end(); ++iter)
        {
            out.push_back(getName()+iter->first);
        }
    } else {
        out.push_back(getName());
    }
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename TField>
void TGaugePropMILC<FImpl>::setupHelper() {
    envTmpLat(TField, "field");

    // Create an output field for each gamma
    auto initializeOutput = [this](std::string ext){
        envCreate(TField, getName()+ext, 1, envGetGrid(TField));
        auto &sol = envGet(TField, getName()+ext);
        sol = Zero();
    };

    // Create an output field for each gamma
    auto initializeVectorOutput = [this](std::string ext){
        auto &src = envGet(std::vector<TField>, par().source);
        envCreate(std::vector<TField>, getName()+ext, 1, std::vector<TField>(src.size(),envGetGrid(TField)));

        auto &sol = envGet(std::vector<TField>, getName()+ext);
        for (auto & s:sol) {
            s = Zero();
        }
    };


    for (auto iter = _mapGammas.begin(); iter != _mapGammas.end(); ++iter)
    {
        if (envHasType(TField,par().source)) {
            initializeOutput(iter->first);
        } else {
            initializeVectorOutput(iter->first);
        }
    }
    if (_mapGammas.empty()) {
        if (envHasType(TField,par().source)) {
            initializeOutput("");
        } else {
            initializeVectorOutput("");
        }
    }
}

template <typename FImpl>
void TGaugePropMILC<FImpl>::setup(void)
{
    _hasGuess = !par().guess.empty();

    if (envHasType(PropagatorField,par().source) || envHasType(std::vector<PropagatorField>,par().source)) {

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

       setupHelper<PropagatorField>();

    } else if (envHasType(FermionField,par().source) || envHasType(std::vector<FermionField>,par().source)) {
        setupHelper<FermionField>();
    } else {
        HADRONS_ERROR(Logic,"Type of source '" + par().source + "' not recognized.");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugePropMILC<FImpl>::solveField(FermionField &sol, const FermionField &src, const FermionField *guess)
{
    auto &solver  = envGet(Solver, par().solver);
    
    if (guess != nullptr) {
        solver(sol, src,*guess);
    } else {
        solver(sol, src);
    }
}

template <typename FImpl>
void TGaugePropMILC<FImpl>::solveField(PropagatorField &sol, 
                                            const PropagatorField &src, const PropagatorField *guess)
{
    auto &solver  = envGet(Solver, par().solver);
    
    envGetTmp(FermionField, fermIn);
    envGetTmp(FermionField, fermOut);
    envGetTmp(FermionField, fermGuess);

    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        PropToFerm<FImpl>(fermIn, src, c);
        if (guess != nullptr) {
            PropToFerm<FImpl>(fermGuess,*guess,c);
            solver(fermOut, fermIn,fermGuess);
        } else {
            solver(fermOut, fermIn);
        }
        FermToProp<FImpl>(sol, fermOut, c);
    }
}

template <typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TGaugePropMILC<FImpl>::executeHelper(TField &sol, const TField &src, StagGamma& gamma, const TField *guess)
{
    envGetTmp(TField,field);

    const TField* temp;

    if (_mapGammas.empty()) {
        temp = &src;
    } else {
        gamma(field,src);
        temp = &field;
    }

    if (_hasGuess) {
        solveField(sol,*temp,guess);
    } else {
        solveField(sol,*temp);
    }
}

template <typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TGaugePropMILC<FImpl>::executeHelper(std::vector<TField> &sol, const std::vector<TField> &src, StagGamma& gamma, const std::vector<TField> *guess)
{
    envGetTmp(TField,field);

    for (int i = 0; i < src.size(); i++) {

        const TField* temp;
        if (_mapGammas.empty()) {
            temp = &src[i];
        } else {
            gamma(field,src[i]);
            temp = &field;
        }

        if (_hasGuess) {
            const TField* guessTemp = &(guess->at(i));
            solveField(sol[i],*temp,guessTemp);
        } else {
            solveField(sol[i],*temp);
        }
    }
}

template <typename FImpl>
template<typename TField>
void TGaugePropMILC<FImpl>::executeHelper(const TField &src)
{
    StagGamma gamma;

    auto solveFunc = [this, &src, &gamma](std::string ext) {
        auto& sol = envGet(TField,getName()+ext);

        if (!par().guess.empty()) {
            if (!envHasType(TField,par().guess+ext)) {
               HADRONS_ERROR(Argument, "guess parameter '" + par().guess + ext + "' must have same data structure as source, '"+par().source+"'");
            }
            auto & guess = envGet(TField,par().guess+ext);
            executeHelper(sol,src,gamma,&guess);
        } else {
            executeHelper(sol,src,gamma);
        }

    };

    if (!(_mapGammas.empty() || par().spinTaste.gauge.empty())) {
        auto& Umu = envGet(GaugeField,par().spinTaste.gauge);
        gamma.setGaugeField(Umu);
    }

    for (auto iter = _mapGammas.begin(); iter != _mapGammas.end(); ++iter)
    {
        // Apply gamma to source
        gamma.setSpinTaste(iter->second);
        LOG(Message) << "Solve for '" << par().source << "' with spin-taste: '" << iter->first << "'" << std::endl;

        solveFunc(iter->first);
    }

    if (_mapGammas.empty()) {
        gamma.setSpinTaste(StagGamma::StagAlgebra::G1,StagGamma::StagAlgebra::G1);
        solveFunc("");
    }
}

template <typename FImpl>
void TGaugePropMILC<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    
    if (envHasType(PropagatorField,par().source)) {
        const auto &src = envGet(PropagatorField,par().source);
        executeHelper(src);
    } else if (envHasType(std::vector<PropagatorField>, par().source)) {
        const auto &src = envGet(std::vector<PropagatorField>, par().source);
        executeHelper(src);
    } else if (envHasType(FermionField,par().source)) {
        const auto &src = envGet(FermionField,par().source);
        executeHelper(src);
    } else {
        const auto &src = envGet(std::vector<FermionField>, par().source);
        executeHelper(src);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_GaugeProp_hpp_
