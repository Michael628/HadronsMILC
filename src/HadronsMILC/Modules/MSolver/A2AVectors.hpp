/*
 * A2AVectorsMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: fionnoh <fionnoh@gmail.com>
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
#ifndef HadronsMILC_MSolver_A2AVectors_hpp_
#define HadronsMILC_MSolver_A2AVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/A2AVectors.hpp>
#include "../../A2AVectors.hpp"
#include "../../DilutedNoise.hpp"
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Create all-to-all V & W vectors                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsMILCPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsMILCPar,
                                  std::string, noise,
                                  std::string, action,
                                  std::string, lowModes,
                                  std::string, solver,
                                  std::string, highOutput,
                                  bool,        highMultiFile);
};

template <typename FImpl, typename Pack>
class TA2AVectorsMILC : public Module<A2AVectorsMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsMILC<FImpl> A2A;

public:
    // constructor
    TA2AVectorsMILC(const std::string name);
    // destructor
    virtual ~TA2AVectorsMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);

    // setup
    virtual void setup(void);

    // execute
    virtual void execute(void);
private:
    unsigned int Nh_{0};
    bool noiseVector_{false};
    std::string subsolver_,solver_;
};

MODULE_REGISTER_TMP(StagA2AVectors, 
   ARG(TA2AVectorsMILC<STAGIMPL,MassShiftEigenPack<STAGIMPL> >), MSolver);

/******************************************************************************
 *                       TA2AVectorsMILC implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2AVectorsMILC<FImpl,Pack>::TA2AVectorsMILC(const std::string name)
: Module<A2AVectorsMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectorsMILC<FImpl,Pack>::getInput(void)
{
    std::vector<std::string> in {par().action};

    if (!par().noise.empty()) {
        in.push_back(par().noise);
    }
    solver_ = par().solver;
    subsolver_ = par().solver;
    if (!par().lowModes.empty()) {
        in.push_back(par().lowModes);
        if (env().hasObject(par().solver+"_subtract")) {
            in.push_back(par().solver+"_subtract");
            subsolver_ += "_subtract";
        } else {
            if(!par().solver.empty()) {
                in.push_back(par().solver);
            }
        }
    } else {
        if(!par().solver.empty()) {
            in.push_back(par().solver);
        }
    }
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectorsMILC<FImpl,Pack>::getOutput(void)
{
    std::vector<std::string> out = {};

    if (!par().noise.empty()) {
        out.push_back(getName() + "_v");
        out.push_back(getName() + "_w");
    }

    return out;
}

/******************************************************************************
 *              TA2AVectorsMILC setup                                         *
 ******************************************************************************/
template <typename FImpl, typename Pack>
void TA2AVectorsMILC<FImpl,Pack>::setup(void)
{
    auto   &action   = envGet(FMat, par().action);
    bool   lowGuess  = envHasType(std::vector<FermionField>,par().lowModes);
    int    Ls        = env().getObjectLs(par().action);


    if (lowGuess) {
        auto   &solver   = envGet(Solver, subsolver_);
        envTmp(A2A, "a2a", 1, action, solver);
    } else {
        auto   &solver   = envGet(Solver, solver_);
        envTmp(A2A, "a2a", 1, action, solver);
    }

    if (Ls > 1) {
       HADRONS_ERROR(Argument, "Ls > 1 not implemented");
    }

    envTmpLat(FermionField, "temp");

    if(envHasType(std::vector<FermionField>,par().noise)) {
        noiseVector_ = true;
    }

    if (!noiseVector_) {
        auto &noise = envGet(SpinColorDiagonalNoiseMILC<FImpl>, par().noise);

        Nh_ = noise.fermSize();

    } else {
        auto &noise = envGet(std::vector<FermionField>, par().noise);

        Nh_ = noise.size();
    }

    envCreate(std::vector<FermionField>, getName() + "_w", 1, 
          Nh_, envGetGrid(FermionField));
    auto &w     = envGet(std::vector<FermionField>, getName() + "_w");
    for (auto & vec: w) {
        vec = Zero();
    }

    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              Nh_, envGetGrid(FermionField));
        auto &v     = envGet(std::vector<FermionField>, getName() + "_v");
        for (auto & vec: v) {
            vec = Zero();
        }
}

/******************************************************************************
 *              TA2AVectorsMILC execution                                     *
 ******************************************************************************/
template <typename FImpl, typename Pack>
void TA2AVectorsMILC<FImpl,Pack>::execute(void)
{
    envGetTmp(A2A, a2a);
    auto        &action     = envGet(FMat, par().action);

    bool hasLow = !par().lowModes.empty();
    bool lowGuess = envHasType(std::vector<FermionField>,par().lowModes);

    if (hasLow)
    {
       LOG(Message) << "Computing all-to-all vectors "
                    << " using lowModes '" << par().lowModes << "' (low modes) and noise '"
                    << par().noise << "' (" << Nh_ 
                    << " noise vectors)" << std::endl;
    } else {
       LOG(Message) << "Computing all-to-all vectors "
                    << " using noise '" << par().noise << "' (" << Nh_ 
                    << " noise vectors)" << std::endl;
    }

    std::vector<FermionField> *noise_vec;
    SpinColorDiagonalNoiseMILC<FImpl> *noise;
    
    auto &v     = envGet(std::vector<FermionField>, getName() + "_v");
    auto &w     = envGet(std::vector<FermionField>, getName() + "_w");

    if (noiseVector_) {
        noise_vec     = &envGet(std::vector<FermionField>, par().noise);
    } else {
        noise = &envGet(SpinColorDiagonalNoiseMILC<FImpl>, par().noise);
    }


    startTimer("W high mode");
    for (int ih = 0; ih < Nh_; ih++) {
        if (noiseVector_) {
            w[ih] = noise_vec->at(ih);
        } else {
            w[ih] = noise->getFerm(ih);
        }
    }
    if (hasLow && !lowGuess) {
        LOG(Message) << "Projecting low contribution from stochastic high mode sources" << std::endl;

        auto &epack = envGet(Pack, par().lowModes);
        a2a.removeLowModeProjection(w,epack.evec,epack.eval);
    }
    stopTimer("W high mode");

    std::vector<FermionField> *guess;
    if (lowGuess) {
        guess = &envGet(std::vector<FermionField>, par().lowModes);
    }

    startTimer("V high mode");
    for (int ih = 0; ih < Nh_; ih++)
    {
        LOG(Message) << "V vector (solve) i = " << ih
                     << " (" << ((hasLow) ? "high " : "") 
                     << "stochastic mode)" << std::endl;

        if (lowGuess) {
            const FermionField &g = guess->at(ih);
            a2a.makeHighModeV(v[ih],w[ih],g);
            envGetTmp(FermionField,temp);
            action.M(g,temp);
            w[ih] = w[ih] - temp;
            RealD n_temp = 1.0/::sqrt(norm2(v[ih])*norm2(g));
            LOG(Message) << "check orthog <v,guess>: " << TensorRemove(innerProduct(v[ih],g))*n_temp << std::endl;
            n_temp = 1.0/::sqrt(norm2(temp)*norm2(w[ih]));
            LOG(Message) << "check orthog: <w,proj>:  " << TensorRemove(innerProduct(w[ih],temp))*n_temp << std::endl;
        }
        else {
            a2a.makeHighModeV(v[ih],w[ih]);
        }

    }
   stopTimer("V high mode");

    if (!par().highOutput.empty())
    {
       startTimer("V I/O");
       A2AVectorsIo::write(par().highOutput + "_v", v, par().highMultiFile, vm().getTrajectory());
       stopTimer("V I/O");
       startTimer("W I/O");
       A2AVectorsIo::write(par().highOutput + "_w", w, par().highMultiFile, vm().getTrajectory());
       stopTimer("W I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectorsMILC_hpp_
