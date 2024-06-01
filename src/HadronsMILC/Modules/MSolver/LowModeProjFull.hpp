/*
 * LowModeProjMILCFull.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MSolver_LowModeProj_hpp_
#define HadronsMILC_MSolver_LowModeProj_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Calculate Low Mode Average Prop                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class LowModeProjMILCFullPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LowModeProjMILCFullPar,
                                  std::string, action,
                                  bool, projector,
                                  std::string, lowModes);
};

template <typename FImpl, typename Pack>
class TLowModeProjMILCFull : public Module<LowModeProjMILCFullPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLowModeProjMILCFull(const std::string name);
    // destructor
    virtual ~TLowModeProjMILCFull(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);

    // setup
    virtual void setup(void);

    // execute
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StagLMAFull, ARG(TLowModeProjMILCFull<STAGIMPL,MassShiftEigenPack<STAGIMPL> >), MSolver);

/******************************************************************************
 *                       TLowModeProjMILCFull implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TLowModeProjMILCFull<FImpl,Pack>::TLowModeProjMILCFull(const std::string name)
: Module<LowModeProjMILCFullPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TLowModeProjMILCFull<FImpl,Pack>::getInput(void)
{
    std::vector<std::string> in {par().action,par().lowModes};
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLowModeProjMILCFull<FImpl,Pack>::getOutput(void)
{
    std::vector<std::string> out {getName(),getName()+"_subtract"};
    return out;
}

template <typename FImpl, typename Pack>
DependencyMap TLowModeProjMILCFull<FImpl,Pack>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().action, getName()});
    dep.insert({par().lowModes, getName()});
    dep.insert({par().action, getName()+"_subtract"});
    dep.insert({par().lowModes, getName()+"_subtract"});

    return dep;
}

/******************************************************************************
 *              TLowModeProjMILCFull setup                                         *
 ******************************************************************************/

// C++11 does not support template lambdas so it is easier
// to make a macro with the solver body
#define SOLVER_BODY                                                                 \
    auto &rbTemp = envGet(FermionField,"rbTemp");                                   \
    auto &temp = envGet(FermionField,"temp");                             \
    auto &tempDag = envGet(FermionField,"tempDag");                             \
    auto &rbFerm = envGet(FermionField,"rbFerm");                                   \
    auto &rbFermNeg = envGet(FermionField,"rbFermNeg");                             \
    auto &MrbFermNeg = envGet(FermionField,"MrbFermNeg");                           \
                                                \
    int cb = epack.evec[0].Checkerboard();                                            \
    int cbNeg = (cb==Even) ? Odd : Even;                                              \
                                                \
    RealD norm = 1./::sqrt(norm2(epack.evec[0]));                                     \
                                                \
    sol = Zero(); \
    rbTemp = Zero();                                                                  \
    rbTemp.Checkerboard() = cb;                                                       \
    temp = Zero();                                                               \
    tempDag = Zero();                                                               \
    for (int k=epack.evec.size()-1;k >= 0;k--) {                                      \
        const FermionField& e = epack.evec[k];                                        \
        const ComplexD lam_D    = ComplexD(0.0,epack.eval[k].imag());                 \
        const ComplexD lam     = epack.eval[k];                                  \
        mat.Meooe(e, rbTemp); \
        rbTemp = rbTemp/lam_D; \
        setCheckerboard(temp,e); \
        setCheckerboard(temp,rbTemp); \
        if (cb == Even) { \
            setCheckerboard(tempDag,e); \
            setCheckerboard(tempDag,-rbTemp); \
        } else { \
            setCheckerboard(tempDag,-e); \
            setCheckerboard(tempDag,rbTemp); \
        } \
        const ComplexD ip    = TensorRemove(innerProduct(temp,source))/lam;       \
        const ComplexD ipDag = TensorRemove(innerProduct(tempDag,source))/conjugate(lam);   \
        axpy(sol,ip,temp,sol); \
        axpy(sol,ipDag,tempDag,sol); \
    }                                                                                 \
    sol *= norm;                         

template <typename FImpl, typename Pack>
void TLowModeProjMILCFull<FImpl,Pack>::setup(void)
{
    auto        &action     = envGet(FMat, par().action);
    int         Ls          = env().getObjectLs(par().action);
    bool project = par().projector;
    if (Ls > 1) {
       HADRONS_ERROR(Argument, "Ls > 1 not implemented");
    }

    envCache(FermionField, "rbFerm", 1, envGetRbGrid(FermionField));
    envCache(FermionField, "rbFermNeg", 1, envGetRbGrid(FermionField));
    envCache(FermionField, "MrbFermNeg", 1, envGetRbGrid(FermionField));
    envCache(FermionField, "rbTemp", 1, envGetRbGrid(FermionField));
    envCache(FermionField, "rbTempNeg", 1, envGetRbGrid(FermionField));

    auto &rbFerm     = envGet(FermionField, "rbFerm");
    auto &rbFermNeg  = envGet(FermionField, "rbFermNeg");
    auto &MrbFermNeg = envGet(FermionField, "MrbFermNeg");

    rbFerm     = Zero();
    rbFermNeg  = Zero();
    MrbFermNeg = Zero();

    LOG(Message) << "Setting up low mode projector "
                 << "for action '" << par().action 
                 << "' using eigenvectors from '" << par().lowModes 
                 << "'" << std::endl;

    auto &mat       = envGet(FMat, par().action);
    auto &epack   = envGet(Pack, par().lowModes);

    auto makeSolver    = [&mat, &epack, project, this] (bool subGuess) {
        return [&mat, &epack, subGuess, project, this] (FermionField &sol, const FermionField &source) {
            SOLVER_BODY;
        };
    };

    auto solver = makeSolver(false);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName(), Ls, solver, mat);
    envCreate(Solver, getName()+"_subtract", Ls, solver_subtract, mat);
}

#undef SOLVER_BODY

/******************************************************************************
 *              TLowModeProjMILCFull execution                                     *
 ******************************************************************************/

template <typename FImpl, typename Pack>
void TLowModeProjMILCFull<FImpl,Pack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_LowModeProj_hpp_
