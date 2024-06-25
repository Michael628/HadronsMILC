/*
 * LowModeProjMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

class LowModeProjMILCPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LowModeProjMILCPar,
                                  std::string, action,
                                  bool, projector,
                                  unsigned int, eigStart,
                                  int, nEigs,
                                  std::string, lowModes);
};

template <typename FImpl, typename Pack>
class TLowModeProjMILC : public Module<LowModeProjMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLowModeProjMILC(const std::string name);
    // destructor
    virtual ~TLowModeProjMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);

    // setup
    virtual void setup(void);

    // execute
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StagLMA, ARG(TLowModeProjMILC<STAGIMPL,MassShiftEigenPack<STAGIMPL> >), MSolver);

/******************************************************************************
 *                       TLowModeProjMILC implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TLowModeProjMILC<FImpl,Pack>::TLowModeProjMILC(const std::string name)
: Module<LowModeProjMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TLowModeProjMILC<FImpl,Pack>::getInput(void)
{
    std::vector<std::string> in {par().action,par().lowModes};
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLowModeProjMILC<FImpl,Pack>::getOutput(void)
{
    std::vector<std::string> out {getName(),getName()+"_subtract"};
    return out;
}

template <typename FImpl, typename Pack>
DependencyMap TLowModeProjMILC<FImpl,Pack>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().action, getName()});
    dep.insert({par().lowModes, getName()});
    dep.insert({par().action, getName()+"_subtract"});
    dep.insert({par().lowModes, getName()+"_subtract"});

    return dep;
}

/******************************************************************************
 *              TLowModeProjMILC setup                                         *
 ******************************************************************************/


template <typename FImpl, typename Pack>
void TLowModeProjMILC<FImpl,Pack>::setup(void)
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

    auto eigStart = par().eigStart;
    auto nEigs = par().nEigs;

    if (nEigs < 1) {
        nEigs = epack.evec.size();
    }
    
    if (eigStart > nEigs || eigStart > epack.evec.size() || nEigs - eigStart > epack.evec.size() - eigStart) {
        HADRONS_ERROR(Argument,"Requested eigs (parameters eigStart and nEigs) out of bounds.")
    }

    auto makeSolver    = [&mat, &epack, project, eigStart, nEigs, this] (bool subGuess) {
        return [&mat, &epack, subGuess, project, eigStart, nEigs, this] (FermionField &sol, const FermionField &source) {
            auto &rbTemp = envGet(FermionField,"rbTemp");
            auto &rbTempNeg = envGet(FermionField,"rbTempNeg");
            auto &rbFerm = envGet(FermionField,"rbFerm");
            auto &rbFermNeg = envGet(FermionField,"rbFermNeg");
            auto &MrbFermNeg = envGet(FermionField,"MrbFermNeg");

            int cb = epack.evec[0].Checkerboard();
            int cbNeg = (cb==Even) ? Odd : Even;

            RealD norm = 1./::sqrt(norm2(epack.evec[0]));

            rbTemp = Zero();
            rbTemp.Checkerboard() = cb;
            rbTempNeg = Zero();
            rbTempNeg.Checkerboard() = cb;

            rbFerm.Checkerboard() = cb;
            rbFermNeg.Checkerboard() = cbNeg;
            MrbFermNeg.Checkerboard() = cb;
            {
              pickCheckerboard(cb,rbFerm,source);
              pickCheckerboard(cbNeg,rbFermNeg,source);
            }
            mat.MeooeDag(rbFermNeg, MrbFermNeg);

            for (int k=(eigStart+nEigs-1) ; k >= int(eigStart) ; k--) {
                const FermionField& e = epack.evec[k];

                const RealD mass     = epack.eval[k].real();
                const RealD lam_D    = epack.eval[k].imag();
                const RealD invlam_D = 1./lam_D;
                const RealD invmag   = 1./(pow(mass,2)+pow(lam_D,2));

                if (!project) {
                    const ComplexD ip    = TensorRemove(innerProduct(e,rbFerm))*invmag;
                    const ComplexD ipNeg = TensorRemove(innerProduct(e,MrbFermNeg))*invmag;
                    axpy(rbTemp,    mass*ip+ipNeg,   e,rbTemp);
                    axpy(rbTempNeg, mass*ipNeg*invlam_D*invlam_D-ip, e,rbTempNeg);
                } else {
                    const ComplexD ip    = TensorRemove(innerProduct(e,rbFerm));
                    const ComplexD ipNeg = TensorRemove(innerProduct(e,MrbFermNeg));
                    axpy(rbTemp,    ip,   e,rbTemp);
                    axpy(rbTempNeg, ipNeg*invlam_D*invlam_D, e,rbTempNeg);
                }
            }

            mat.Meooe(rbTempNeg, rbFermNeg);
            {
              setCheckerboard(sol,rbTemp);
              setCheckerboard(sol,rbFermNeg);
            }
            sol *= norm;
            if (subGuess) {
                if (project) {
                    sol = source - sol;
                } else {
                   HADRONS_ERROR(Argument, "Subtracted solver only supported for projector=true");
               }
            }
        };
    };

    auto solver = makeSolver(false);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName(), Ls, solver, mat);
    envCreate(Solver, getName()+"_subtract", Ls, solver_subtract, mat);
}

#undef SOLVER_BODY

/******************************************************************************
 *              TLowModeProjMILC execution                                     *
 ******************************************************************************/

template <typename FImpl, typename Pack>
void TLowModeProjMILC<FImpl,Pack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_LowModeProj_hpp_
