/*
 * MixedPrecisionRBPrecCGMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef HadronsMILC_MSolver_MixedPrecisionRBPrecCG_hpp_
#define HadronsMILC_MSolver_MixedPrecisionRBPrecCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Mixed precision schur red-black preconditioned CG             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MixedPrecisionRBPrecCGMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MixedPrecisionRBPrecCGMILCPar,
                                    std::string , innerAction,
                                    std::string , outerAction,
                                    unsigned int, maxInnerIteration,
                                    unsigned int, maxOuterIteration,
                                    double      , residual,
                                    std::string , isEven,
                                    std::string , innerGuesser,
                                    std::string , outerGuesser);
};

template <typename FImplInner, typename FImplOuter>
class TMixedPrecisionRBPrecCGMILC: public Module<MixedPrecisionRBPrecCGMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
    HADRONS_DEFINE_SCHUR_SOLVE(schurSolve_t,FImplOuter);
private:
    template <typename Field>
    class OperatorFunctionWrapper: public OperatorFunction<Field>
    {
    public:
        using OperatorFunction<Field>::operator();
        OperatorFunctionWrapper(LinearFunction<Field> &fn): fn_(fn) {};
        virtual ~OperatorFunctionWrapper(void) = default;
        virtual void operator()(LinearOperatorBase<Field> &op, 
                                const Field &in, Field &out)
        {
            fn_(in, out);
        }
    private:
        LinearFunction<Field> &fn_;
    };
    class GuessWrapper: public LinearFunction<FermionFieldOuter> {
    public:
      GuessWrapper(const FermionFieldOuter& guess)
      :guess_(guess){}
      using LinearFunction<FermionFieldOuter>::operator();
      virtual void operator()(const FermionFieldOuter &src, FermionFieldOuter &guess) { guess = guess_; };
    private:
        const FermionFieldOuter& guess_;
    };
public:
    // constructor
    TMixedPrecisionRBPrecCGMILC(const std::string name);
    // destructor
    virtual ~TMixedPrecisionRBPrecCGMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
public:
    bool isEven_ {false};
};

MODULE_REGISTER_TMP(MixedPrecisionRBPrecCGMILC, ARG(TMixedPrecisionRBPrecCGMILC<STAGIMPLF,STAGIMPL>), MSolver);

/******************************************************************************
 *                 TMixedPrecisionRBPrecCGMILC implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter>
TMixedPrecisionRBPrecCGMILC<FImplInner, FImplOuter>::TMixedPrecisionRBPrecCGMILC(const std::string name)
: Module<MixedPrecisionRBPrecCGMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter>
std::vector<std::string> TMixedPrecisionRBPrecCGMILC<FImplInner, FImplOuter>::getInput(void)
{
    std::vector<std::string> in = {par().innerAction, par().outerAction};
    
    if (!par().innerGuesser.empty())
    {
        in.push_back(par().innerGuesser);
    }
    if (!par().outerGuesser.empty())
    {
        in.push_back(par().outerGuesser);
    }

    return in;
}

template <typename FImplInner, typename FImplOuter>
std::vector<std::string> TMixedPrecisionRBPrecCGMILC<FImplInner, FImplOuter>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};

    return out;
}

template <typename FImplInner, typename FImplOuter>
DependencyMap TMixedPrecisionRBPrecCGMILC<FImplInner, FImplOuter>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().innerAction, getName()});
    dep.insert({par().innerAction, getName() + "_subtract"});
    dep.insert({par().outerAction, getName()});
    dep.insert({par().outerAction, getName() + "_subtract"});
    if (!par().innerGuesser.empty())
    {
        dep.insert({par().innerGuesser, getName(),             });
        dep.insert({par().innerGuesser, getName() + "_subtract"});
    }
    if (!par().outerGuesser.empty())
    {
        dep.insert({par().outerGuesser, getName(),             });
        dep.insert({par().outerGuesser, getName() + "_subtract"});
    }

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
// C++11 does not support template lambdas so it is easier
// to make a macro with the solver body
#define SOLVER_BODY                                                                                   \
int cb    = this->isEven_?Even:Odd;                                                                   \
ZeroGuesser<FermionFieldInner> iguesserDefault;                                                       \
ZeroGuesser<FermionFieldOuter> oguesserDefault;                                                       \
LinearFunction<FermionFieldInner> &iguesser = (iguesserPt == nullptr) ? iguesserDefault : *iguesserPt;\
LinearFunction<FermionFieldOuter> &oguesser = (oguesserPt == nullptr) ? oguesserDefault : *oguesserPt;\
FSchurOpInner simat(imat);                                                                            \
FSchurOpOuter somat(omat);                                                                            \
MixedPrecisionConjugateGradient<FermionFieldOuter, FermionFieldInner>                                 \
    mpcg(par().residual, par().maxInnerIteration,                                                     \
         par().maxOuterIteration,                                                                     \
         getGrid<FermionFieldInner>(true, Ls),                                                        \
         simat, somat);                                                                               \
mpcg.useGuesser(iguesser);                                                                            \
OperatorFunctionWrapper<FermionFieldOuter> wmpcg(mpcg);                                               \
schurSolve_t<FermionFieldOuter> schurSolver(wmpcg,false,false,cb);                                    \
schurSolver.subtractGuess(subGuess);                                                                  \
schurSolver(omat, source, sol, oguesser);

#define SOLVER_IMPROVE_BODY                                                                           \
int cbNeg = this->isEven_?Odd:Even;                                                                   \
FermionFieldOuter test(envGetGrid(FermionFieldOuter));                                                \
FermionFieldOuter solRb(envGetRbGrid(FermionFieldOuter)), solRbNeg(envGetRbGrid(FermionFieldOuter));  \
LOG(Message) << "Improving residual of complementary checkerboard." << std::endl;                     \
schurSolve_t<FermionFieldOuter> schurSolverNeg(wmpcg,false,false,cbNeg);                              \
schurSolverNeg.subtractGuess(subGuess);                                                               \
pickCheckerboard(cb,solRb,sol);                                                                       \
pickCheckerboard(cbNeg,solRbNeg,sol);                                                                 \
GuessWrapper guessNeg(solRbNeg);                                                                      \
schurSolverNeg(omat, source, sol, guessNeg);                                                          \
setCheckerboard(sol,solRb);                                                                           \
omat.M(sol,test);                                                                                     \
test = test-source;                                                                                   \
RealD ns = norm2(source);                                                                             \
RealD nr = norm2(test);                                                                               \
LOG(Message) << "..........Combining checkerboards." << std::endl;                                    \
LOG(Message) << "Final combined true residual: "<< std::sqrt(nr/ns) << std::endl;

template <typename FImplInner, typename FImplOuter>
void TMixedPrecisionRBPrecCGMILC<FImplInner, FImplOuter>::setup(void)
{
    LOG(Message) << "Setting up Schur red-black preconditioned mixed-precision "
                 << "CG for inner/outer action '" << par().innerAction 
                 << "'/'" << par().outerAction << "', residual "
                 << par().residual << ", and maximum inner/outer iteration " 
                 << par().maxInnerIteration << "/" << par().maxOuterIteration
                 << std::endl;

    auto                              Ls          = env().getObjectLs(par().innerAction);
    auto                              &imat       = envGet(FMatInner, par().innerAction);
    auto                              &omat       = envGet(FMatOuter, par().outerAction);
    LinearFunction<FermionFieldInner> *iguesserPt = nullptr; 
    LinearFunction<FermionFieldOuter> *oguesserPt = nullptr;

    if (!par().isEven.empty()) {
        if (!(std::istringstream(par().isEven) >> std::boolalpha >> isEven_)) {
            HADRONS_ERROR(Logic,"parameter isEven='" + par().isEven + "' must be 'true' or 'false'");
        }
    }

    if (!par().innerGuesser.empty())
    {
        iguesserPt = &envGet(LinearFunction<FermionFieldInner>, par().innerGuesser);
    }
    if (!par().outerGuesser.empty())
    {
        oguesserPt = &envGet(LinearFunction<FermionFieldOuter>, par().outerGuesser);
    }
    auto makeSolver = [&imat, &omat, iguesserPt, oguesserPt, Ls, this](bool subGuess)
    {
        return [&imat, &omat, iguesserPt, oguesserPt, subGuess, Ls, this]
            (FermionFieldOuter &sol, const FermionFieldOuter &source) 
        {
            SOLVER_BODY;

            SOLVER_IMPROVE_BODY;
        };
    };
    auto makeVecSolver = [&imat, &omat, iguesserPt, oguesserPt, Ls, this](bool subGuess)
    {
        return [&imat, &omat, iguesserPt, oguesserPt, subGuess, Ls, this]
            (std::vector<FermionFieldOuter> &sol, const std::vector<FermionFieldOuter> &source) 
        {
            SOLVER_BODY;

            LOG(Warning) << "Vector solver does not improve residual of complementary checkerboard. Desired residual not reached." << std::endl;
        };
    };

    auto makeGuessSolver = [&imat, &omat, iguesserPt, Ls, this](bool subGuess) {
        return [&imat, &omat, iguesserPt, subGuess, Ls, this](FermionFieldOuter &sol,
                                     const FermionFieldOuter &source, const FermionFieldOuter &guess) {

            LinearFunction<FermionFieldOuter> *oguesserPt = nullptr;

            FermionFieldOuter guessRb(envGetRbGrid(FermionFieldOuter));
            {
                int cb = this->isEven_?Even:Odd;
                conformable(source,guess);
                pickCheckerboard(cb,guessRb,guess);
                oguesserPt = new GuessWrapper(guessRb);
            }

            SOLVER_BODY;

            SOLVER_IMPROVE_BODY;
        };
    };

    auto solver    = makeSolver(false);
    auto vecSolver = makeVecSolver(false);
    auto guessSolver = makeGuessSolver(false);
    envCreate(Solver, getName(), Ls, solver, vecSolver, guessSolver, omat);
    auto solver_subtract    = makeSolver(true);
    auto vecSolver_subtract = makeVecSolver(true);
    auto guessSolver_subtract = makeGuessSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, vecSolver_subtract, guessSolver_subtract, omat);
}

#undef SOLVER_BODY
#undef SOLVER_IMPROVE_BODY

// execution ///////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter>
void TMixedPrecisionRBPrecCGMILC<FImplInner, FImplOuter>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MixedPrecisionRBPrecCG_hpp_
