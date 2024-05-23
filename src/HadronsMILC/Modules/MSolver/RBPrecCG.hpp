/*
 * RBPrecCGMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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

#ifndef HadronsMILC_MSolver_RBPrecCG_hpp_
#define HadronsMILC_MSolver_RBPrecCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Schur red-black preconditioned CG                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class RBPrecCGMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RBPrecCGMILCPar ,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    std::string , isEven,
                                    std::string , guesser);
};

template <typename FImpl>
class TRBPrecCGMILC: public Module<RBPrecCGMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    HADRONS_DEFINE_SCHUR_SOLVE(schurSolve_t,FImpl);
private:
    class GuessWrapper: public LinearFunction<FermionField> {
    public:
      GuessWrapper(const FermionField& guess)
      :guess_(guess){}
      using LinearFunction<FermionField>::operator();
      virtual void operator()(const FermionField &src, FermionField &guess) { guess = guess_; };
    private:
        const FermionField& guess_;
    };
public:
    // constructor
    TRBPrecCGMILC(const std::string name);
    // destructor
    virtual ~TRBPrecCGMILC(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool isEven_ {false};
};

MODULE_REGISTER_TMP(StagRBPrecCG, ARG(TRBPrecCGMILC<STAGIMPL>), MSolver);

/******************************************************************************
 *                      TRBPrecCGMILC template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRBPrecCGMILC<FImpl>::TRBPrecCGMILC(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRBPrecCGMILC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    if (!par().guesser.empty())
    {
        in.push_back(par().guesser);
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TRBPrecCGMILC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};
    
    return out;
}

template <typename FImpl>
DependencyMap TRBPrecCGMILC<FImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().action, getName()});
    dep.insert({par().action, getName() + "_subtract"});
    if (!par().guesser.empty())
    {
        dep.insert({par().guesser, getName(),             });
        dep.insert({par().guesser, getName() + "_subtract"});
    }

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
// C++11 does not support template lambdas so it is easier
// to make a macro with the solver body
#define SOLVER_BODY                                                                          \
int cb = this->isEven_?Even:Odd;                                                             \
ZeroGuesser<FermionField>    defaultGuesser;                                                 \
LinearFunction<FermionField> &guesser = (guesserPt == nullptr) ? defaultGuesser : *guesserPt;\
ConjugateGradient<FermionField> cg(par().residual,                                           \
                                   par().maxIteration);                                      \
schurSolve_t<FermionField> schurSolver(cg,false,false,cb);                                   \
schurSolver.subtractGuess(subGuess);                                                         \
schurSolver(mat, source, sol, guesser);

#define SOLVER_IMPROVE_BODY                                                                  \
int cbNeg = this->isEven_?Odd:Even;                                                          \
FermionField test(envGetGrid(FermionField));                                                 \
FermionField solRb(envGetRbGrid(FermionField)), solRbNeg(envGetRbGrid(FermionField));        \
LOG(Message) << "Improving residual of complementary checkerboard." << std::endl;            \
schurSolve_t<FermionField> schurSolverNeg(cg,false,false,cbNeg);                             \
schurSolverNeg.subtractGuess(subGuess);                                                      \
pickCheckerboard(cb,solRb,sol);                                                              \
pickCheckerboard(cbNeg,solRbNeg,sol);                                                        \
GuessWrapper guessNeg(solRbNeg);                                                             \
schurSolverNeg(mat, source, sol, guessNeg);                                                  \
setCheckerboard(sol,solRb);                                                                  \
mat.M(sol,test);                                                                             \
test = test-source;                                                                          \
RealD ns = norm2(source);                                                                    \
RealD nr = norm2(test);                                                                      \
LOG(Message) << "..........Combining checkerboards." << std::endl;                           \
LOG(Message) << "Final combined true residual: "<< std::sqrt(nr/ns) << std::endl;

template <typename FImpl>
void TRBPrecCGMILC<FImpl>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up Schur red-black preconditioned CG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << ", maximum iteration " 
                 << par().maxIteration << std::endl;

    auto                         Ls         = env().getObjectLs(par().action);
    auto                         &mat       = envGet(FMat, par().action);
    LinearFunction<FermionField> *guesserPt = nullptr;

    if (!par().isEven.empty()) {
        if (!(std::istringstream(par().isEven) >> std::boolalpha >> isEven_)) {
            HADRONS_ERROR(Logic,"parameter isEven='" + par().isEven + "' must be 'true' or 'false'");
        }
    }

    if (!par().guesser.empty())
    {
        guesserPt = &envGet(LinearFunction<FermionField>, par().guesser);
    }

    auto makeSolver = [&mat, guesserPt, this](bool subGuess) {
        return [&mat, guesserPt, subGuess, this]
        (FermionField &sol, const FermionField &source) 
        {
            SOLVER_BODY;

            SOLVER_IMPROVE_BODY;
        };
    };
    auto makeVecSolver = [&mat, guesserPt, this](bool subGuess) {
        return [&mat, guesserPt, subGuess, this]
        (std::vector<FermionField> &sol, const std::vector<FermionField> &source) 
        {
            SOLVER_BODY;

            LOG(Warning) << "Vector solver does not improve residual of complementary checkerboard. Desired residual not reached." << std::endl;

        };
    };
    auto makeGuessSolver = [&mat, this](bool subGuess) {
        return [&mat, subGuess, this](FermionField &sol,
                                     const FermionField &source, const FermionField &guess) {

            LinearFunction<FermionField> *guesserPt = nullptr;
            FermionField guessRb(envGetRbGrid(FermionField));
            {
                int cb = this->isEven_?Even:Odd;
                conformable(source,guess);
                pickCheckerboard(cb,guessRb,guess);
                guesserPt = new GuessWrapper(guessRb);
            }

            SOLVER_BODY;

            SOLVER_IMPROVE_BODY;
        };
    };    
    auto solver    = makeSolver(false);
    auto vecSolver = makeVecSolver(false);
    auto guessSolver = makeGuessSolver(false);
    envCreate(Solver, getName(), Ls, solver, vecSolver, guessSolver, mat);
    auto solver_subtract    = makeSolver(true);
    auto vecSolver_subtract = makeVecSolver(true);
    auto guessSolver_subtract = makeGuessSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, vecSolver_subtract, guessSolver_subtract, mat);
}

#undef SOLVER_BODY
#undef SOLVER_IMPROVE_BODY

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRBPrecCGMILC<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_RBPrecCG_hpp_