/*
 * CGMILC.hpp, part of HadronsMILC
 *
 * Standard (single-precision) conjugate-gradient solver module for MILC
 * staggered actions.
 */

#ifndef HadronsMILC_MSolver_CG_hpp_
#define HadronsMILC_MSolver_CG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Standard CG solver (CGNE)                          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class CGMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CGMILCPar, std::string, action, unsigned int,
                                  maxIteration, double, residual, std::string,
                                  guesser);
};

template <typename FImpl> class TCGMILC : public Module<CGMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );
  SOLVER_TYPE_ALIASES(FImpl, );

private:
  class GuessWrapper : public LinearFunction<FermionField> {
  public:
    GuessWrapper(const FermionField &guess) : guess_(guess) {}
    using LinearFunction<FermionField>::operator();
    virtual void operator()(const FermionField &src, FermionField &guess) {
      guess = guess_;
    };

  private:
    const FermionField &guess_;
  };

public:
  // constructor
  TCGMILC(const std::string name);
  // destructor
  virtual ~TCGMILC(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  virtual DependencyMap getObjectDependencies(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(StagCGMILC, ARG(TCGMILC<STAGIMPL>), MSolver);

/******************************************************************************
 *                           TCGMILC implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TCGMILC<FImpl>::TCGMILC(const std::string name) : Module<CGMILCPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TCGMILC<FImpl>::getInput(void) {
  std::vector<std::string> in = {par().action};

  if (!par().guesser.empty()) {
    in.push_back(par().guesser);
  }

  return in;
}

template <typename FImpl>
std::vector<std::string> TCGMILC<FImpl>::getOutput(void) {
  std::vector<std::string> out = {getName(), getName() + "_subtract"};

  return out;
}

template <typename FImpl>
DependencyMap TCGMILC<FImpl>::getObjectDependencies(void) {
  DependencyMap dep;

  dep.insert({par().action, getName()});
  dep.insert({par().action, getName() + "_subtract"});
  if (!par().guesser.empty()) {
    dep.insert({par().guesser, getName()});
    dep.insert({par().guesser, getName() + "_subtract"});
  }

  return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
// C++11 does not support template lambdas so it is easier
// to make a macro with the solver body
#define SOLVER_BODY                                                            \
  GridBase *g = sol.Grid();                                                    \
  FermionField tmp(g);                                                         \
  MdagMLinearOperator<FMat, FermionField> hermOp(mat);                         \
  ConjugateGradient<FermionField> cg(par().residual, par().maxIteration);      \
  ZeroGuesser<FermionField> defaultGuesser;                                    \
  LinearFunction<FermionField> &guesser =                                      \
      (guesserPt == nullptr) ? defaultGuesser : *guesserPt;                    \
  mat.Mdag(source, tmp);                                                       \
  guesser(source, sol);                                                        \
  cg(hermOp, tmp, sol);                                                        \
  mat.M(sol, tmp);                                                             \
  tmp = tmp - source;                                                          \
  RealD ns = norm2(source);                                                    \
  RealD nr = norm2(tmp);                                                       \
  LOG(Message) << "Final true residual: " << std::sqrt(nr / ns) << std::endl;  \
  if (subGuess) {                                                              \
    guesser(source, sol);                                                      \
    sol = sol - tmp;                                                           \
  }

template <typename FImpl> void TCGMILC<FImpl>::setup(void) {
  if (par().maxIteration == 0) {
    HADRONS_ERROR(Argument, "zero maximum iteration");
  }

  LOG(Message) << "setting up standard CG for action '" << par().action
               << "' with residual " << par().residual << ", maximum iteration "
               << par().maxIteration << std::endl;

  auto Ls = env().getObjectLs(par().action);
  auto &mat = envGet(FMat, par().action);
  LinearFunction<FermionField> *guesserPt = nullptr;

  if (!par().guesser.empty()) {
    guesserPt = &envGet(LinearFunction<FermionField>, par().guesser);
  }
  auto makeSolver = [&mat, guesserPt, this](bool subGuess) {
    return [&mat, guesserPt, subGuess, this](
               FermionField &sol, const FermionField &source) { SOLVER_BODY; };
  };

  auto makeGuessSolver = [&mat, this](bool subGuess) {
    return [&mat, subGuess, this](FermionField &sol, const FermionField &source,
                                  const FermionField &guess) {
      LinearFunction<FermionField> *guesserPt = nullptr;
      GuessWrapper guessWrapper(guess);
      guesserPt = &guessWrapper;

      SOLVER_BODY;
    };
  };

  auto solver = makeSolver(false);
  auto guessSolver = makeGuessSolver(false);
  envCreate(Solver, getName(), Ls, solver, guessSolver, mat);
  auto solver_subtract = makeSolver(true);
  auto guessSolver_subtract = makeGuessSolver(true);
  envCreate(Solver, getName() + "_subtract", Ls, solver_subtract,
            guessSolver_subtract, mat);
}

#undef SOLVER_BODY

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl> void TCGMILC<FImpl>::execute(void) {}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MSolver_CG_hpp_
