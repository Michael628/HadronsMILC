/*
 * EigenPackCBPairs.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
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
 * (same GNU GPL as the repo headers; see the LICENSE file).
 */

/*  END LEGAL */
#ifndef HadronsMILC_MUtilities_EigenPackCBPairs_hpp_
#define HadronsMILC_MUtilities_EigenPackCBPairs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <EigenPack.hpp>
#include <A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *   On-demand checkerboarded evec pair source (function-producing module)   *
 ******************************************************************************/
/*  Produces the callable consumed by MContraction::StagA2AMesonField's
    'cbPairs' option. For a range [start, start+n) of a checkerboarded
    eigenpack's eigenvectors (eigenvalue 2m + i*lam, so Im(eval) = lam), the
    source fills module-managed grow-only buffers with

        even[k] = the even-CB part of evec[start+k]
        odd[k]  = the odd-CB part of evec[start+k]
                      = Meooe(even part)/(i*lam)

    -- the same e_!cb relation as MesonFieldLegacy's swapEvecCheckerFn and
    EigenPackFullPairs' pair construction, but ON DEMAND and non-destructive:
    the eigenpack's evec array is never touched, so the legacy SwapFn
    cache/swap/restore machinery is bypassed entirely and a re-executed
    eigensolver is picked up automatically on the next get() (no staleness).
    The product is created envCreateDerived under its framework interface
    (the ImprovedStaggered/FMat pattern), so consumers only need
    A2AMatrix.hpp's A2ALowModePairSourceMILC<FermionField> to envGet it --
    no concrete-type coupling.
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EigenPackCBPairsPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(EigenPackCBPairsPar,
                                  std::string, eigenPack,
                                  std::string, action);
};

template <typename FImpl, typename Pack>
class EigenPackCBPairSource
    : public A2ALowModePairSourceMILC<typename FImpl::FermionField> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  EigenPackCBPairSource(Pack &pack, FMat &action)
      : pack_(pack), action_(action) {}

  std::pair<const FermionField *, const FermionField *>
  get(const int start, const int n) override {
    if (n <= 0 || start < 0 || start + n > (int)pack_.evec.size()) {
      HADRONS_ERROR(Size, "EigenPackCBPairs: evec range [" +
                              std::to_string(start) + ", " +
                              std::to_string(start + n) +
                              ") invalid for pack of " +
                              std::to_string(pack_.evec.size()) + " evecs");
    }
    if ((int)even_.size() < n) {
      even_.resize(n, pack_.evec[0].Grid());
      odd_.resize(n, pack_.evec[0].Grid());
    }
    if (!temp_) {
      temp_ = std::make_unique<FermionField>(pack_.evec[0].Grid());
    }
    for (int k = 0; k < n; ++k) {
      const FermionField &e = pack_.evec[start + k];
      int cb = e.Checkerboard();
      int cbNeg = (cb == Even) ? Odd : Even;
      ComplexD eval_D = ComplexD(0., pack_.eval[start + k].imag());

      // e_!cb = Meooe(e_cb)/(i*lam). The Meooe output must be labeled cbNeg
      // BEFORE the call, and re-labeled after the expression-template
      // multiply (ET assignments do not propagate the checkerboard
      // attribute -- EigenPackFullPairs' documented lesson).
      *temp_ = Zero();
      temp_->Checkerboard() = cbNeg;
      action_.Meooe(e, *temp_);
      *temp_ = (1.0 / eval_D) * *temp_;
      temp_->Checkerboard() = cbNeg;

      if (cb == Even) {
        even_[k] = e;
        odd_[k] = *temp_;
      } else {
        odd_[k] = e;
        even_[k] = *temp_;
      }
      even_[k].Checkerboard() = Even;
      odd_[k].Checkerboard() = Odd;
    }
    return {even_.data(), odd_.data()};
  }

private:
  Pack &pack_;
  FMat &action_;
  std::vector<FermionField> even_, odd_; // grow-only, reused across get() calls
  std::unique_ptr<FermionField> temp_;   // lazily created Meooe scratch
};

template <typename FImpl, typename Pack>
class TEigenPackCBPairs : public Module<EigenPackCBPairsPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TEigenPackCBPairs(const std::string name);
  // destructor
  virtual ~TEigenPackCBPairs(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  virtual DependencyMap getObjectDependencies(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(EigenPackCBPairs,
                    ARG(TEigenPackCBPairs<STAGIMPL,
                                          MassShiftEigenPack<STAGIMPL>>),
                    MUtilities);

/******************************************************************************
 *                  TEigenPackCBPairs implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TEigenPackCBPairs<FImpl, Pack>::TEigenPackCBPairs(const std::string name)
    : Module<EigenPackCBPairsPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TEigenPackCBPairs<FImpl, Pack>::getInput(void) {
  std::vector<std::string> in{par().eigenPack, par().action};

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TEigenPackCBPairs<FImpl, Pack>::getOutput(void) {
  std::vector<std::string> out{getName()};

  return out;
}

template <typename FImpl, typename Pack>
DependencyMap TEigenPackCBPairs<FImpl, Pack>::getObjectDependencies(void) {
  DependencyMap dep;

  dep.insert({par().eigenPack, getName()});
  dep.insert({par().action, getName()});

  return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TEigenPackCBPairs<FImpl, Pack>::setup(void) {
  int Ls = env().getObjectLs(par().eigenPack);

  // Create under the framework interface type: consumers envGet
  // A2ALowModePairSourceMILC<FermionField> without knowing this header
  // (the ImprovedStaggered/FMat envCreateDerived pattern). ARG shields the
  // template comma from the macro preprocessor.
  envCreateDerived(A2ALowModePairSourceMILC<FermionField>,
                   ARG(EigenPackCBPairSource<FImpl, Pack>), getName(), Ls,
                   envGet(Pack, par().eigenPack), envGet(FMat, par().action));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TEigenPackCBPairs<FImpl, Pack>::execute(void) {
  auto &epack = envGet(Pack, par().eigenPack);

  LOG(Message) << "CB pair source '" << getName() << "' over "
               << epack.evec.size() << " evecs of '" << par().eigenPack
               << "' (on-demand Meooe fill)" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MUtilities_EigenPackCBPairs_hpp_
