/*
 * PropToFermions.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
 */

/*  END LEGAL */
#ifndef HadronsMILC_MUtilities_PropToFermions_hpp_
#define HadronsMILC_MUtilities_PropToFermions_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *          PropagatorField -> per-color FermionFields adapter               *
 ******************************************************************************/
/*  Extracts the color columns of a scalar PropagatorField into a
    std::vector<FermionField> (one entry per color index, FImpl::Dimension
    entries) via PropToFerm -- the inverse of FermToProp assembly. This
    adapts PropagatorField producers (e.g. MFermion::StagLMAMesonFieldProp,
    whose per-(t,gamma) outputs are scalar propagators) to consumers that
    take a noise-vector-like object, such as the MContraction::
    StagA2AMesonField probe contraction (vector<FermionField> only):
    entry c of the output equals the fermion reconstructed from the
    color-diluted source's column noiseIndex + c, so probe contraction
    column c is directly comparable to the reference chain's noise column
    noiseIndex + c.

    prop        name of the PropagatorField environment object to adapt
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class PropToFermionsPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(PropToFermionsPar, std::string, prop);
};

template <typename FImpl>
class TPropToFermions : public Module<PropToFermionsPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TPropToFermions(const std::string name);
  // destructor
  virtual ~TPropToFermions(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(PropToFermions, TPropToFermions<STAGIMPL>, MUtilities);

/******************************************************************************
 *                       TPropToFermions implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPropToFermions<FImpl>::TPropToFermions(const std::string name)
    : Module<PropToFermionsPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPropToFermions<FImpl>::getInput(void) {
  std::vector<std::string> in = {par().prop};

  return in;
}

template <typename FImpl>
std::vector<std::string> TPropToFermions<FImpl>::getOutput(void) {
  std::vector<std::string> out = {getName()};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropToFermions<FImpl>::setup(void) {
  // one FermionField per color index, sized from the grid only -- dry-run
  // safe (no input-content reads). The count is bound to a local first:
  // passing FImpl::Dimension straight through envCreate's forwarding
  // references would odr-use the constexpr static member (no out-of-line
  // definition in StaggeredImpl) -- EigenPackFullPairs local-variable
  // pattern
  const unsigned int nColor = FImpl::Dimension;
  envCreate(std::vector<FermionField>, getName(), 1, nColor,
            envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropToFermions<FImpl>::execute(void) {
  auto &prop = envGet(PropagatorField, par().prop);
  auto &fermions = envGet(std::vector<FermionField>, getName());

  LOG(Message) << "Extracting " << FImpl::Dimension
               << " color FermionFields from propagator '" << par().prop
               << "'" << std::endl;
  for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
    PropToFerm<FImpl>(fermions[c], prop, c);
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MUtilities_PropToFermions_hpp_
