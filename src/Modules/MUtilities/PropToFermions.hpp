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
/*  Extracts the color columns of PropagatorField(s) into a flat
    std::vector<FermionField> via PropToFerm -- the inverse of FermToProp
    assembly. Scalar input (a bare PropagatorField): FImpl::Dimension
    entries, one per color index. Vector input (a std::vector<
    PropagatorField>, e.g. MFermion::StagLMAMesonFieldProp's nNoise > 1
    per-(t,gamma) outputs): nNoise*Dimension entries, NOISE-MAJOR (entry
    n*Dimension + c) so the flat vector matches the meson-field table's
    noise-window column layout and the reference chain's per-noise-field
    columns one-to-one. This adapts propagator producers to consumers
    that take a noise-vector-like object, such as the MContraction::
    StagA2AMesonField probe contraction (vector<FermionField> only):
    with a color-diluted source, entry n*Dimension + c of the output
    equals the fermion reconstructed from table column
    noiseIndex + n*Dimension + c.

    prop        name of the PropagatorField or std::vector<PropagatorField>
                environment object to adapt
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
  // container dispatch (the MContraction::Meson setup pattern): a scalar
  // PropagatorField adapts to FImpl::Dimension fermions (the original
  // behavior); a std::vector<PropagatorField> (nNoise > 1 producers,
  // e.g. StagLMAMesonFieldProp) adapts to nNoise*Dimension fermions,
  // noise-major. The vector path cannot know nNoise at setup time
  // without reading the input object's size, so the output vector is
  // created empty and resized at execute (the RandomWall/MNoise
  // resize-at-execute precedent -- the memory profiler under-counts by
  // design there). The count is bound to a local first: passing
  // FImpl::Dimension straight through envCreate's forwarding references
  // would odr-use the constexpr static member (no out-of-line definition
  // in StaggeredImpl) -- EigenPackFullPairs local-variable pattern
  if (envHasType(PropagatorField, par().prop)) {
    const unsigned int nColor = FImpl::Dimension;
    envCreate(std::vector<FermionField>, getName(), 1, nColor,
              envGetGrid(FermionField));
  } else if (envHasType(std::vector<PropagatorField>, par().prop)) {
    envCreate(std::vector<FermionField>, getName(), 1, 0,
              envGetGrid(FermionField));
  } else {
    HADRONS_ERROR(Argument, "prop parameter '" + par().prop +
                                "' must be a PropagatorField or a "
                                "std::vector<PropagatorField>");
  }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropToFermions<FImpl>::execute(void) {
  auto &fermions = envGet(std::vector<FermionField>, getName());

  if (envHasType(PropagatorField, par().prop)) {
    auto &prop = envGet(PropagatorField, par().prop);

    LOG(Message) << "Extracting " << FImpl::Dimension
                 << " color FermionFields from propagator '" << par().prop
                 << "'" << std::endl;
    for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
      PropToFerm<FImpl>(fermions[c], prop, c);
    }
  } else if (envHasType(std::vector<PropagatorField>, par().prop)) {
    auto &props = envGet(std::vector<PropagatorField>, par().prop);
    const unsigned int nNoise = props.size();
    const unsigned int nColor = FImpl::Dimension;

    fermions.resize(nNoise * nColor, envGetGrid(FermionField));
    LOG(Message) << "Extracting " << nNoise << " x " << nColor
                 << " color FermionFields (noise-major) from " << nNoise
                 << " propagator(s) of '" << par().prop << "'" << std::endl;
    for (unsigned int n = 0; n < nNoise; ++n) {
      for (unsigned int c = 0; c < nColor; ++c) {
        PropToFerm<FImpl>(fermions[n * nColor + c], props[n], c);
      }
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MUtilities_PropToFermions_hpp_
