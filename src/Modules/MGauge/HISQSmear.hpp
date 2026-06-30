/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MGauge/HISQSmear.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Michael Lynch <michaellynch628@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#ifndef HadronsMILC_MGauge_HISQSmear_hpp_
#define HadronsMILC_MGauge_HISQSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Grid/qcd/utils/HighlyImprovedStaggeredFermionImpl.h>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              HISQSmear                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class HISQSmearPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(HISQSmearPar, std::string, gauge, std::string,
                                  boundary);
};

template <typename FImpl> class THISQSmear : public Module<HISQSmearPar> {
public:
  typedef typename FImpl::GaugeField GaugeField;
  typedef typename FImpl::GaugeLinkField GaugeLinkField;

public:
  // constructor
  THISQSmear(const std::string name);
  // destructor
  virtual ~THISQSmear(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);

protected:
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(HISQSmear, THISQSmear<STAGIMPL>, MGauge);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(HISQSmearF, THISQSmear<STAGIMPLF>, MGauge);
#endif

/******************************************************************************
 *                  THISQSmear implementation                                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
THISQSmear<FImpl>::THISQSmear(const std::string name)
    : Module<HISQSmearPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> THISQSmear<FImpl>::getInput(void) {
  std::vector<std::string> in = {par().gauge};

  return in;
}

template <typename FImpl>
std::vector<std::string> THISQSmear<FImpl>::getOutput(void) {
  std::vector<std::string> out = {getName() + "_fat", getName() + "_long"};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl> void THISQSmear<FImpl>::setup(void) {
  envCreateLat(GaugeField, getName() + "_fat");
  envCreateLat(GaugeField, getName() + "_long");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl> void THISQSmear<FImpl>::execute(void) {
  LOG(Message) << "Smearing gauge field '" << par().gauge
               << "' (HISQ: fat7 + project + asqtad)" << std::endl;

  auto &U = envGet(GaugeField, par().gauge);
  auto *grid = envGetGrid(GaugeField);

  // Staggered impl params: boundary phases drive calcStagPhases/rephase.
  // Default is anti-periodic in time (the HISQ impl's own 1-arg default),
  // matching the usual MILC physics convention.
  StaggeredImplParams stagParams;
  if (!par().boundary.empty()) {
    stagParams.boundary_phases = strToVec<Complex>(par().boundary);
  } else {
    stagParams.boundary_phases = std::vector<Complex>{
        Complex(1.), Complex(1.), Complex(1.), Complex(-1.)};
  }

  // HISQ smearing implementation. The constructor computes the KS phases
  // (and folds in the boundary phases) into stagPhases, used by rephase().
  HighlyImprovedStaggeredFermionImpl<FImpl> hisq(grid, stagParams);

  // Two-level HISQ smear (cf. Grid/tests/smearing/Test_fatLinks.cc):
  //   rephase  -> bake KS phases + boundary into the thin gauge
  //   fat7     -> level-1 fat links (no Naik)
  //   project  -> U(3) unitary projection of the fat links
  //   asqtad   -> level-2 smear: Ufat = fat links, Ulong = Naik (long) links
  auto &Ufat = envGet(GaugeField, getName() + "_fat");
  auto &Ulong = envGet(GaugeField, getName() + "_long");
  GaugeField R(grid), V(grid), W(grid);
  hisq.rephase(R, U);
  hisq.smear(V, R);
  hisq.project(W, V);
  hisq.smear(Ufat, Ulong, W);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MGauge_HISQSmear_hpp_
