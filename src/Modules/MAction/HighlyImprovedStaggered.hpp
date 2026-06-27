/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MAction/HighlyImprovedStaggered.hpp

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef HadronsMILC_MAction_HighlyImprovedStaggered_hpp_
#define HadronsMILC_MAction_HighlyImprovedStaggered_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/HighlyImprovedStaggeredFermionImpl.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       HighlyImprovedStaggeredMILC                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class HighlyImprovedStaggeredMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HighlyImprovedStaggeredMILCPar,
                                    std::string, gauge,
                                    double     , mass,
                                    std::string, boundary);
};

template <typename FImpl>
class THighlyImprovedStaggeredMILC: public Module<HighlyImprovedStaggeredMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    THighlyImprovedStaggeredMILC(const std::string name);
    // destructor
    virtual ~THighlyImprovedStaggeredMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(HighlyImprovedStaggeredMILC, THighlyImprovedStaggeredMILC<STAGIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(HighlyImprovedStaggeredMILCF, THighlyImprovedStaggeredMILC<STAGIMPLF>, MAction);
#endif

/******************************************************************************
 *              THighlyImprovedStaggeredMILC implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
THighlyImprovedStaggeredMILC<FImpl>::THighlyImprovedStaggeredMILC(const std::string name)
: Module<HighlyImprovedStaggeredMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> THighlyImprovedStaggeredMILC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> THighlyImprovedStaggeredMILC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void THighlyImprovedStaggeredMILC<FImpl>::setup(void)
{
    LOG(Message) << "Setting up HighlyImprovedStaggered fermion matrix." << std::endl;
    LOG(Message) << "Using m=" << par().mass << std::endl;
    LOG(Message) << "Using thin links: " << par().gauge << std::endl;
    
    auto &U      = envGet(GaugeField, par().gauge);
    auto &grid   = *envGetGrid(FermionField);
    auto &gridRb = *envGetRbGrid(FermionField);

    // Staggered impl params: boundary phases drive calcStagPhases/rephase.
    // Default is anti-periodic in time (the HISQ impl's own 1-arg default),
    // matching the usual MILC physics convention.
    StaggeredImplParams stagParams;
    if (!par().boundary.empty())
    {
        stagParams.boundary_phases = strToVec<Complex>(par().boundary);
    }
    else
    {
        stagParams.boundary_phases = std::vector<Complex>{Complex(1.), Complex(1.), Complex(1.), Complex(-1.)};
    }

    // HISQ smearing implementation. The constructor computes the KS phases
    // (and folds in the boundary phases) into stagPhases, used by rephase().
    HighlyImprovedStaggeredFermionImpl<FImpl> hisq(&grid, stagParams);

    // Two-level HISQ smear (cf. Grid/tests/smearing/Test_fatLinks.cc):
    //   rephase  -> bake KS phases + boundary into the thin gauge
    //   fat7     -> level-1 fat links (no Naik)
    //   project  -> U(3) unitary projection of the fat links
    //   asqtad   -> level-2 smear: Ufat = fat links, Ulong = Naik (long) links
    // Unlike TImprovedStaggeredMILC, which assumed its inputs already carried
    // the KS phases and antiperiodic boundary, this module applies them here.
    GaugeField R(&grid), V(&grid), W(&grid), Ufat(&grid), Ulong(&grid);
    hisq.rephase(R, U);
    hisq.smear(V, R);
    hisq.project(W, V);
    hisq.smear(Ufat, Ulong, W);

    // ImportGaugeSimple stores the fat & Naik links raw (the DhopImproved
    // kernel weights both 1- and 3-link hops at 1.0, with no c1/c2/u0 scaling).
    // Bake the symmetric-hop 0.5 normalisation into the links, per direction,
    // mirroring ImportGauge's own U*(0.5*c1/u0) idiom. With c1=c2=u0=1 this
    // reproduces the standard ImprovedStaggeredFermion operator.
    typename FImpl::GaugeLinkField link(&grid);
    for (int mu = 0; mu < Nd; ++mu)
    {
        link = PeekIndex<LorentzIndex>(Ufat, mu);
        PokeIndex<LorentzIndex>(Ufat, link*RealD(0.5), mu);

        link = PeekIndex<LorentzIndex>(Ulong, mu);
        PokeIndex<LorentzIndex>(Ulong, link*RealD(0.5), mu);
    }

    // Build the fermion matrix. boundary_phases on the matrix ImplParams is
    // inert under ImportGaugeSimple (only DoubleStore/ImportGauge apply it),
    // so the default constructor params are used. mass is passed MILC-style
    // (2x) to match the existing TImprovedStaggeredMILC convention.
    envCreateDerived(FMat, ImprovedStaggeredFermion<FImpl>, getName(), 1,
                     grid, gridRb, 2.*par().mass);

    auto &fmat = envGetDerived(FMat, ImprovedStaggeredFermion<FImpl>, getName());
    fmat.ImportGaugeSimple(Ulong, Ufat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void THighlyImprovedStaggeredMILC<FImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MAction_HighlyImprovedStaggered_hpp_
