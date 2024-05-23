/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MAction/ImprovedStaggered.hpp

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

#ifndef HadronsMILC_MAction_ImprovedStaggeredMILC5D_hpp_
#define HadronsMILC_MAction_ImprovedStaggeredMILC5D_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ImprovedStaggeredMILC5D                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class ImprovedStaggeredMILC5DPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ImprovedStaggeredMILC5DPar,
                                    std::string, gauge,
                                    std::string, gaugefat,
                                    std::string, gaugelong,
                                    double     , mass,
                                    double     , c1,
                                    double     , c2,
                                    double     , tad,
                                    unsigned int, Ls,
                                    std::string, boundary,
                                    std::string, string,
                                    std::string, twist);
};

template <typename FImpl>
class TImprovedStaggeredMILC5D: public Module<ImprovedStaggeredMILC5DPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TImprovedStaggeredMILC5D(const std::string name);
    // destructor
    virtual ~TImprovedStaggeredMILC5D(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ImprovedStaggeredMILC5D, TImprovedStaggeredMILC5D<STAGIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(ImprovedStaggeredMILC5DF, TImprovedStaggeredMILC5D<STAGIMPLF>, MAction);
#endif

/******************************************************************************
 *                 TImprovedStaggeredMILC5D implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TImprovedStaggeredMILC5D<FImpl>::TImprovedStaggeredMILC5D(const std::string name)
: Module<ImprovedStaggeredMILC5DPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TImprovedStaggeredMILC5D<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge, par().gaugefat, par().gaugelong };
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TImprovedStaggeredMILC5D<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TImprovedStaggeredMILC5D<FImpl>::setup(void)
{
    LOG(Message) << "Setting up ImprovedStaggered fermion matrix." << std::endl;
    LOG(Message) << "Using m=" << par().mass << std::endl;
    LOG(Message) << "Using c1=" << par().c1 << std::endl;
    LOG(Message) << "Using c2=" << par().c2 << std::endl;
    LOG(Message) << "Using tadpole u0=" << par().tad << std::endl;
    LOG(Message) << "Using thin links: " << par().gauge << std::endl;
    LOG(Message) << "Using fat links: " << par().gaugefat << std::endl;
    LOG(Message) << "Using long links: " << par().gaugelong << std::endl;
                 
    unsigned int Ls = par().Ls;
    auto &U      = envGet(GaugeField, par().gauge);
    auto &Ufat   = envGet(GaugeField, par().gaugefat);
    auto &Ulong  = envGet(GaugeField, par().gaugelong);
    auto &grid   = *envGetGrid(FermionField);
    auto &gridRb = *envGetRbGrid(FermionField);
    auto &grid5D   = *envGetGrid(FermionField, Ls);
    auto &gridRb5D = *envGetRbGrid(FermionField,Ls);
    typename ImprovedStaggeredFermion5D<FImpl>::ImplParams implParams;
    
    if (!par().boundary.empty())
    {
        implParams.boundary_phases = strToVec<Complex>(par().boundary);
    }
    if (!par().twist.empty())
    {
        implParams.twist_n_2pi_L   = strToVec<Real>(par().twist);
    }

    envCreateDerived(FMat, ImprovedStaggeredFermion5D<FImpl>, getName(), Ls,
                     grid5D, gridRb5D, grid, gridRb,
                     2.*par().mass, 2.*par().c1, 2.*par().c2, par().tad, implParams);

    auto &fmat = envGetDerived(FMat, ImprovedStaggeredFermion5D<FImpl>, getName());
    fmat.ImportGaugeSimple(Ulong, Ufat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TImprovedStaggeredMILC5D<FImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_ImprovedStaggeredMILC5D_hpp_
