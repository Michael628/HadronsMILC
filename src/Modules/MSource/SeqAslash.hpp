/*
 * SeqAslashMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
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

#ifndef HadronsMILC_MSource_SeqAslash_hpp_
#define HadronsMILC_MSource_SeqAslash_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * i * A_mu g_mu * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - emField: input photon field (string)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                        Sequential Aslash source                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqAslashMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqAslashMILCPar,
                                    std::string,    q,
                                    unsigned int,   tA,
                                    unsigned int,   tB,
                                    std::string,    gammas,
                                    std::string,    emField,
                                    std::string,    mom);
};
// gammas: name of an MFermion::SpinTaste module publishing exactly 4
//         std::vector<StagGamma> entries -- the directional J.A loop
//         contracts entry mu with PeekIndex<LorentzIndex>(A, mu)

template <typename FImpl>
class TSeqAslashMILC: public Module<SeqAslashMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    typedef PhotonR::GaugeField     EmField;
public:
    // constructor
    TSeqAslashMILC(const std::string name);
    // destructor
    virtual ~TSeqAslashMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeSource(PropagatorField &src, const PropagatorField &q);
};

MODULE_REGISTER_TMP(StagSeqAslash, TSeqAslashMILC<STAGIMPL>, MSource);

/******************************************************************************
 *                         TSeqAslashMILC implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqAslashMILC<FImpl>::TSeqAslashMILC(const std::string name)
: Module<SeqAslashMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqAslashMILC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().emField};

    if (par().gammas.empty())
    {
        HADRONS_ERROR(Argument, "SeqAslash requires the 'gammas' "
                      "SpinTaste module name");
    }
    in.push_back(par().gammas);

    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqAslashMILC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqAslashMILC<FImpl>::setup(void)
{
    if (par().gammas.empty())
    {
        HADRONS_ERROR(Argument, "SeqAslash requires the 'gammas' "
                      "SpinTaste module name");
    }

    envTmpLat(PropagatorField, "field");

    if (envHasType(PropagatorField, par().q))
    {
        envCreateLat(PropagatorField, getName());
    }
    else if (envHasType(std::vector<PropagatorField>, par().q))
    {
        auto &q = envGet(std::vector<PropagatorField>, par().q);

        envCreate(std::vector<PropagatorField>, getName(), 1, q.size(),
                envGetGrid(PropagatorField));
    }
    else
    {
        HADRONS_ERROR_REF(ObjectType, "object '" + par().q 
                          + "' has an incompatible type ("
                          + env().getObjectType(par().q)
                          + ")", env().getObjectAddress(par().q))
    }
    envTmp(Lattice<iScalar<vInteger>>, "t", 1, envGetGrid(LatticeComplex));
    envTmpLat(LatticeComplex, "ph");
    envTmpLat(LatticeComplex, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqAslashMILC<FImpl>::makeSource(PropagatorField &src, 
                                   const PropagatorField &q)
{
    envGetTmp(LatticeComplex,ph);
    envGetTmp(Lattice<iScalar<vInteger>>, t);
    auto &stoch_photon = envGet(EmField, par().emField);

    Complex           i(0.0,1.0);
    std::vector<Real> p;

    envGetTmp(LatticeComplex, coor);
    p  = strToVec<Real>(par().mom);
    ph = Zero();
    for(unsigned int mu = 0; mu < env().getNd(); mu++)
    {
        LatticeCoordinate(coor, mu);
        ph = ph + (p[mu]/env().getDim(mu))*coor;
    }
    ph = exp((Real)(2*M_PI)*i*ph);
    LatticeCoordinate(t, Tp);
    
    Complex ci(0.0,1.0);
    src = Zero();

    const auto &gammas = envGet(std::vector<StagGamma>, par().gammas);
    if (gammas.size() != 4)
    {
        HADRONS_ERROR(Argument, "gammas module '" + par().gammas +
                      "' must publish exactly 4 gammas for J.A (got " +
                      std::to_string(gammas.size()) + ")");
    }

    envGetTmp(PropagatorField,field);
    field = Zero();

    for (unsigned int mu = 0; mu < gammas.size(); mu++)
    {
        gammas[mu](field, q);

        src = src + where((t >= par().tA) and (t <= par().tB),
                          ci*PeekIndex<LorentzIndex>(stoch_photon, mu) *(ph*field), 0.*field);
    }
}

template <typename FImpl>
void TSeqAslashMILC<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating Aslash sequential source(s) at t= " << par().tA 
		             << " using gammas module '" << par().gammas << "' and photon field '" << par().emField 
                     << "'" << std::endl; 
    }
    else
    {
        LOG(Message) << "Generating Aslash sequential source(s) for "
                     << par().tA << " <= t <= " << par().tB 
		             << " using gammas module '" << par().gammas << "' and photon field '" << par().emField 
                     << "'" << std::endl;
    }
    
    if (envHasType(PropagatorField, par().q))
    {
        auto  &src = envGet(PropagatorField, getName()); 
        auto  &q   = envGet(PropagatorField, par().q);

        LOG(Message) << "Using propagator '" << par().q << "'" << std::endl;
        makeSource(src, q);
    }
    else
    {
        auto  &src = envGet(std::vector<PropagatorField>, getName()); 
        auto  &q   = envGet(std::vector<PropagatorField>, par().q);

        for (unsigned int i = 0; i < q.size(); ++i)
        {
            LOG(Message) << "Using element " << i << " of propagator vector '" 
                         << par().q << "'" << std::endl;
            makeSource(src[i], q[i]);
        }
    }

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SeqAslash_hpp_
