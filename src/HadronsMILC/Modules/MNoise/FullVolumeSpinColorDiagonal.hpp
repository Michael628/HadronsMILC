/*
 * FullVolumeSpinColorDiagonal.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
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
#ifndef HadronsMILC_MNoise_FullVolumeSpinColorDiagonal_hpp_
#define HadronsMILC_MNoise_FullVolumeSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include "../../DilutedNoise.hpp"

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             Generate full volume spin-color diagonal noise                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class FullVolumeSpinColorDiagonalMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FullVolumeSpinColorDiagonalMILCPar,
                                    unsigned int, nsrc);
};

template <typename FImpl>
class TFullVolumeSpinColorDiagonalMILC: public Module<FullVolumeSpinColorDiagonalMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TFullVolumeSpinColorDiagonalMILC(const std::string name);
    // destructor
    virtual ~TFullVolumeSpinColorDiagonalMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StagFullVolumeSpinColorDiagonal, TFullVolumeSpinColorDiagonalMILC<STAGIMPL>, MNoise);

/******************************************************************************
 *              TFullVolumeSpinColorDiagonalMILC implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TFullVolumeSpinColorDiagonalMILC<FImpl>::TFullVolumeSpinColorDiagonalMILC(const std::string name)
: Module<FullVolumeSpinColorDiagonalMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TFullVolumeSpinColorDiagonalMILC<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TFullVolumeSpinColorDiagonalMILC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_vec"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TFullVolumeSpinColorDiagonalMILC<FImpl>::setup(void)
{
    envCreateDerived(SpinColorDiagonalNoiseMILC<FImpl>, 
                     FullVolumeNoiseMILC<FImpl>,
                     getName(), 1, envGetGrid(FermionField), par().nsrc);

    envCreate(std::vector<FermionField>, getName() + "_vec", 1, 0, envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TFullVolumeSpinColorDiagonalMILC<FImpl>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoiseMILC<FImpl>, getName());
    LOG(Message) << "Generating full volume, spin-color diagonal noise" << std::endl;
    noise.generateNoise(rng4d());

    auto &noisevec = envGet(std::vector<FermionField>,getName()+"_vec");

    int nferm = noise.fermSize();
    int nsc   = (int)(nferm/(noise.size()));

    noisevec.resize(nferm,envGetGrid(FermionField));
    for (int i=0;i<nferm;i++) {
        noisevec[i] = noise.getFerm(i);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MNoise_FullVolumeSpinColorDiagonal_hpp_
