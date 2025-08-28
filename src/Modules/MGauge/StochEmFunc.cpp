/*
 * StochEm.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#include "StochEmFunc.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

/******************************************************************************
*                  TStochEmFunc implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TStochEmFunc::TStochEmFunc(const std::string name)
: Module<StochEmFuncPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TStochEmFunc::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TStochEmFunc::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TStochEmFunc::setup(void)
{
    weightDone_ = env().hasCreatedObject("_" + getName() + "_weight");
    envCacheLat(EmComp, "_" + getName() + "_weight");

    std::vector<Real> improvements = strToVec<Real>(par().improvement);
    envCreate(PhotonR, getName(), 1, envGetGrid(EmField), par().gauge, par().zmScheme, improvements);
}

// execution ///////////////////////////////////////////////////////////////////
void TStochEmFunc::execute(void)
{
    LOG(Message) << "Generating stochastic EM potential..." << std::endl;

    auto    &w = envGet(EmComp, "_" + getName() + "_weight");
    auto    &photon = envGet(PhotonR, getName());
    if (!weightDone_)
    {
        LOG(Message) << "Caching stochastic EM potential weight (gauge: "
                     << par().gauge << ", zero-mode scheme: "
                     << par().zmScheme << ")..." << std::endl;
        photon.StochasticWeight(w);
    }
}
