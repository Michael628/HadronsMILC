/*
 * StochEm.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: James Harrison <jch1g10@soton.ac.uk>
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
#ifndef HadronsMILC_MGauge_StochEmFunc_hpp_
#define HadronsMILC_MGauge_StochEmFunc_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         StochEmFunc                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class StochEmFuncPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochEmFuncPar,
                                    PhotonR::Gauge,    gauge,
                                    PhotonR::ZmScheme, zmScheme,
                                    std::string,       improvement);
};

class TStochEmFunc: public Module<StochEmFuncPar>
{
public:
    typedef PhotonR::GaugeField     EmField;
    typedef PhotonR::GaugeLinkField EmComp;
public:
    // constructor
    TStochEmFunc(const std::string name);
    // destructor
    virtual ~TStochEmFunc(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool    weightDone_;
};

MODULE_REGISTER(StochEmFunc, TStochEmFunc, MGauge);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_StochEm_hpp_
