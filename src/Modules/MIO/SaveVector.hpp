/*
 * SaveVector.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
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
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MIO_SaveVector_hpp_
#define Hadrons_MIO_SaveVector_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Module to save a single field to disk                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class SaveVectorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveVectorPar,
                                    std::string, field,
                                    std::string, output,
                                    bool, multiFile);
};

template <typename Field>
class TSaveVector: public Module<SaveVectorPar>
{
public:
    // constructor
    TSaveVector(const std::string name);
    // destructor
    virtual ~TSaveVector(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveStagVector, TSaveVector<STAGIMPL::FermionField>, MIO);
MODULE_REGISTER_TMP(SaveColourMatrixFieldVector, TSaveVector<GIMPL::GaugeLinkField>, MIO);

/******************************************************************************
 *                 TSaveVector implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TSaveVector<Field>::TSaveVector(const std::string name)
: Module<SaveVectorPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TSaveVector<Field>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename Field>
std::vector<std::string> TSaveVector<Field>::getOutput(void)
{
    std::vector<std::string> out;
    
    return out;
}

template <typename Field>
std::vector<std::string> TSaveVector<Field>::getOutputFiles(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TSaveVector<Field>::setup(void)
{}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TSaveVector<Field>::execute(void)
{
    auto &v     = envGet(std::vector<Field>, par().field);

    A2AVectorsIo::write(par().output, v, par().multiFile, vm().getTrajectory());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_SaveVector_hpp_a
