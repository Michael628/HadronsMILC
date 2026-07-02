/*
 * LoadMilc.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2024
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Michael Lynch <michaellynch628@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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
#ifndef HadronsMILC_MIO_LoadMilc_hpp_
#define HadronsMILC_MIO_LoadMilc_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Load a MILC v5 gauge configuration

 file                    Namestem of the MILC gauge file to read. The current
                         trajectory number is appended as "<file>.<traj>".
 exitOnChecksumMismatch  If true, a MILC sum29/sum31 checksum mismatch is a
                         fatal error; if false (default) it is only logged.
                         Mirrors MILC's own non-fatal checksum warning and
                         NerscIO::exitOnReadPlaquetteMismatch().

 The MILC v5 gauge-configuration reader itself (MilcHeader/MilcIO) lives in
 GridMilc/io/MilcIO.h; this module is just the Hadrons wrapper around it.
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 * LoadMilc module parameters
 ******************************************************************************/
class LoadMilcPar : Serializable {
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadMilcPar,
                                    std::string, file,
                                    bool,        exitOnChecksumMismatch);
    // warn by default on a sum29/sum31 mismatch (mirrors MILC's non-fatal
    // warning); opt-in to a fatal error via exitOnChecksumMismatch=true.
    LoadMilcPar(void)
        : exitOnChecksumMismatch(false) {}
};

/******************************************************************************
 * TLoadMilc module
 ******************************************************************************/
template <typename GImpl> class TLoadMilc : public Module<LoadMilcPar> {
public:
    GAUGE_TYPE_ALIASES(GImpl, );

public:
    // constructor
    TLoadMilc(const std::string name);
    // destructor
    virtual ~TLoadMilc(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadMilc, TLoadMilc<GIMPL>, MIO);

/******************************************************************************
 *                       TLoadMilc implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLoadMilc<GImpl>::TLoadMilc(const std::string name)
    : Module<LoadMilcPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TLoadMilc<GImpl>::getInput(void) {
    return {};
}

template <typename GImpl>
std::vector<std::string> TLoadMilc<GImpl>::getOutput(void) {
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl> void TLoadMilc<GImpl>::setup(void) {
    envCreateLat(GaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl> void TLoadMilc<GImpl>::execute(void) {
    MilcHeader  header;
    std::string fileName
        = par().file + "." + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading MILC configuration from file '" << fileName << "'"
                 << std::endl;

    auto &U = envGet(GaugeField, getName());
    MilcIO::readConfiguration(U, header, fileName, par().exitOnChecksumMismatch);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MIO_LoadMilc_hpp_
