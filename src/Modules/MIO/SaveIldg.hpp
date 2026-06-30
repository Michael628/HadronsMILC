/*
 * SaveIldg.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MIO_SaveIldg_hpp_
#define HadronsMILC_MIO_SaveIldg_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Save an ILDG configuration

 gauge         Name of the gauge field object to write. This also covers
               smeared gauge fields such as the fat/long links produced by
               MGauge::HISQSmear, whose GaugeField type is identical to the
               standard PeriodicGimpl gauge field.
 fileStem      Namestem of the file to write the gauge field to. The current
               trajectory number is appended as "<fileStem>.<traj>".
 ensembleLabel Label of the ensemble (written as the ILDG description).
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MIO)

class SaveIldgPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(SaveIldgPar, std::string, gauge, std::string,
                                  fileStem, std::string, ensembleLabel);
};

template <typename GImpl> class TSaveIldg : public Module<SaveIldgPar> {
public:
  GAUGE_TYPE_ALIASES(GImpl, );

public:
  // constructor
  TSaveIldg(const std::string name);
  // destructor
  virtual ~TSaveIldg(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  virtual std::vector<std::string> getOutputFiles(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveIldg, TSaveIldg<GIMPL>, MIO);

/******************************************************************************
 *                       TSaveIldg implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TSaveIldg<GImpl>::TSaveIldg(const std::string name)
    : Module<SaveIldgPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TSaveIldg<GImpl>::getInput(void) {
  return {par().gauge};
}

template <typename GImpl>
std::vector<std::string> TSaveIldg<GImpl>::getOutput(void) {
  return {};
}

template <typename GImpl>
std::vector<std::string> TSaveIldg<GImpl>::getOutputFiles(void) {
  return {par().fileStem + "." + std::to_string(vm().getTrajectory())};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl> void TSaveIldg<GImpl>::setup(void) {}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl> void TSaveIldg<GImpl>::execute(void) {
  std::string fileName =
      par().fileStem + "." + std::to_string(vm().getTrajectory());
  LOG(Message) << "Saving ILDG configuration to file '" << fileName << "'"
               << std::endl;

  auto &U = envGet(GaugeField, par().gauge);
  makeFileDir(fileName, U.Grid());

  std::string description =
      par().ensembleLabel.empty() ? fileName : par().ensembleLabel;

  IldgWriter writer(U.Grid()->IsBoss());
  writer.open(fileName);
  writer.writeConfiguration(U, vm().getTrajectory(), fileName, description);
  writer.close();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MIO_SaveIldg_hpp_
