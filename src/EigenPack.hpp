/*
 * EigenPack.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
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
#ifndef HadronsMILC_EigenPack_hpp_
#define HadronsMILC_EigenPack_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

template <typename Field> class AdaptorEigenPackMILC {
public:
  AdaptorEigenPackMILC(std::vector<Field> &_evec,
                       const std::vector<RealD> &_eval, RealD _mass = 0.0)
      : evec(_evec), mass(_mass) {
    // Store shifted M eigenvalues instead of massless M^dagM eigenvalues
    eval.resize(_eval.size(), 0.0);

    Real m = 2 * _mass;

    for (int i = 0; i < eval.size(); i++) {
      eval[i] = ComplexD(m, sqrt(_eval[i]));
    }
  }

public:
  std::vector<Field> &evec;
  std::vector<ComplexD> eval;
  RealD mass;
};

template <typename FImpl>
using MassShiftEigenPack = AdaptorEigenPackMILC<typename FImpl::FermionField>;

END_HADRONS_NAMESPACE

#endif // Hadrons_EigenPack_hpp_
