/*
 * RandomWallMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
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

#ifndef HadronsMILC_MSource_RandomWall_hpp_
#define HadronsMILC_MSource_RandomWall_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/*

 Random Wall source
 -----------------------------

 * options:
 - tStep:     step between source timeslices (integer)
 - t0:        first source timeslice, must be < tStep (integer)
 - nSrc:      number of noise sources (integer)
 - reuset0:   reuse noise vectors at t=t0 and shift in time ("true"/"false")
 - colorDiag: store propagator fields if true, fermion fields if false
 - noise:     name of an external spin-color diagonal noise module, either
              time-diluted (MNoise::TimeDilutedSpinColorDiagonal) or
              full-volume (MNoise::FullVolumeSpinColorDiagonal); if empty,
              time-diluted noise is generated internally. Full-volume noise
              is masked to a single time slice per source (zeros elsewhere),
              matching the internal-noise output bit-for-bit when the noise
              module name in this run equals the RandomWall module name in
              the internal-noise run (same RNG stream) with nsrc >= nSrc.

 */

/******************************************************************************
 *                         Random Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class RandomWallMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomWallMILCPar,
                                    unsigned int, tStep,
                                    unsigned int, t0,
                                    unsigned int, nSrc,
                                    std::string,  reuset0,
                                    bool,         colorDiag,
                                    std::string,  noise);
};

template <typename FImpl>
class TRandomWallMILC: public Module<RandomWallMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    HADRONS_DEFINE_setProp_setFerm(FImpl);
public:
    // constructor
    TRandomWallMILC(const std::string name);
    // destructor
    virtual ~TRandomWallMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool reuset0_ = false;
};

MODULE_REGISTER_TMP(StagRandomWall, TRandomWallMILC<STAGIMPL>, MSource);

/******************************************************************************
 *                 TRandomWallMILC implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRandomWallMILC<FImpl>::TRandomWallMILC(const std::string name)
: Module<RandomWallMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRandomWallMILC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    if (!par().noise.empty()) {
        in.push_back(par().noise);
    }
    return in;
}

template <typename FImpl>
std::vector<std::string> TRandomWallMILC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_shift"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomWallMILC<FImpl>::setup(void)
{
    if (par().noise.empty()) {
        envTmp(TimeDilutedNoiseMILC<FImpl>, "noise", 1, envGetGrid(FermionField), par().nSrc);
    }
    envTmpLat(PropagatorField, "shiftedField");
    envTmp(Lattice<iScalar<vInteger>>, "t", 1, envGetGrid(PropagatorField));
    envTmpLat(FermionField,"ferm");

    // Scalar collapse: with a single output vector (nSrc == 1 and a
    // single time slice, nt/min(tStep,nt) == 1) publish the module name
    // as a SCALAR field instead of a length-1 std::vector — GaugeProp
    // outputs and guesses mirror the source container type, so a
    // vector-of-1 wall would force vector-typed guesses downstream
    // while scalar producers (StagLMAMesonFieldProp) offer scalar
    // guesses only. Both paths kept (a3f6e74 precedent); _shift stays a
    // std::vector<Integer> (Meson.hpp envGets it unconditionally).
    if (par().tStep < 1) {
        HADRONS_ERROR(Logic, "Parameter tStep must be >= 1 (got " +
                            std::to_string(par().tStep) + ")");
    }
    const int nt       = static_cast<int>(env().getDim().back());
    const int tStep    = par().tStep;
    const int nSources = par().nSrc;
    const int nSlices  = nt/std::min(tStep,nt);
    const bool scalar  = (nSources*nSlices == 1);

    if (par().colorDiag) {
        if (scalar) {
            envCreate(PropagatorField, getName(), 1, envGetGrid(PropagatorField));
        } else {
            envCreate(std::vector<PropagatorField>, getName(), 1, 0, envGetGrid(PropagatorField));
        }
    } else {
        if (scalar) {
            envCreate(FermionField, getName(), 1, envGetGrid(FermionField));
        } else {
            envCreate(std::vector<FermionField>, getName(), 1, 0, envGetGrid(FermionField));
        }
    }

    envCreate(std::vector<Integer>, getName()+"_shift", 1, 0, 0);

    if (!par().reuset0.empty()) {
        if (!(std::istringstream(par().reuset0) >> std::boolalpha >> reuset0_)) {
            HADRONS_ERROR(Logic,"parameter reuset0='" + par().reuset0 + "' must be 'true' or 'false'");
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomWallMILC<FImpl>::execute(void)
{    
    SpinColorDiagonalNoiseMILC<FImpl> *noise;
    bool fullVolumeNoise = false;
    if (par().noise.empty()) {
        noise = env().template getObject<TimeDilutedNoiseMILC<FImpl> >(getName() + "_tmp_noise");
        LOG(Message) << "Generating " << par().nSrc << " time-diluted, spin-color diagonal noise sources at every " << par().tStep << " time step(s)" << std::endl;
        noise->generateNoise(rng4d());
    } else {
        noise = env().template getObject<SpinColorDiagonalNoiseMILC<FImpl> >(par().noise);
        if (static_cast<int>(par().nSrc) > noise->size()) {
            HADRONS_ERROR(Logic,"external noise '" + par().noise + "' provides "
                         + std::to_string(noise->size()) + " source(s),"
                         + " cannot feed nSrc = " + std::to_string(par().nSrc));
        }
        fullVolumeNoise = (dynamic_cast<FullVolumeNoiseMILC<FImpl> *>(noise) != nullptr);
        if (fullVolumeNoise) {
            LOG(Message) << "Reading full-volume, spin-color diagonal noise from '"
                         << par().noise << "', keeping single time slices only" << std::endl;
        }
    }
    envGetTmp(PropagatorField,shiftedField);
    envGetTmp(FermionField,ferm);
    envGetTmp(Lattice<iScalar<vInteger>>, t);

    bool colorDiag = par().colorDiag;

    auto &time_shift = envGet(std::vector<Integer>,getName()+"_shift");

    int nt       = envGetGrid(PropagatorField)->GlobalDimensions()[Tp];
    int tStep    = par().tStep;
    int nSources = par().nSrc;
    int t0       = par().t0;

    if (t0 >= tStep) {
        HADRONS_ERROR(Logic,"Parameter t0 >= tStep");
    }

    int nSlices = nt/std::min(tStep,nt);
    int nVecs   = nSources*nSlices;

    time_shift.resize(nVecs,0);

    if (fullVolumeNoise) {
        LatticeCoordinate(t, Tp);
    }

    // time-diluted noise: propagator index is i*nt + t;
    // full-volume noise:  propagator index is i, masked to the t time slice
    auto getNoiseProp = [&](int i, int tSlice) -> PropagatorField & {
        if (!fullVolumeNoise) {
            return noise->getProp(i*nt + tSlice);
        } else {
            PropagatorField &src = noise->getProp(i);
            shiftedField = where((t == tSlice), src, 0.*src);
            return shiftedField;
        }
    };

    if (reuset0_) {
        LOG(Message) << "Reusing noise vectors at t=0 and shifting by " << par().tStep << std::endl;
    }


    if (colorDiag) {
        if (nVecs == 1) {
            // scalar collapse (see setup): single source, single slice.
            // j == 0 in either reuset0 mode — both compute
            // getNoiseProp(0, t0) — so the branch is unconditional.
            auto &noisevec = envGet(PropagatorField,getName());
            noisevec = getNoiseProp(0, t0);
            time_shift[0] = t0;
            return;
        }
        auto &noisevec = envGet(std::vector<PropagatorField>,getName());
        noisevec.resize(nVecs,envGetGrid(PropagatorField));

        for (int i=0;i<nSources;i++) {
            if (reuset0_) {
                noisevec[i*nSlices] = getNoiseProp(i, t0);
            }
            for (int j=0;j<nSlices;j++) {
                int idx = i*nSlices+j;
                if (!reuset0_) {
                    noisevec[idx] = getNoiseProp(i, j*tStep+t0);
                } else {
                    if (j != 0) {
                        noisevec[idx] = Cshift(noisevec[idx-1],Tp,tStep);
                    }
                }
                time_shift[idx] = j*tStep+t0;
            }
        }
    } else {
        if (nVecs == 1) {
            // scalar collapse (see setup): single source, single slice;
            // color-summed exactly as the vector path's element 0
            auto &noisevec = envGet(FermionField,getName());
            noisevec = Zero();
            PropagatorField &srcProp = getNoiseProp(0, t0);
            for (int k=0;k<FImpl::Dimension;k++) {
                setFerm(ferm,srcProp,k);
                noisevec += ferm;
            }
            time_shift[0] = t0;
            return;
        }
        auto &noisevec = envGet(std::vector<FermionField>,getName());
        noisevec.resize(nVecs,envGetGrid(FermionField));

        for (int i=0;i<nSources;i++) {
            if (reuset0_) {
                PropagatorField &srcProp = getNoiseProp(i, t0);
                noisevec[i*nSlices] = Zero();
                for (int k=0;k<FImpl::Dimension;k++) {
                    setFerm(ferm,srcProp,k);
                    noisevec[i*nSlices] += ferm;
                }
            }
            for (int j=0;j<nSlices;j++) {
                int idx = i*nSlices+j;
                noisevec[idx] = Zero();
                if (!reuset0_) {
                    PropagatorField &srcProp = getNoiseProp(i, j*tStep+t0);
                    for (int k=0;k<FImpl::Dimension;k++) {
                        setFerm(ferm,srcProp,k);
                        noisevec[idx] += ferm;
                    }
                } else {
                    if (j != 0) {
                        noisevec[idx] = Cshift(noisevec[idx-1],Tp,tStep);
                    }
                }
                time_shift[idx] = j*tStep+t0;
            }
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_RandomWall_hpp_
