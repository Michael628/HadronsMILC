/*
 * MesonFieldMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
 * Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef HadronsMILC_MContraction_A2AASlashMesonField_hpp_
#define HadronsMILC_MContraction_A2AASlashMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include "../../A2AVectors.hpp"
#include <Hadrons/EigenPack.hpp>
#include "../../A2AMatrix.hpp"
#include "../../../a2a/A2AWorker.h"
#include "../../../spin/StagGamma.h"

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
*                     All-to-all meson field creation                        *
******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class ASlashMesonFieldMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ASlashMesonFieldMILCPar,
        int, block,
        std::string, lowModes,
        std::string, left,
        std::string, action,
        std::string, right,
        std::string, output,
        SpinTasteParams, spinTaste,
        std::string, EmFunc,
        int,         nEmFields,
        std::string, EmSeedString,
        std::vector<std::string>, mom);
    ASlashMesonFieldMILCPar(): nEmFields(0) {}
};

class ASlashMesonFieldMILCMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ASlashMesonFieldMILCMetadata,
        std::vector<RealF>, momentum,
        std::string, amu_seed,
        int, amu_index);

    ASlashMesonFieldMILCMetadata(): momentum{}, amu_index(-1), amu_seed("N/A") {}
};

template <typename T, typename FImpl>
class ASlashMesonFieldKernelMILC: public A2AKernelMILC<T, typename FImpl::FermionField>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    ASlashMesonFieldKernelMILC(GridBase *grid, const std::vector<ComplexField> &mom, const std::vector<ComplexField> &Amu, int orthogDir)
    {
        _worker = std::make_unique<A2AWorkerLocal<FImpl> >(grid,mom,Amu,orthogDir);
        _vol = 1.;
        for (auto &d: grid->GlobalDimensions())
        {
            _vol *= d;
        }
    }

    virtual ~ASlashMesonFieldKernelMILC(void) = default;
    virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left_e, const FermionField *left_o, 
        const FermionField *right_e, const FermionField *right_o)
    {
        MesonFunction<FImpl>(m, left_e, left_o, right_e, right_o);
    }

    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej, int cbDiv = 1)
    {

        return _vol/cbDiv*(_worker->getFlops())*blockSizei*blockSizej;
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
// return _vol*(12.0*sizeof(T))*blockSizei*blockSizej
// +  _vol*(2.0*sizeof(T)*_mom.size())*blockSizei*blockSizej*_gamma.size();
        return -1.0;
    }

    virtual double kernelTime() {return _worker->_t_kernel;}
    virtual double globalSumTime() {return _worker->_t_gsum;}
private:
template<typename TFImpl, typename ... Args>
    IfNotStag<TFImpl,void> MesonFunction(Args && ... args){
        assert(0);
    }

template<typename TFImpl, typename ... Args>
    IfStag<TFImpl,void> MesonFunction(Args && ... args){
        _worker->StagMesonField(args...);
    }

private:
    double _vol;
    std::unique_ptr<A2AWorkerBase<FImpl> > _worker;
};

template <typename FImpl, typename Pack>
class TASlashMesonFieldMILC : public Module<ASlashMesonFieldMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AMatrixBlockComputationMILC<Complex, 
    FermionField, 
    ASlashMesonFieldMILCMetadata, 
    HADRONS_A2AM_IO_TYPE> Computation;
    typedef ASlashMesonFieldKernelMILC<Complex, FImpl> Kernel;
    typedef PhotonR::GaugeLinkField EmComp;
    typedef PhotonR::GaugeField EmField;

public:
// constructor
    TASlashMesonFieldMILC(const std::string name);
// destructor
    virtual ~TASlashMesonFieldMILC(void){};

// dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
// setup
    virtual void setup(void);
// execution
    virtual void execute(void);
private:
    bool                               _hasPhase{false};
    bool                               _hasAmu{false};
    std::string                        _momphName, _EmName;
    std::vector<std::vector<Real>>     _mom;
};

MODULE_REGISTER(StagA2AASlashMesonField, ARG(TASlashMesonFieldMILC<STAGIMPL,MassShiftEigenPack<STAGIMPL> >), MContraction);

/******************************************************************************
*                  TASlashMesonFieldMILC implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TASlashMesonFieldMILC<FImpl,Pack>::TASlashMesonFieldMILC(const std::string name)
: Module<ASlashMesonFieldMILCPar>(name)
, _momphName(name + "_momph"), _EmName(name + "_em")
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TASlashMesonFieldMILC<FImpl,Pack>::getInput(void)
{
    std::vector<std::string> in = {};
    if (!par().left.empty()) 
        in.push_back(par().left);
    if (!par().right.empty()) 
        in.push_back(par().right);

    if (!par().lowModes.empty()) {
        if (!par().action.empty())
            in.push_back(par().action);
        in.push_back(par().lowModes);
    }

    if (!par().EmFunc.empty()) {
        in.push_back(par().EmFunc);
    }

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TASlashMesonFieldMILC<FImpl,Pack>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TASlashMesonFieldMILC<FImpl,Pack>::setup(void)
{
    _mom.clear();

    for (auto &pstr: par().mom)
    {
        auto p = strToVec<Real>(pstr);

        if (p.size() != env().getNd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                + " components instead of " 
                + std::to_string(env().getNd() - 1));
        }
        _mom.push_back(p);
    }
    int nmom = _mom.size();
    bool allzero = true;
    if (par().mom.size() == 1) {
        for (auto p: _mom[0]) {
            if (p != 0) allzero = false;
        }
    }
    if (allzero) nmom = 0;

    envCache(std::vector<ComplexField>, _momphName, 1, 
        nmom, envGetGrid(ComplexField));

    envCache(std::vector<ComplexField>, _EmName, 1, 
        par().nEmFields, envGetGrid(ComplexField));
    envTmpLat(ComplexField, "coor");
    envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
        env().getNd() - 1, _mom.size(), par().nEmFields, par().block, this);
    envTmp(std::vector<FermionField>, "dummy", 1, 0, envGetGrid(FermionField));
    envTmpLat(PhotonR::GaugeField,"AField");


}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TASlashMesonFieldMILC<FImpl,Pack>::execute(void)
{
    bool hasLowModes = (!par().lowModes.empty());
    bool isCheckerBoarded = (!par().action.empty());

    std::vector<FermionField> *left, *right;

    envGetTmp(std::vector<FermionField>, dummy);
    if (!par().left.empty()) {
        left  = &(envGet(std::vector<FermionField>, par().left));
    } else {
        left = &dummy;
    }
    if (!par().right.empty()) { 
        right = &(envGet(std::vector<FermionField>, par().right));
    } else {
        right = &dummy;
    }

    int nt         = env().getDim().back();
    int N_i        = left->size();
    int N_j        = right->size();

    if (hasLowModes)
    {
        auto &lowModes = envGet(Pack, par().lowModes);
        if (N_j != 0 || N_i == 0) {
            N_i += (isCheckerBoarded?2:1)*lowModes.evec.size();
            LOG(Message) << "N_i: " << N_i << std::endl;
        } else if (N_i != 0 || N_j == 0) {
            N_j += (isCheckerBoarded?2:1)*lowModes.evec.size();
            LOG(Message) << "N_j: " << N_j << std::endl;
        }
    }
    int block      = par().block;

    // if (N_i < block || N_j < block)
    // {
        // HADRONS_ERROR(Range, "blockSize must not exceed size of input vector.");
    // }

    LOG(Message) << "Computing all-to-all meson fields" << std::endl;
    if (hasLowModes)
        LOG(Message) << "Low Modes: '" << par().lowModes << std::endl;

    if (!(par().left.empty() && par().right.empty())) {
        if (!par().left.empty())
            LOG(Message) << "Left: '" << par().left << "'" << std::endl;
        if (!par().right.empty())
            LOG(Message) << "Right: '" << par().right << "'" << std::endl;
    }

    LOG(Message) << "Momenta:" << std::endl;

    for (auto &p: _mom)
    {
        LOG(Message) << "  " << p << std::endl;
    }

    LOG(Message) << " Amu Fields: " << par().nEmFields << std::endl;

    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
    << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) 
    << "/momentum/bilinear)" << std::endl;

    auto &ph = envGet(std::vector<ComplexField>, _momphName);
    auto &Amu = envGet(std::vector<ComplexField>, _EmName);

    startTimer("Momentum phases");
    for (unsigned int j = 0; j < ph.size(); ++j)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(ComplexField, coor);
        ph[j] = Zero();
        for(unsigned int mu = 0; mu < _mom[j].size(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph[j] = ph[j] + (_mom[j][mu]/env().getDim(mu))*coor;
        }
        ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
    }
    stopTimer("Momentum phases");

    if (!_hasAmu && !par().EmFunc.empty())
    {
        StagGamma gamma;

        envGetTmp(EmField,AField);
        envGetTmp(ComplexField, coor);

        startTimer("Stochastic Amu");

        auto &photon = envGet(PhotonR, par().EmFunc);
        auto    &w = envGet(EmComp, "_" + par().EmFunc + "_weight");

        auto &rng = rng4d();
        if (!par().EmSeedString.empty())
            rng.SeedUniqueString(par().EmSeedString);

        auto gamma_vals = StagGamma::ParseSpinTasteString(par().spinTaste.gammas,par().spinTaste.applyG5);
        if (gamma_vals.size() != env().getNd()) {
        HADRONS_ERROR(Argument,"spinTaste parameter must provide 4 gammas for J.A")
        }

        for (unsigned int j = 0; j < par().nEmFields; ++j)
        {

            Amu[j] = Zero();
            photon.StochasticField(AField, rng, w);
            for (unsigned int k = 0; k < env().getNd(); ++k)
            {
                gamma.setSpinTaste(gamma_vals[k]);
                coor = PeekIndex<LorentzIndex>(AField,k);
                gamma.applyPhase(coor,coor);
                Amu[j] += coor;
            }
        }
        _hasAmu = true;
        stopTimer("Stochastic Amu");
    }

    auto AmuIOnameFn = [this](const unsigned int m, const unsigned int g)
    {
        std::stringstream ss;

        ss << par().EmSeedString << "_" << int(g) << "_";

        for (unsigned int mu = 0; mu < _mom[m].size(); ++mu)
        {
            ss << _mom[m][mu] << ((mu == _mom[m].size() - 1) ? "" : "_");
        }

        return ss.str();
    };

    auto AmuFilenameFn = [this, &AmuIOnameFn](const unsigned int m, const unsigned int g)
    {
        return par().output + "." + std::to_string(vm().getTrajectory()) 
        + "/" + AmuIOnameFn(m, g) + ".h5";
    };

    auto AmuMetadataFn = [this](const unsigned int m, const unsigned int g)
    {
        ASlashMesonFieldMILCMetadata md;

        for (auto pmu: _mom[m])
        {
            md.momentum.push_back(pmu);
        }

        md.amu_seed = par().EmSeedString.empty() ? "<blank> see output file for seed string!": par().EmSeedString;
        md.amu_index = g;

        return md;
    };

    envGetTmp(Computation, computation);

    int orthogDir = env().getNd() - 1;

    if(hasLowModes) {
        auto &lowModes = envGet(Pack, par().lowModes);

        if (isCheckerBoarded) {

            auto &action      = envGet(FMat, par().action);
            std::function<void(int)> swapEvecCheckerFn = [this,&action, &lowModes](int index)
            {
                ComplexD eval_D = ComplexD(0.0,lowModes.eval[index].imag());
                int cb = lowModes.evec[index].Checkerboard();
                int cbNeg = (cb==Even) ? Odd : Even;

                FermionField temp(lowModes.evec[index].Grid());
                temp.Checkerboard() = cbNeg;
                action.Meooe(lowModes.evec[index], temp);
                lowModes.evec[index].Checkerboard() = cbNeg;
                lowModes.evec[index] = (1.0/eval_D) * temp;
            };

            if (par().nEmFields > 0) {
                Kernel kernel(envGetGrid(FermionField), ph, Amu, orthogDir);
                computation.execute(*left, *right, kernel, AmuIOnameFn, AmuFilenameFn, AmuMetadataFn, &lowModes.evec, lowModes.eval, &swapEvecCheckerFn);
            }
        } else{
            if (par().nEmFields > 0) {
                Kernel kernel(envGetGrid(FermionField), ph, Amu, orthogDir);
                computation.execute(*left, *right, kernel, AmuIOnameFn, AmuFilenameFn, AmuMetadataFn, &lowModes.evec, lowModes.eval);
            }
        }
    } else {
        if (par().nEmFields > 0) {
            Kernel kernel(envGetGrid(FermionField), ph, Amu, orthogDir);
            computation.execute(*left, *right, kernel, AmuIOnameFn, AmuFilenameFn, AmuMetadataFn);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_ASlashMesonFieldMILC_hpp_
