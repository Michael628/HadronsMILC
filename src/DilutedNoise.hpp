/*
 * DilutedNoiseMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
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
#ifndef HadronsMILC_DilutedNoise_hpp_
#define HadronsMILC_DilutedNoise_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Grid/util/Sha.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Abstract container for spin color diagonal noise              *
 ******************************************************************************/
template <typename FImpl>
class SpinColorDiagonalNoiseMILC
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    SpinColorDiagonalNoiseMILC(GridCartesian *g);
    SpinColorDiagonalNoiseMILC(GridCartesian *g, const int nNoise);
    virtual ~SpinColorDiagonalNoiseMILC(void) = default;
    // access
    FermionField &                      getFerm(const int i);
    PropagatorField &                   getProp(const int i);
    int                                 size(void) const;

    void                                generateNoise(GridParallelRNG &rng);
    void                                resize(const int nNoise);
    virtual int                         dilutionSize(void) const = 0;

    GridCartesian *                     getGrid(void) const      { return grid_;}
    int                                 fermSize(void) const     { return dilutionSize(); }
    std::vector<LatticeComplex> &       getNoise(void)           { return noise_; }
    const std::vector<LatticeComplex> & getNoise(void) const     { return noise_; }
    template <typename T = FImpl>
    IfNotStag<T,int>                    getNsc(void) const       { return Ns*FImpl::Dimension; }
    template <typename T = FImpl>
    IfStag<T,int>                       getNsc(void) const       { return FImpl::Dimension; }
private:
    template <typename T = FImpl>
    IfNotStag<T,void>                   setFerm(const int i);

    template <typename T = FImpl>
    IfStag<T,void>                      setFerm(const int i);
    virtual void                        setProp(const int i) = 0;
protected:
    int                                 getNd(void) const;
    FermionField & getFerm(void) 
    { 
        return ferm_; 
    }
    PropagatorField & getProp(void) 
    { 
        return prop_; 
    }
    template <typename T = FImpl>
    IfNotStag<T,void>                   setPropagator(LatticeComplex* eta);
    template <typename T = FImpl>
    IfStag<T,void>                      setPropagator(LatticeComplex* eta);
protected:
    int                            currentProp_;
private:
    FermionField                   ferm_;
    GridCartesian *                grid_;
    std::vector<LatticeComplex>    noise_;
    PropagatorField                prop_;
};

template <typename FImpl>
class TimeDilutedNoiseMILC: public SpinColorDiagonalNoiseMILC<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    TimeDilutedNoiseMILC(GridCartesian *g);
    TimeDilutedNoiseMILC(GridCartesian *g, const int nNoise);
    virtual ~TimeDilutedNoiseMILC(void) = default;
    int dilutionSize(void) const;
    Lattice<iScalar<vInteger>> & getTLat() { return tLat_; }
protected:
    Lattice<iScalar<vInteger>> tLat_;
private:
    void setProp(const int i);
};

template <typename FImpl>
class FullVolumeNoiseMILC: public SpinColorDiagonalNoiseMILC<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    FullVolumeNoiseMILC(GridCartesian *g, const int nNoise);
    virtual ~FullVolumeNoiseMILC(void) = default;
    int dilutionSize(void) const;
private:
    void setProp(const int i);
};

template <typename FImpl>
class CheckerboardNoiseMILC: public SpinColorDiagonalNoiseMILC<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    CheckerboardNoiseMILC(GridCartesian *g, const int nNoise, const int nSparse);
    virtual ~CheckerboardNoiseMILC(void) = default;
    int dilutionSize(void) const;
private:
    void setProp(const int i);
    int nSparse_, nSrc_ec_;
    LatticeInteger coor_, coorTot_;
};

template <typename FImpl>
class SparseNoiseMILC: public SpinColorDiagonalNoiseMILC<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    SparseNoiseMILC(GridCartesian *g, const int nNoise, const int nSparseL, const int nSparseT);
    virtual ~SparseNoiseMILC(void) = default;
    int dilutionSize(void) const;
private:
    void setProp(const int i);
    int nSparseL_,nSparseT_;
    LatticeInteger coor_;
};
/******************************************************************************
 *               SpinColorDiagonalNoiseMILC template implementation               *
 ******************************************************************************/
template <typename FImpl>
SpinColorDiagonalNoiseMILC<FImpl>::SpinColorDiagonalNoiseMILC(GridCartesian *g)
: grid_(g), ferm_(g), prop_(g)
{}

template <typename FImpl>
SpinColorDiagonalNoiseMILC<FImpl>::SpinColorDiagonalNoiseMILC(GridCartesian *g,
                                                      const int nNoise)
: SpinColorDiagonalNoiseMILC(g)
{
    resize(nNoise);
}

template <typename FImpl>
template <typename T>
IfNotStag<T,void> SpinColorDiagonalNoiseMILC<FImpl>::setFerm(const int i)
{
    int nc  = FImpl::Dimension;
    std::div_t divs;
    divs = std::div(i, nc);

    PropToFerm<FImpl>(ferm_, prop_, divs.quot, divs.rem);
}

template <typename FImpl>
template <typename T>
IfStag<T,void> SpinColorDiagonalNoiseMILC<FImpl>::setFerm(const int i)
{
    int nc  = FImpl::Dimension;
    std::div_t divs;
    divs = std::div(i, nc);

    PropToFerm<FImpl>(ferm_, prop_, divs.rem);
}

template <typename FImpl>
typename SpinColorDiagonalNoiseMILC<FImpl>::FermionField & 
SpinColorDiagonalNoiseMILC<FImpl>::getFerm(const int i)
{
    auto nsc   = this->getNsc();
    std::div_t divs;
    divs = std::div(i, nsc);
    setProp(divs.quot);
    setFerm(divs.rem);
    return getFerm();
}

template <typename FImpl>
template <typename T>
IfStag<T,void> SpinColorDiagonalNoiseMILC<FImpl>::setPropagator(LatticeComplex * eta)
{
    prop_ = Zero();
    for (int i=0; i<this->getNsc();i++) {
        pokeColour(prop_,eta[i],i,i);
    }
}

template <typename FImpl>
template <typename T>
IfNotStag<T,void> SpinColorDiagonalNoiseMILC<FImpl>::setPropagator(LatticeComplex * eta)
{
    int nc  = FImpl::Dimension;
    std::div_t divs;
    divs = std::div(nc, this->getNsc());

    prop_ = Zero();
    for (int i=0; i<divs.quot; i++) {
        auto propTmp = peekSpin(prop_,i,i);
        for (int j=0; j<divs.rem; j++) {
            pokeColour(propTmp,eta[divs.rem*i+j],j,j);
        }
        pokeSpin(prop_,propTmp,i,i);
    }
}

template <typename FImpl>
typename SpinColorDiagonalNoiseMILC<FImpl>::PropagatorField & 
SpinColorDiagonalNoiseMILC<FImpl>::getProp(const int i)
{
    setProp(i);
    return getProp();
}

template <typename FImpl>
int SpinColorDiagonalNoiseMILC<FImpl>::size(void) const
{
    return noise_.size()/this->getNsc();
}

template <typename FImpl>
int SpinColorDiagonalNoiseMILC<FImpl>::getNd(void) const
{
    return grid_->GlobalDimensions().size();
}

template <typename FImpl>
void SpinColorDiagonalNoiseMILC<FImpl>::resize(const int nNoise)
{
    noise_.resize(this->getNsc()*nNoise, grid_);
}

template <typename FImpl>
void SpinColorDiagonalNoiseMILC<FImpl>::generateNoise(GridParallelRNG &rng)
{
    Complex        shift(1., 1.);
    LatticeComplex eta(grid_);

    eta = Zero();
    
    for (int n = 0; n < noise_.size(); ++n)
    {
        bernoulli(rng, eta);
        eta = (2.*eta - shift)*(1./::sqrt(2.));
        noise_[n] = eta;
    }
}

/******************************************************************************
 *                  TimeDilutedNoiseMILC template implementation                  *
 ******************************************************************************/
template <typename FImpl>
TimeDilutedNoiseMILC<FImpl>::
TimeDilutedNoiseMILC(GridCartesian *g, int nNoise)
: SpinColorDiagonalNoiseMILC<FImpl>(g, nNoise), tLat_(g)
{}

template <typename FImpl>
int TimeDilutedNoiseMILC<FImpl>::dilutionSize() const
{
    auto nt = this->getGrid()->GlobalDimensions()[Tp];
    return nt*this->getNsc()*this->size();
}

template <typename FImpl>
void TimeDilutedNoiseMILC<FImpl>::setProp(const int i)
{
    auto noise = this->getNoise();
    auto nd    = this->getNd();
    auto nt    = this->getGrid()->GlobalDimensions()[Tp];
    auto nsc   = this->getNsc();

    std::vector<LatticeComplex> eta(nsc,this->getGrid());
    LatticeCoordinate(tLat_, nd - 1);

    std::div_t divs = std::div(i, nt);
    int t = divs.rem;

    for (int j=0;j<nsc;j++)
        eta[j] = where((tLat_ == t), noise[divs.quot*nsc+j], 0.*noise[divs.quot*nsc+j]);

    this->setPropagator(&(eta[0]));
}

/******************************************************************************
 *                   FullVolumeNoiseMILC template implementation                  *
 ******************************************************************************/
template <typename FImpl>
FullVolumeNoiseMILC<FImpl>::
FullVolumeNoiseMILC(GridCartesian *g, int nNoise)
: SpinColorDiagonalNoiseMILC<FImpl>(g, nNoise)
{}

template <typename FImpl>
int FullVolumeNoiseMILC<FImpl>::dilutionSize() const
{
    return this->getNsc()*this->size();
}

template <typename FImpl>
void FullVolumeNoiseMILC<FImpl>::setProp(const int i)
{
    auto noise = this->getNoise();
    this->setPropagator(&noise[i*this->getNsc()]);
}

/******************************************************************************
 *                CheckerboardNoiseMILC template implementation                   *
 ******************************************************************************/
template <typename FImpl>
CheckerboardNoiseMILC<FImpl>::
CheckerboardNoiseMILC(GridCartesian *g, int nNoise, int nSparse)
: SpinColorDiagonalNoiseMILC<FImpl>(g, nNoise), nSparse_(nSparse),
  coor_(g), coorTot_(g)
{
    if(nNoise%nSparse_==0)
    {
         nSrc_ec_ = nNoise/nSparse_;
    }
    else
    {
         nSrc_ec_ = (nNoise - nNoise%nSparse_)/nSparse_;
    }
}

template <typename FImpl>
int CheckerboardNoiseMILC<FImpl>::dilutionSize() const
{
    return this->getNsc()*this->size();
}

template <typename FImpl>
void CheckerboardNoiseMILC<FImpl>::setProp(const int i)
{
    auto nd    = this->getNd();
    auto noise = this->getNoise();
    auto nsc   = this->getNsc();
    unsigned int j;

    std::vector<LatticeComplex> eta(nsc,this->getGrid());

    j   = i/nSrc_ec_;

    coorTot_ = 0.;
    for(int d = 0; d < nd; ++d) 
    {
        LatticeCoordinate(coor_, d);
        coorTot_ = coorTot_ + coor_;
    }
    coor_ = j;
    coorTot_ = coorTot_ + coor_;

    for (int k=0;k<nsc;k++) {
        eta[k] = where(mod(coorTot_,nSparse_), 0.*noise[i*nsc+k], noise[i*nsc+k]);
        eta[k] *= sqrt(1./nSrc_ec_);
    }
    this->setPropagator(&(eta[0]));
}

/******************************************************************************
 *                SparseNoiseMILC template implementation                   *
 ******************************************************************************/
template <typename FImpl>
SparseNoiseMILC<FImpl>::
SparseNoiseMILC(GridCartesian *g, int nNoise, int nSparseL, int nSparseT)
: SpinColorDiagonalNoiseMILC<FImpl>(g, nNoise), nSparseL_(nSparseL), nSparseT_(nSparseT), coor_(g)
{}

template <typename FImpl>
int SparseNoiseMILC<FImpl>::dilutionSize() const
{
    auto nd  = this->getNd();
    return this->getNsc()*this->size()*pow(nSparseL_, nd-1)*nSparseT_;
}

template <typename FImpl>
void SparseNoiseMILC<FImpl>::setProp(const int i)
{
    auto nd    = this->getNd();
    auto noise = this->getNoise();
    auto nsc   = this->getNsc();

    std::vector<LatticeComplex> eta(nsc,this->getGrid());

    std::div_t divs = std::div(i, pow(nSparseL_, nd-1)*nSparseT_);
    std::div_t subdivs;

    for (int j=0;j<nsc;j++) {
        eta[j] = noise[divs.quot*nsc+j];
    }

    unsigned int sparseIndex = divs.rem;
    for(int d = 0; d < nd; ++d) {
        LatticeCoordinate(coor_, d);
        if (d < nd-1)
            subdivs = std::div(sparseIndex,nSparseL_);
        else
            subdivs = std::div(sparseIndex,nSparseT_);

        for (int j=0;j<nsc;j++) {
            eta[j] = where(coor_ == ((uint32_t) subdivs.rem), eta[j], 0.*eta[j]);
        }
        sparseIndex = subdivs.quot;
    }
    this->setPropagator(&(eta[0]));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_DilutedNoise_hpp_
