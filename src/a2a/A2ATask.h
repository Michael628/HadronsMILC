#pragma once

#include "A2AView.h"
#include "../spin/StagGamma.h"

#ifndef MF_SUM_ARRAY_MAX
#define MF_SUM_ARRAY_MAX 16
#endif

NAMESPACE_BEGIN(Grid);

#define A2A_TYPEDEFS \
typedef typename FImpl::SiteSpinor vobj; \
typedef typename FImpl::ComplexField ComplexField; \
typedef typename FImpl::FermionField FermionField; \
typedef typename ComplexField::vector_object cobj; \
typedef typename FImpl::ImplParams FImplParams; \
typedef CartesianStencil<vColourMatrix,vColourMatrix,FImplParams> GaugeStencil; \
typedef CartesianStencilView<vColourMatrix,vColourMatrix,FImplParams> GaugeStencilView; \
typedef typename vobj::scalar_type scalar_type; \
typedef typename vobj::vector_type vector_type; \
typedef iSinglet<scalar_type> Scalar_s; \
typedef iSinglet<vector_type> Scalar_v; \
typedef decltype(coalescedRead(Scalar_v())) calcScalar; \
typedef decltype(coalescedRead(vobj())) calcSpinor; \
typedef decltype(coalescedRead(vColourMatrix())) calcColourMatrix; \
  typedef typename FImpl::StencilImpl FermStencil; \
typedef typename FImpl::StencilView FermStencilView; \
typedef LatticeView<vColourMatrix> GaugeView; \
typedef LatticeView<cobj> ComplexView; \
typedef LatticeView<vobj> FermView;\
typedef typename A2ATaskBase<FImpl>::ContractType ContractType; \
typedef std::function<void(scalar_type*,cobj *)> SimdFunc; \
typedef std::function<void(cobj *, int, int)> VectorFunc;

#define COMMON_VARS \
int sizeL = this->_left_view->size(); \
int sizeR = this->_right_view->size(); \
\
int orthogDir = this->_orthog_dir; \
const int simdSize             = this->_grid->Nsimd(); \
const int reducedOrthogDimSize = this->_grid->_rdimensions[orthogDir]; \
\
const int nBlocks      = this->_grid->_slice_nblock[orthogDir]; \
\
const int localSpatialVolume     = this->_grid->_ostride[orthogDir]; \
\
FermView     *viewL_p    = this->_left_view->getView(); \
FermView     *viewR_p    = this->_right_view->getView(); \
int localOrthogDimSize   = this->_grid->_ldimensions[orthogDir]; \
\
int Nt     = this->_grid->GlobalDimensions()[orthogDir]; \
\
int pd = this->_grid->_processors[orthogDir]; \
int pc = this->_grid->_processor_coor[orthogDir]; \
\
int orthogSimdSize = this->_grid->_simd_layout[orthogDir]; \
Coordinate *icoor_p = this->_i_coor_container_device; \
Integer *ocoor_p = this->_o_coor_map_device; \
bool cbEven = this->_cb_left == Even; \
bool oddShifts = this->_odd_shifts;

template <typename FImpl>
class A2ATaskBase {
public:
    GRID_SERIALIZABLE_ENUM(ContractType, undef,
        Full,      0,
        RightHalf, 1,
        LeftHalf,  2,
        BothHalf,  3);

    A2A_TYPEDEFS;

protected:
    std::shared_ptr<A2AFieldView<vobj> > _left_view, _right_view;

    std::vector<Coordinate> _i_coor_container;
    Coordinate *_i_coor_container_device;

    std::vector<Integer> _o_coor_map;
    Integer *_o_coor_map_device;

    GridBase *_grid, *_full_grid;

    ContractType _contract_type;

    bool _odd_shifts;
    int _orthog_dir;
    const int _cb_left;

public:
    A2ATaskBase(GridBase *grid, int orthogDir, int cb = Even):
    _grid(grid), _full_grid(grid), _cb_left(cb),
    _orthog_dir(orthogDir), _odd_shifts(false), _contract_type(ContractType::undef) {

        _i_coor_container.resize(grid->Nsimd(), Coordinate(grid->_ndimension));
        for(int p = 0; p < grid->Nsimd(); p++) {
            grid->iCoorFromIindex(_i_coor_container[p],p);
        }

        size_t size = _i_coor_container.size()*sizeof(Coordinate);
        _i_coor_container_device = (Coordinate *)acceleratorAllocDevice(size);
        acceleratorCopyToDevice(_i_coor_container.data(),_i_coor_container_device,size);
    }

    virtual ~A2ATaskBase() {
        acceleratorFreeDevice(_i_coor_container_device);
        if (_o_coor_map.size() > 0) acceleratorFreeDevice(_o_coor_map_device);
    }

    virtual double getFlops() = 0;

    // Populates `_o_coor_map` property with full grid indices
    void generateCoorMap() {

        assert(_grid->CheckerBoarded(_orthog_dir) != 1);

        if (_o_coor_map.size() == 0) {
            _o_coor_map.resize(_grid->oSites(),0);
            _o_coor_map_device = (Integer *)acceleratorAllocDevice(_o_coor_map.size()*sizeof(Integer));

            int nBlocks = _grid->_slice_nblock[_orthog_dir];
            int vecsPerSlicePerBlock = _grid->_slice_block[_orthog_dir];
            int blockStride  = _grid->_slice_stride[_orthog_dir];
            int rtStride   = _grid->_ostride[_orthog_dir];
            int cb = _cb_left;

            thread_for(ss,_full_grid->oSites(),{
            int cbos;
            Coordinate coor;

            _full_grid->oCoorFromOindex(coor,ss);
            cbos=_grid->CheckerBoard(coor);
              
            if (cbos==cb) {
              int ssh=_grid->oIndex(coor);
              _o_coor_map[ssh]=ss;
            }
            });
            acceleratorCopyToDevice(_o_coor_map.data(),_o_coor_map_device,_o_coor_map.size()*sizeof(Integer));
        }
    }

    // Updates object to compute meson field with new `left` vectors
    // Updated properties:
    // `_contract_type` - Updated to reflect `left` checkerboard state
    // `_left_view`     - Allocates views from `left` pointer parameter
    // `_grid`          - Sets to checkerboarded grid if `_left_view` 
    //                    or `_right_view` are checkerboarded, otherwise set to full grid
    virtual void setLeft(const FermionField *left,int size) {

        bool checkerL = left[0].Grid()->_isCheckerBoarded;

        // Toggle LeftHalf bit
        if (_contract_type == ContractType::undef) {
            _contract_type = checkerL ? ContractType::LeftHalf : ContractType::Full;
        } else {
            if (checkerL)
                _contract_type = _contract_type | ContractType::LeftHalf;
            else
                _contract_type = _contract_type & ContractType::RightHalf;
        }

        switch(_contract_type) {
        case ContractType::LeftHalf:
        case ContractType::BothHalf:
        case ContractType::Full:
            _grid = left[0].Grid();
        default:
            break;
        }

        if (checkerL) generateCoorMap();

        _left_view = std::make_shared<A2AFieldView<vobj> >();
        _left_view->openViews(left,size);
    }

    // Updates object to compute meson field with new `right` vectors.
    // See corresponding `setLeft` method.
    virtual void setRight(const FermionField *right,int size) {

        bool checkerR = right[0].Grid()->_isCheckerBoarded;

        // Toggle RightHalf bit
        if (_contract_type == ContractType::undef) {
            _contract_type = checkerR ? ContractType::RightHalf : ContractType::Full;
        } else {
            if (checkerR)
                _contract_type = _contract_type | ContractType::RightHalf;
            else
                _contract_type = _contract_type & ContractType::LeftHalf;

        }

        switch(_contract_type) {
        case ContractType::RightHalf:
        case ContractType::BothHalf:
        case ContractType::Full:
            _grid = right[0].Grid();
        default:
            break;
        }

        if (checkerR) generateCoorMap();

        _right_view = std::make_shared<A2AFieldView<vobj> >();
        _right_view->openViews(right,size);
    }

    // Updates object to compute meson field with `other._left_view` vectors.
    // Should only be used for mixed cb/full calculation.
    virtual void setLeft(A2ATaskBase<FImpl> &other)  {
        assert(!(other.getType() & ContractType::LeftHalf));
        _left_view = other.getLeftView(); 

        _contract_type = other.getType();
    }

    virtual void setRight(A2ATaskBase<FImpl> &other) {
        assert(!(other.getType() & ContractType::RightHalf));
        _right_view = other.getRightView();

        _contract_type = other.getType();
    }

    std::shared_ptr<A2AFieldView<vobj> > getLeftView()  { return _left_view; }
    std::shared_ptr<A2AFieldView<vobj> > getRightView() { return _right_view; }
    ContractType getType() { return _contract_type; }
    

    virtual int getNgamma() = 0;
    virtual void vectorSumHalf(cobj *a,int b,int c)  = 0;
    virtual void vectorSumFull(cobj *a,int b,int c)  = 0;
    virtual void vectorSumMixed(cobj *a,int b,int c) = 0;

    virtual void execute(scalar_type *result_p) {

        COMMON_VARS;

        VectorFunc vectorSum;
        SimdFunc simdSum;

        int multFact;

        switch(this->_contract_type) {
        case ContractType::BothHalf:
            vectorSum = [this](cobj *a, int b,int c){this->vectorSumHalf(a,b,c);};
            simdSum = [this](scalar_type *a, cobj *b){this->simdSumHalf(a,b);};
            multFact = 4;
        break;
        case ContractType::LeftHalf:
        case ContractType::RightHalf:
            vectorSum = [this](cobj *a, int b,int c){this->vectorSumMixed(a,b,c);};
            simdSum = [this](scalar_type *a, cobj *b){this->simdSumMixed(a,b);};
            multFact = 2;
        break;
        case ContractType::Full:
            vectorSum = [this](cobj *a, int b, int c){this->vectorSumFull(a,b,c);};
            simdSum = [this](scalar_type *a, cobj *b){this->simdSumFull(a,b);};
            multFact = 1;
        break;
        }

        int nGamma = this->getNgamma();
        int gammaStride = sizeR*sizeL*reducedOrthogDimSize;
        int MFrvol = gammaStride * nGamma;

        cobj *shm_p = (cobj *)acceleratorAllocDevice(MFrvol*sizeof(cobj));

        // Loop over gammas in batches of MF_SUM_ARRAY_MAX
        for (int mu=0;mu<nGamma;mu+=MF_SUM_ARRAY_MAX) {

            int nGammaBlock = std::min(nGamma-mu,MF_SUM_ARRAY_MAX);

            vectorSum(shm_p + mu*gammaStride, mu, nGammaBlock);
        }

        for (int mu=0;mu<nGamma;mu++) {
            simdSum(result_p + mu*multFact*sizeR*sizeL*Nt, shm_p + mu*gammaStride);
        }

        acceleratorFreeDevice(shm_p);                
    }

    // Sums over SIMD vectorized results stored in `shm` parameter.
    // Case: Neither `_left_view` nor `_right_view` are checkerboarded.
    void simdSumFull(scalar_type *result, cobj *shm) {

        COMMON_VARS;

        auto shm_p = shm;
        auto result_p = result;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,1,{

            ExtractBuffer<Scalar_s> extracted(simdSize);
            scalar_type temp;

            int shmem_idx = reducedOrthogDimSize*(l_index+sizeL*r_index);

            for (int rt=0;rt<reducedOrthogDimSize;rt++) {

                extract(shm_p[shmem_idx+rt],extracted);

                for (int simdOffset=0;simdOffset<orthogSimdSize;simdOffset++) {

                    temp = scalar_type(0.0);
                    for(int idx=0;idx<simdSize;idx++){
                        if (icoor_p[idx][orthogDir] == simdOffset) {
                            temp +=  TensorRemove(extracted[idx]);
                        }
                        acceleratorSynchronise();
                    }

                    // Calculate local time from reduced time
                    int lt = rt+simdOffset*reducedOrthogDimSize;
                    int gt = lt + pc*localOrthogDimSize;

                    int ij_dx = r_index + sizeR*( l_index + sizeL*gt );

                    result_p[ij_dx] = temp;
                }
            }
        });
    }

    // Sums over SIMD vectorized results stored in `shm` parameter.
    // Case: Either `_left_view` or `_right_view` is checkerboarded, not both.
    void simdSumMixed(scalar_type *result, cobj *shm) {

        COMMON_VARS;

        auto shm_p = shm;
        auto result_p = result;

        bool fullLeft = this->_contract_type == ContractType::RightHalf;

        int multRight = fullLeft ? 2 : 1;
        int multLeft = fullLeft ? 1 : 2;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,1,{

            ExtractBuffer<Scalar_s> extracted(simdSize);
            scalar_type temp;

            int shmem_idx = reducedOrthogDimSize*(l_index+sizeL*r_index);

            for (int rt=0;rt<reducedOrthogDimSize;rt++) {

                extract(shm_p[shmem_idx+rt],extracted);

                for (int simdOffset=0;simdOffset<orthogSimdSize;simdOffset++) {

                    temp = scalar_type(0.0);
                    for(int idx=0;idx<simdSize;idx++){
                        if (icoor_p[idx][orthogDir] == simdOffset) {
                            temp +=  TensorRemove(extracted[idx]);
                        }
                        acceleratorSynchronise();
                    }

                    // Calculate local time from reduced time
                    int lt = rt+simdOffset*reducedOrthogDimSize;
                    int gt = lt + pc*localOrthogDimSize;

                    int ij_dx = multRight*(r_index + sizeR*multLeft*(l_index + sizeL*gt ) );

                    result_p[ij_dx] += temp;

                    if (fullLeft) {
                        if (cbEven && oddShifts) {
                            result_p[ij_dx+1]       -= temp;
                        } else if (cbEven) {
                            result_p[ij_dx+1]       += temp;
                        } else if (oddShifts) {
                            result_p[ij_dx+1]       += temp;
                        } else {
                            result_p[ij_dx+1]       -= temp;
                        }
                        acceleratorSynchronise();
                    } else {
                        if (cbEven && oddShifts) {
                            result_p[ij_dx+sizeR]   += temp;
                        } else if (cbEven) {
                            result_p[ij_dx+sizeR]   += temp;
                        } else if (oddShifts) {
                            result_p[ij_dx+sizeR]   -= temp;
                        } else {
                            result_p[ij_dx+sizeR]   -= temp;
                        }
                        acceleratorSynchronise();
                    }
                    acceleratorSynchronise();
                }
            }
        });
    }

    // Sums over SIMD vectorized results stored in `shm` parameter.
    // Case: Both `_left_view` and `_right_view` are checkerboarded.
    void simdSumHalf(scalar_type *result, cobj *shm) {

        COMMON_VARS;

        auto shm_p = shm;
        auto result_p = result;

        int sizeROut = 2*sizeR;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,1,{

            ExtractBuffer<Scalar_s> extracted(simdSize);
            scalar_type temp;

            int shmem_idx = reducedOrthogDimSize*(l_index+sizeL*r_index);

            for (int rt=0;rt<reducedOrthogDimSize;rt++) {

                extract(shm_p[shmem_idx+rt],extracted);

                for (int simdOffset=0;simdOffset<orthogSimdSize;simdOffset++) {

                    temp = scalar_type(0.0);
                    for(int idx=0;idx<simdSize;idx++){
                        if (icoor_p[idx][orthogDir] == simdOffset) {
                            temp +=  TensorRemove(extracted[idx]);
                        }
                        acceleratorSynchronise();
                    }

                    // Calculate local time and global time from reduced time
                    int lt = rt+simdOffset*reducedOrthogDimSize;
                    int gt = lt + pc*localOrthogDimSize;

                    int ij_dx = 2*( r_index + sizeR*2*( l_index + sizeL*gt ) );

                    result_p[ij_dx] += temp;

                    if (cbEven && oddShifts) {
                        result_p[ij_dx+1]          -= temp;
                        result_p[ij_dx+sizeROut+1] -= temp;
                        result_p[ij_dx+sizeROut]   += temp;
                    } else if (cbEven) {
                        result_p[ij_dx+1]          += temp;
                        result_p[ij_dx+sizeROut+1] += temp;
                        result_p[ij_dx+sizeROut]   += temp;

                    } else if (oddShifts) {
                        result_p[ij_dx+1]          += temp;
                        result_p[ij_dx+sizeROut+1] -= temp;
                        result_p[ij_dx+sizeROut]   -= temp;
                    } else {
                        result_p[ij_dx+1]          -= temp;
                        result_p[ij_dx+sizeROut+1] += temp;
                        result_p[ij_dx+sizeROut]   -= temp;
                    }
                    acceleratorSynchronise();
                }
            }
        });
    }
};

template <typename FImpl>
class A2ATaskLocal: public A2ATaskBase<FImpl> {
public: 
    A2A_TYPEDEFS;

protected:
    std::vector<StagGamma::SpinTastePair> _gammas;
    std::vector<ComplexField> _phase;
    std::shared_ptr<A2AFieldView<cobj> > _phase_view;

public:
    A2ATaskLocal(GridBase *grid, int orthogDir, A2ATaskLocal<FImpl> &other, const std::vector<StagGamma::SpinTastePair> &gammas = {}, int cb = Even):
    A2ATaskBase<FImpl>(grid,orthogDir,cb),_gammas(gammas)  {
        _phase_view = other.getPhaseView();
    }

    A2ATaskLocal(GridBase *grid, int orthogDir, const std::vector<StagGamma::SpinTastePair> &gammas, int cb = Even):
    A2ATaskBase<FImpl>(grid,orthogDir,cb),_gammas(gammas)  {

        int nGamma = _gammas.size();

        _phase.resize(nGamma,this->_full_grid);

        ComplexField temp(this->_full_grid);
        temp = 1.0;
        StagGamma spinTaste;

        _phase_view = std::make_shared<A2AFieldView<cobj> >();
        _phase_view->reserve(nGamma);

        for (int mu = 0; mu < nGamma; mu++) {

            spinTaste.setSpinTaste(_gammas[mu]);

            spinTaste.applyPhase(_phase[mu],temp); // store spin-taste phase
        }
        _phase_view->openViews(_phase.data(),nGamma);
    }

    A2ATaskLocal(GridBase *grid, int orthogDir, const std::vector<ComplexField> &gammas, int cb = Even):
    A2ATaskBase<FImpl>(grid,orthogDir,cb)  {

        int nGamma = gammas.size();

        _phase_view = std::make_shared<A2AFieldView<cobj> >();
        _phase_view->reserve(nGamma);
        _phase_view->openViews(gammas.data(),nGamma);
    }

    virtual ~A2ATaskLocal() { 
        if (_phase_view) _phase_view->closeViews();
    }

    std::shared_ptr<A2AFieldView<cobj> > getPhaseView() { return _phase_view; }

    virtual int getNgamma() { return _phase_view->size(); }

    virtual double getFlops() {
        // One complex multiply takes 6 floating point ops (4 mult, 2 add) 
        // --> complex inner product is 3 complex mult, 2 complex add = 3*6 + 2*2 = 22 double precision floating ops

        // For each vector and at each lattice site:
        //  - one inner product
        //  - For each gamma and momentum
        //    - multiply by gamma phase
        //    - sum
        return (22.0+(6.0+2.0)*(this->getNgamma()));
    }

    virtual void vectorSumHalf(cobj *shm_p, int mu_offset, int N) {

        COMMON_VARS;

        ComplexView  *viewG_p    = this->_phase_view->getView() + mu_offset;

        assert(orthogDir == Tdir); // This kernel assumes lattice is coalesced over time slices
        assert(nBlocks == 1);

        int gammaStride = sizeR*sizeL*reducedOrthogDimSize;
        int nGamma = N;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,simdSize,{

            int ss, shmem_base = reducedOrthogDimSize*(l_index+sizeL*r_index);
            calcScalar gamma_phase, temp_site, sum[MF_SUM_ARRAY_MAX];

            for (int rt=0;rt < reducedOrthogDimSize; rt++) {

                for (int mu=0;mu<nGamma;mu++) {
                    sum[mu] = Zero();
                }

                for (int so=0;so < localSpatialVolume; so++) {

                    ss = rt*localSpatialVolume+so;

                    temp_site = innerProduct(coalescedRead(viewL_p[l_index][ss]),coalescedRead(viewR_p[r_index][ss]));

                    for (int mu = 0; mu < nGamma; mu++) {
                        gamma_phase = coalescedRead(viewG_p[mu][ocoor_p[ss]]);
                        sum[mu] += gamma_phase*temp_site;
                    }
                }

                for (int mu=0;mu<nGamma;mu++) {
                    int shmem_idx = rt+shmem_base+mu*gammaStride;
                    coalescedWrite(shm_p[shmem_idx],sum[mu]);
                }
            }
        });
    }

    virtual void vectorSumFull(cobj *shm_p, int mu_offset, int N) {

        COMMON_VARS;

        ComplexView  *viewG_p    = this->_phase_view->getView() + mu_offset;

        assert(orthogDir == Tdir); // This kernel assumes lattice is coalesced over time slices
        assert(nBlocks == 1);

        int gammaStride = sizeR*sizeL*reducedOrthogDimSize;
        int nGamma = N;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,simdSize,{

            int ss, shmem_base = reducedOrthogDimSize*(l_index+sizeL*r_index);
            calcScalar gamma_phase, temp_site, sum[MF_SUM_ARRAY_MAX];

            for (int rt=0;rt < reducedOrthogDimSize; rt++) {

                for (int mu=0;mu<nGamma;mu++) {
                    sum[mu] = Zero();
                }

                for (int so=0;so < localSpatialVolume; so++) {

                    ss = rt*localSpatialVolume+so;

                    temp_site = innerProduct(coalescedRead(viewL_p[l_index][ss]),coalescedRead(viewR_p[r_index][ss]));

                    for (int mu = 0; mu < nGamma; mu++) {
                        gamma_phase = coalescedRead(viewG_p[mu][ss]);
                        sum[mu] += gamma_phase*temp_site;
                    }
                }

                for (int mu=0;mu<nGamma;mu++) {
                    int shmem_idx = rt+shmem_base+mu*gammaStride;
                    coalescedWrite(shm_p[shmem_idx],sum[mu]);
                }
            }
        });
    }

    virtual void vectorSumMixed(cobj *shm_p, int mu_offset, int N) {

        COMMON_VARS;

        ComplexView  *viewG_p    = this->_phase_view->getView() + mu_offset;

        assert(orthogDir == Tdir); // This kernel assumes lattice is coalesced over time slices
        assert(nBlocks == 1);

        int gammaStride = sizeR*sizeL*reducedOrthogDimSize;
        int nGamma = N;

        bool checkerL = this->_contract_type == ContractType::LeftHalf;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,simdSize,{

            int ss, shmem_base = reducedOrthogDimSize*(l_index+sizeL*r_index);
            calcScalar gamma_phase, temp_site, sum[MF_SUM_ARRAY_MAX];

            for (int rt=0;rt < reducedOrthogDimSize; rt++) {

                for (int mu=0;mu<nGamma;mu++) {
                    sum[mu] = Zero();
                }

                for (int so=0;so < localSpatialVolume; so++) {

                    ss = rt*localSpatialVolume+so;

                    if (checkerL) {
                        temp_site = innerProduct(coalescedRead(viewL_p[l_index][ss]),coalescedRead(viewR_p[r_index][ocoor_p[ss]]));
                    } else {
                        temp_site = innerProduct(coalescedRead(viewL_p[l_index][ocoor_p[ss]]),coalescedRead(viewR_p[r_index][ss]));
                    }
                    acceleratorSynchronise();

                    for (int mu = 0; mu < nGamma; mu++) {
                        gamma_phase = coalescedRead(viewG_p[mu][ocoor_p[ss]]);
                        sum[mu] += gamma_phase*temp_site;
                    }
                }

                for (int mu=0;mu<nGamma;mu++) {
                    int shmem_idx = rt+shmem_base+mu*gammaStride;
                    coalescedWrite(shm_p[shmem_idx],sum[mu]);
                }
            }
        });
    }
};

template <typename FImpl>
class A2ATaskOnelink: public A2ATaskBase<FImpl> {
public: 
    A2A_TYPEDEFS;

protected:
    const std::vector<StagGamma::SpinTastePair> &_gammas = {};
    std::vector<int> _shift_dirs, _shift_displacements;
    LatticeGaugeField *_U;

    std::vector<LatticeColourMatrix>                     _link,        _link_shifted;
    std::shared_ptr<A2AFieldView<vColourMatrix> >        _link_view,   _link_shifted_view;
    std::shared_ptr<A2AStencilView<vobj,FImplParams> >   _right_stencil_view;

public:

    A2ATaskOnelink(GridBase *grid, int orthogDir, const std::vector<StagGamma::SpinTastePair> &gammas, int cb = Even)
    : A2ATaskBase<FImpl>(grid,orthogDir,cb),_gammas(gammas)  {

        this->_odd_shifts = true;

        StagGamma spinTaste;

        for (int i = 0; i < _gammas.size(); i++) {

          spinTaste.setSpinTaste(_gammas[i]);

            int shift = (spinTaste._spin ^ spinTaste._taste);
            // Assume 1-link for now -- break loop when you find a shift direction
            for (int j = 0; j < StagGamma::gmu.size(); j++) {
              if (StagGamma::gmu[j] & shift) {
                _shift_dirs.push_back(j);
                _shift_dirs.push_back(j);
                _shift_displacements.push_back(1);
                _shift_displacements.push_back(-1);
                break;
              }
            }
        }
    }

    A2ATaskOnelink(GridBase *grid, int orthogDir, A2ATaskOnelink<FImpl> &other, const std::vector<StagGamma::SpinTastePair> &gammas, LatticeGaugeField* U, int cb = Even)
    : A2ATaskOnelink<FImpl>(grid,orthogDir,gammas,cb) {
        _U = U;
        _link_view = other.getLinkView();
        _link_shifted_view = other.getLinkView(true);

    }

    A2ATaskOnelink(GridBase *grid, int orthogDir, const std::vector<StagGamma::SpinTastePair> &gammas, LatticeGaugeField* U, int cb = Even)
    : A2ATaskOnelink<FImpl>(grid,orthogDir,gammas,cb)  {

        _U = U;

        double t0 = usecond();

        int nGamma = _gammas.size();

        _link.resize(nGamma,_U->Grid());
        _link_shifted.resize(nGamma,_U->Grid());

        StagGamma spinTaste;

        for (int i = 0; i < nGamma; i++) {

            spinTaste.setSpinTaste(_gammas[i]);

            _link[i] = PeekIndex<LorentzIndex>(*_U,_shift_dirs[2*i]); // Store full lattice links in shift direction

            _link_shifted[i] = Cshift(_link[i],_shift_dirs[2*i],-1);

            spinTaste.applyPhase(_link[i],_link[i]); // store spin-taste phase
            spinTaste.applyPhase(_link_shifted[i],_link_shifted[i]); // store spin-taste phase

        
        }
        _link_view = std::make_shared<A2AFieldView<vColourMatrix> >();
        _link_shifted_view = std::make_shared<A2AFieldView<vColourMatrix> >();

        _link_view->openViews(_link.data(),nGamma);
        _link_shifted_view->openViews(_link_shifted.data(),nGamma);

        double t1 = usecond();
        std::cout << GridLogPerformance << " MesonField onelink timings: build link fields + comms: " << (t1-t0)/1000 << " ms" << std::endl;   
    }

    // Updates object to compute meson field with new `right` vectors.
    // See corresponding `setLeft` method.
    virtual void setRight(const FermionField *right,int size) {

        bool checkerR = right[0].Grid()->_isCheckerBoarded;

        // Toggle RightHalf bit
        if (this->_contract_type == ContractType::undef) {
            this->_contract_type = checkerR ? ContractType::RightHalf : ContractType::Full;
        } else {
            if (checkerR)
                this->_contract_type = this->_contract_type | ContractType::RightHalf;
            else
                this->_contract_type = this->_contract_type & ContractType::LeftHalf;

        }

        int cb = (this->_cb_left == Odd) ? Even : Odd;

        switch(this->_contract_type) {
        case ContractType::Full:
            cb = Even;
        case ContractType::BothHalf:
        case ContractType::RightHalf:
            this->_grid = right[0].Grid();
            break;
        case ContractType::LeftHalf:
            cb = Even;
        default:
            break;
        }

        if (checkerR) this->generateCoorMap();

        auto ptr = std::make_unique<FermStencil>(this->_grid, _shift_dirs.size(), cb,
                                                                _shift_dirs, _shift_displacements);

        _right_stencil_view = std::make_shared<A2AStencilView<vobj,FImplParams> >();

        _right_stencil_view->addStencil(ptr);
        _right_stencil_view->openViews(right,size);

        this->_right_view = std::make_shared<A2AFieldView<vobj> >();
        this->_right_view->openViews(right,size);
    }

    virtual void setRight(A2ATaskBase<FImpl> &other) {
        assert(!(other.getType() & ContractType::RightHalf));
        this->_right_view = other.getRightView();
        _right_stencil_view = dynamic_cast<A2ATaskOnelink<FImpl> &>(other).getRightStencilView();

        this->_contract_type = other.getType();
    }

    virtual ~A2ATaskOnelink() { 
        if (_right_stencil_view) _right_stencil_view->closeViews(); 
        if (_link_view) _link_view->closeViews(); 
        if (_link_shifted_view) _link_shifted_view->closeViews(); 
    }

    std::shared_ptr<A2AStencilView<vobj,FImplParams> > getRightStencilView() { return _right_stencil_view; }

    std::shared_ptr<A2AFieldView<vColourMatrix> >   getLinkView(bool shifted = false) 
    { 
        if (shifted) 
            return _link_shifted_view;
        else 
            return _link_view; 
    }

    virtual int getNgamma() { return _link_view->size(); }

    virtual double getFlops() {
      // matrix*vector = 3 inner products
      // current code: innerProduct(left,link_ahead*shift_ahead+adj(link_behind)*shift_behind)
      //  = matrix*vector + matrix* vector --> inner product = 7 inner products and 1 complex sum

      // For each vector, each gamma, and at each lattice site:
      //  - one inner product
      //  - two su(3) matrix*vector ops
      //  - one complex sum
      return ((7*22.0+2.0)*this->getNgamma());
    }

    virtual void vectorSumHalf(cobj *shm_p, int mu_offset, int N) {

        COMMON_VARS;

        assert(orthogDir == Tdir); // This kernel assumes lattice is coalesced over time slices
        assert(nBlocks == 1);

        // Pointers for accelerator indexing
        GaugeView *viewGL_p   = this->_link_view->getView() + mu_offset;
        GaugeView *viewGR_p   = this->_link_shifted_view->getView() + mu_offset;

        FermStencilView  *stencilR_p = this->_right_stencil_view->getView(); // ket stencil for shifted links
        vobj             *bufRight_p = this->_right_stencil_view->getBuffer(); // buffer for shifted kets in halo region

        int haloBuffRightSize = this->_right_stencil_view->getOffset(1);

        int gammaStride = sizeR*sizeL*reducedOrthogDimSize;
        int nGamma = N;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,simdSize,{
            calcColourMatrix link_ahead, link_behind;
            calcSpinor left, shift_ahead, shift_behind;
            calcScalar sum[MF_SUM_ARRAY_MAX];

            StencilEntry *SE;
            int ptype, ss;
            int shmem_base = reducedOrthogDimSize*(l_index+sizeL*r_index);

            for (int rt=0;rt < reducedOrthogDimSize; rt++) {

                for (int mu=0;mu<nGamma;mu++) {
                    sum[mu] = Zero();
                }

                for (int so=0;so < localSpatialVolume; so++) {
                    ss = rt*localSpatialVolume+so;

                    left   = coalescedRead(viewL_p[l_index][ss]);

                    for (int mu = 0; mu < nGamma; mu++) {
                        int rightShiftOffset = 2*(mu+mu_offset);
                        link_ahead = coalescedRead(viewGL_p[mu][ocoor_p[ss]]);

                        SE=stencilR_p->GetEntry(ptype,rightShiftOffset,ss);
                        if(SE->_is_local) {
                            shift_ahead = coalescedReadPermute(viewR_p[r_index][SE->_offset],ptype,SE->_permute);
                        } else {
                            shift_ahead = coalescedRead(bufRight_p[r_index*haloBuffRightSize+SE->_offset]);
                        }
                        acceleratorSynchronise();

                        SE=stencilR_p->GetEntry(ptype,rightShiftOffset+1,ss);
                        if(SE->_is_local) {
                            shift_behind = coalescedReadPermute(viewR_p[r_index][SE->_offset],ptype,SE->_permute);
                        } else {
                            shift_behind = coalescedRead(bufRight_p[r_index*haloBuffRightSize+SE->_offset]);
                        }
                        acceleratorSynchronise();

                        link_behind = adj(coalescedRead(viewGR_p[mu][ocoor_p[ss]]));

                        sum[mu] += innerProduct(left,link_ahead*shift_ahead+link_behind*shift_behind);
                    }
                    for (int mu=0;mu<nGamma;mu++) {
                        int shmem_idx = rt+shmem_base+mu*gammaStride;
                        coalescedWrite(shm_p[shmem_idx],sum[mu]);
                    }
                }
            }
        });
    }

    virtual void vectorSumFull(cobj *shm_p, int mu_offset, int N) { 

        COMMON_VARS;

        assert(orthogDir == Tdir); // This kernel assumes lattice is coalesced over time slices
        assert(nBlocks == 1);

        // Pointers for accelerator indexing
        GaugeView *viewGL_p   = this->_link_view->getView() + mu_offset;
        GaugeView *viewGR_p   = this->_link_shifted_view->getView() + mu_offset;

        FermStencilView  *stencilR_p = this->_right_stencil_view->getView(); // ket stencil for shifted links
        vobj             *bufRight_p = this->_right_stencil_view->getBuffer(); // buffer for shifted kets in halo region

        int haloBuffRightSize = this->_right_stencil_view->getOffset(1);

        int gammaStride = sizeR*sizeL*reducedOrthogDimSize;
        int nGamma = N;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,simdSize,{
            calcColourMatrix link_ahead, link_behind;
            calcSpinor left, shift_ahead, shift_behind;
            calcScalar sum[MF_SUM_ARRAY_MAX];

            StencilEntry *SE;
            int ptype, ss;
            int shmem_base = reducedOrthogDimSize*(l_index+sizeL*r_index);

            for (int rt=0;rt < reducedOrthogDimSize; rt++) {

                for (int mu=0;mu<nGamma;mu++) {
                    sum[mu] = Zero();
                }

                for (int so=0;so < localSpatialVolume; so++) {
                    ss = rt*localSpatialVolume+so;

                    left   = coalescedRead(viewL_p[l_index][ss]);

                    for (int mu = 0; mu < nGamma; mu++) {
                        int rightShiftOffset = 2*(mu+mu_offset);
                        link_ahead = coalescedRead(viewGL_p[mu][ss]);

                        SE=stencilR_p->GetEntry(ptype,rightShiftOffset,ss);
                        if(SE->_is_local) {
                            shift_ahead = coalescedReadPermute(viewR_p[r_index][SE->_offset],ptype,SE->_permute);
                        } else {
                            shift_ahead = coalescedRead(bufRight_p[r_index*haloBuffRightSize+SE->_offset]);
                        }
                        acceleratorSynchronise();

                        SE=stencilR_p->GetEntry(ptype,rightShiftOffset+1,ss);
                        if(SE->_is_local) {
                            shift_behind = coalescedReadPermute(viewR_p[r_index][SE->_offset],ptype,SE->_permute);
                        } else {
                            shift_behind = coalescedRead(bufRight_p[r_index*haloBuffRightSize+SE->_offset]);
                        }
                        acceleratorSynchronise();

                        link_behind = adj(coalescedRead(viewGR_p[mu][ss]));

                        sum[mu] += innerProduct(left,link_ahead*shift_ahead+link_behind*shift_behind);
                    }
                    for (int mu=0;mu<nGamma;mu++) {
                        int shmem_idx = rt+shmem_base+mu*gammaStride;
                        coalescedWrite(shm_p[shmem_idx],sum[mu]);
                    }
                }
            }
        });
    }

    virtual void vectorSumMixed(cobj *shm_p, int mu_offset, int N) { 

        COMMON_VARS;

        assert(orthogDir == Tdir); // This kernel assumes lattice is coalesced over time slices
        assert(nBlocks == 1);

        // Pointers for accelerator indexing
        GaugeView *viewGL_p   = this->_link_view->getView() + mu_offset;
        GaugeView *viewGR_p   = this->_link_shifted_view->getView() + mu_offset;

        FermStencilView  *stencilR_p = this->_right_stencil_view->getView(); // ket stencil for shifted links
        vobj             *bufRight_p = this->_right_stencil_view->getBuffer(); // buffer for shifted kets in halo region

        int haloBuffRightSize = this->_right_stencil_view->getOffset(1);

        int gammaStride = sizeR*sizeL*reducedOrthogDimSize;
        int nGamma = N;

        bool checkerL = this->_contract_type == ContractType::LeftHalf;

        accelerator_for2d(l_index,sizeL,r_index,sizeR,simdSize,{
            calcColourMatrix link_ahead, link_behind;
            calcSpinor left, shift_ahead, shift_behind;
            calcScalar sum[MF_SUM_ARRAY_MAX];

            StencilEntry *SE;
            int ptype, ss,ss_left,ss_right,ss_link;
            int shmem_base = reducedOrthogDimSize*(l_index+sizeL*r_index);

            for (int rt=0;rt < reducedOrthogDimSize; rt++) {

                for (int mu=0;mu<nGamma;mu++) {
                    sum[mu] = Zero();
                }

                for (int so=0;so < localSpatialVolume; so++) {
                    ss = rt*localSpatialVolume+so;

                    ss_link = ocoor_p[ss];
                    if (checkerL) {
                        ss_left = ss;
                        ss_right = ss_link;
                    } else {
                        ss_left = ss_link;
                        ss_right = ss;
                    }
                    acceleratorSynchronise();

                    left   = coalescedRead(viewL_p[l_index][ss_left]);


                    for (int mu = 0; mu < nGamma; mu++) {
                        int rightShiftOffset = 2*(mu+mu_offset);
                        link_ahead = coalescedRead(viewGL_p[mu][ss_link]);

                        SE=stencilR_p->GetEntry(ptype,rightShiftOffset,ss_right);
                        if(SE->_is_local) {
                            shift_ahead = coalescedReadPermute(viewR_p[r_index][SE->_offset],ptype,SE->_permute);
                        } else {
                            shift_ahead = coalescedRead(bufRight_p[r_index*haloBuffRightSize+SE->_offset]);
                        }
                        acceleratorSynchronise();

                        SE=stencilR_p->GetEntry(ptype,rightShiftOffset+1,ss_right);
                        if(SE->_is_local) {
                            shift_behind = coalescedReadPermute(viewR_p[r_index][SE->_offset],ptype,SE->_permute);
                        } else {
                            shift_behind = coalescedRead(bufRight_p[r_index*haloBuffRightSize+SE->_offset]);
                        }
                        acceleratorSynchronise();

                        link_behind = adj(coalescedRead(viewGR_p[mu][ss_link]));

                        sum[mu] += innerProduct(left,link_ahead*shift_ahead+link_behind*shift_behind);
                    }
                    for (int mu=0;mu<nGamma;mu++) {
                        int shmem_idx = rt+shmem_base+mu*gammaStride;
                        coalescedWrite(shm_p[shmem_idx],sum[mu]);
                    }
                }
            }
        });
    }
};
#undef A2A_TYPEDEFS
#undef COMMON_VARS

NAMESPACE_END(Grid);
