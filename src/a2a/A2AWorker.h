#pragma once

#include <Grid/Grid_Eigen_Tensor.h>
#include "A2ATask.h"
#include "../spin/StagGamma.h"

NAMESPACE_BEGIN(Grid);

template <typename FImpl>
class A2AWorkerBase 
{
public:
    typedef typename FImpl::ComplexField ComplexField;
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::scalar_type scalar_type;

public:
    GridBase *_grid, *_cb_grid;

    double _flops,_t_kernel,_t_gsum;

    A2ATaskBase<FImpl> *_task_e, *_task_o;

    const FermionField *_l_addr, *_r_addr;

    scalar_type *_cache_device;
    size_t _cache_bytes = 0;

    bool _odd_shifts {false};
public:
    A2AWorkerBase() = delete;
    A2AWorkerBase(GridBase *grid) : _grid(grid),_l_addr(nullptr),_r_addr(nullptr) {}

    virtual ~A2AWorkerBase() {
        if (_cache_bytes != 0) {
            acceleratorFreeDevice(_cache_device);
        }
        delete _task_e;
        delete _task_o;
    }
public:
    template <typename TensorType> // output: rank 5 tensor, e.g. Eigen::Tensor<ComplexD, 5>
    void StagMesonField(TensorType &mat,
        const FermionField *lhs_wi_E, const FermionField *lhs_wi_O,
        const FermionField *rhs_vj_E, const FermionField *rhs_vj_O);

    void setFlops(double flops) {_flops = flops; }
    double getFlops() {return _flops; }
};

template <typename FImpl>
class A2AWorkerLocal: public A2AWorkerBase<FImpl>
{
public:
    typedef typename FImpl::ComplexField ComplexField;
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::scalar_type scalar_type;
public:
    A2AWorkerLocal() = delete;
    A2AWorkerLocal(GridBase *grid, const std::vector<ComplexField> &mom, const std::vector<StagGamma::SpinTastePair> &gammas, int orthogDir)
    : A2AWorkerBase<FImpl>(grid) {
        this->_odd_shifts = false;
        if(mom.size()) {
            assert(0);
        } else {
            this->_task_e = new A2ATaskLocal<FImpl>(grid, orthogDir, gammas, Even);
            this->_task_o = new A2ATaskLocal<FImpl>(grid, orthogDir, dynamic_cast<A2ATaskLocal<FImpl> &>(*this->_task_e), gammas, Odd);
        }
    }
    A2AWorkerLocal(GridBase *grid, const std::vector<ComplexField> &mom, const std::vector<ComplexField> &Amu, int orthogDir)
    : A2AWorkerBase<FImpl>(grid) {
        this->_odd_shifts = false;
        if(mom.size()) {
            assert(0);
        } else {
            this->_task_e = new A2ATaskLocal<FImpl>(grid,orthogDir,Amu,Even);
            this->_task_o = new A2ATaskLocal<FImpl>(grid,orthogDir,dynamic_cast<A2ATaskLocal<FImpl> &>(*this->_task_e),{},Odd);
        }
    }
};

template <typename FImpl>
class A2AWorkerOnelink: public A2AWorkerBase<FImpl>
{
public:
    typedef typename FImpl::ComplexField ComplexField;
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::scalar_type scalar_type;
public:
    A2AWorkerOnelink() = delete;
    A2AWorkerOnelink(GridBase *grid, const std::vector<ComplexField> &mom, const std::vector<StagGamma::SpinTastePair> &gammas, LatticeGaugeField* U, int orthogDir)
    : A2AWorkerBase<FImpl>(grid) {
        this->_odd_shifts = true;
        if(mom.size()) {
            assert(0);
        } else {
            this->_task_e = new A2ATaskOnelink<FImpl>(grid, orthogDir, gammas,U,Even);
            this->_task_o = new A2ATaskOnelink<FImpl>(grid, orthogDir, dynamic_cast<A2ATaskOnelink<FImpl> &>(*this->_task_e), gammas, U, Odd);
        }
    }
    A2AWorkerOnelink(GridBase *grid, const std::vector<ComplexField> &mom, const std::vector<ComplexField> &Amu, LatticeGaugeField* U, int orthogDir)
    : A2AWorkerBase<FImpl>(grid) {
        // A2A onelink EM not implemented yet
        assert(0);
        this->_odd_shifts = true;
        // if(mom.size()) {
            // assert(0);
        // } else {
            // this->_task_e = new A2ATaskOnelink<FImpl>(grid,orthogDir,Amu,Even,U);
            // this->_task_o = new A2ATaskOnelink<FImpl>(grid,orthogDir,dynamic_cast<A2ATaskOnelink<FImpl> &>(*this->_task_e),{},Odd,U);
        // }
    }
};

template <class FImpl>
template <typename TensorType>
void A2AWorkerBase<FImpl>::StagMesonField(TensorType &mat,
    const FermionField *lhs_wi_E,
    const FermionField *lhs_wi_O,
    const FermionField *rhs_vj_E,
    const FermionField *rhs_vj_O)
{
    if (_cache_bytes < mat.size()*sizeof(scalar_type)) {
        if (_cache_bytes != 0) {
            acceleratorFreeDevice(_cache_device);
        }
        _cache_bytes = mat.size()*sizeof(scalar_type);
        std::cout << GridLogPerformance << "cache bytes: " << _cache_bytes << std::endl;
        _cache_device = (scalar_type*)acceleratorAllocDevice(_cache_bytes);
    }

    int sizeL = mat.dimension(3);
    int sizeR = mat.dimension(4);
    bool checkerL = lhs_wi_E[0].Grid()->_isCheckerBoarded;
    bool checkerR = rhs_vj_E[0].Grid()->_isCheckerBoarded;

    if (checkerL) sizeL /= 2;
    if (checkerR) sizeR /= 2;

    // scalar_type *matDevice = _cache_device;
    scalar_type *matDevice = mat.data();

    if (_l_addr != lhs_wi_E) {
        _l_addr = lhs_wi_E;

        _task_e->setLeft(lhs_wi_E,sizeL);
        if (checkerL) _task_o->setLeft(lhs_wi_O,sizeL);
        else if (checkerR) _task_o->setLeft(*_task_e);
    }

    if (_r_addr != rhs_vj_E) {
        _r_addr = rhs_vj_E;

        if (checkerR) {
            if (_odd_shifts) {
                _task_e->setRight(rhs_vj_O,sizeR);
                _task_o->setRight(rhs_vj_E,sizeR);
            } else {
                _task_e->setRight(rhs_vj_E,sizeR);
                _task_o->setRight(rhs_vj_O,sizeR);
            }
        } else {
            _task_e->setRight(rhs_vj_E,sizeR);
            if (checkerL) {
                _task_o->setRight(*_task_e);
            }
        }
    }

    _t_kernel = -usecond();
    if (!(checkerL || checkerR)) {
        _task_e->execute(matDevice);
    } else {
        _task_e->execute(matDevice);
        _task_o->execute(matDevice);
    }

    setFlops(_task_e->getFlops());
    _t_kernel += usecond();

    int resultStride = mat.dimension(2)*mat.dimension(3)*mat.dimension(4);
    int nGamma = mat.dimension(0)*mat.dimension(1);

    // auto matHost = mat.data();
    // size_t matBytes = mat.size()*sizeof(scalar_type);
    // acceleratorCopyFromDevice(matDevice,matHost,matBytes);

    _t_gsum = -usecond();
    int comm_batch = std::pow(2,20);
    int offset = 0;
    int buff = comm_batch;

    // while (offset < resultStride*nGamma) {
    //     buff = std::min(comm_batch,resultStride*nGamma - offset);
    //     this->_grid->GlobalSumVector(matDevice+offset,buff);
    //     offset += buff;
    // }


    // for (int i=0; i<nGamma;i++) {
        // this->_grid->GlobalSumVector(matDevice+i*resultStride,resultStride);
    // }

    //this->_grid->GlobalSumVector(matDevice,nGamma*resultStride);
    // this->_grid->GlobalSumVector(matHost,nGamma*resultStride);
    this->_grid->GlobalSumVector(matDevice,nGamma*resultStride);
    _t_gsum += usecond();

}

NAMESPACE_END(Grid);

