#pragma once
#include <Grid/Grid_Eigen_Tensor.h>
#include <spin/StagGamma.h>

NAMESPACE_BEGIN(Grid);

template <typename FImpl> class A2AWorkerMILC {
public:
  typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::FermionField FermionField;
  typedef typename FImpl::PropagatorField PropagatorField;

  typedef typename FImpl::SiteSpinor vobj;

  typedef typename FImpl::ImplParams FImplParams;

  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;

  typedef decltype(coalescedRead(Scalar_v())) calcScalar;
  typedef decltype(coalescedRead(vobj())) calcSpinor;
  typedef decltype(coalescedRead(vColourMatrix())) calcColourMatrix;

  typedef LatticeView<vobj> FermView;
  typedef typename FImpl::StencilImpl FermStencil;
  typedef typename FImpl::StencilView FermStencilView;

  typedef LatticeView<vColourMatrix> GaugeView;
  typedef CartesianStencil<vColourMatrix, vColourMatrix, FImplParams>
      GaugeStencil;
  typedef CartesianStencilView<vColourMatrix, vColourMatrix, FImplParams>
      GaugeStencilView;

  typedef typename ComplexField::vector_object cobj;
  typedef LatticeView<cobj> ComplexView;

public:
  template <typename Vtype, typename obj> class A2AViewBase {
  public:
    Vector<Vtype> _view;

  public:
    Vtype &operator[](size_t i) { return _view[i]; }

    int size() { return _view.size(); }
    void reserve(int size) { _view.reserve(size); }

    virtual void closeViews() {
      for (int p = 0; p < _view.size(); p++)
        _view[p].ViewClose();
      _view.resize(0);
    }
  };

  template <typename obj>
  class A2AFieldView : public A2AViewBase<LatticeView<obj>, obj> {
  public:
    void addView(const Lattice<obj> &field) {
      this->_view.push_back(field.View(AcceleratorRead));
    }
  };

  template <typename obj>
  class A2AStencilView
      : public A2AViewBase<CartesianStencilView<obj, obj, FImplParams>, obj> {
  protected:
    std::vector<std::unique_ptr<CartesianStencil<obj, obj, FImplParams>>>
        _stencils;
    Vector<obj> _buffer;
    Vector<Integer> _offset;

  public:
    CartesianStencil<obj, obj, FImplParams> &getStencil(int i) {
      return *_stencils[i];
    }
    obj &getBuffer(int i) { return _buffer[i]; }
    Integer &getOffset(int i) { return _offset[i]; }

    void addStencil(
        std::unique_ptr<CartesianStencil<obj, obj, FImplParams>> &stencil) {
      _stencils.push_back(std::move(stencil));
    }

    void openViews() {
      this->reserve(_stencils.size());
      for (auto &stencil : _stencils) {
        this->_view.push_back(stencil->View(AcceleratorRead));
      }
    }

    void append(const Lattice<obj> &field) {
      GridBase *grid = field.Grid();
      SimpleCompressor<obj> compressor;

      auto &stencil = _stencils.back();

      int comm_buf_size;
      obj *buf_p;

      stencil->HaloExchange(field, compressor);

      comm_buf_size = stencil->_unified_buffer_size;
      _offset.push_back(_buffer.size());

      _buffer.resize(_buffer.size() + comm_buf_size);
      buf_p = &_buffer[_offset.back()];

      if (comm_buf_size > 0) {
        obj *comm_buf_p = stencil->CommBuf();
        accelerator_for(i, comm_buf_size, 1, { buf_p[i] = comm_buf_p[i]; });
      }
    }
    virtual void closeViews() {
      for (int p = 0; p < this->_view.size(); p++)
        this->_view[p].ViewClose();
      this->_view.resize(0);
      _buffer.resize(0);
      _offset.resize(0);
      _stencils.resize(0);
    }
  };

  class A2ATaskBase {
  protected:
    A2AFieldView<vobj> &_left, &_right;
    Vector<Integer> &_gamma_indices;
    GridBase *_grid;

    Vector<Coordinate> _i_coor_container;

  public:
    A2ATaskBase(A2AFieldView<vobj> &left, A2AFieldView<vobj> &right,
                Vector<Integer> &gammaIndices, GridBase *grid)
        : _left(left), _right(right), _gamma_indices(gammaIndices),
          _grid(grid) {

      _i_coor_container.resize(grid->Nsimd(), Coordinate(grid->_ndimension));
      for (int p = 0; p < grid->Nsimd(); p++) {
        grid->iCoorFromIindex(_i_coor_container[p], p);
      }
    }
  };

  class A2ATaskHalfHalf : public A2ATaskBase {
  protected:
    bool _odd_shifts;
    int _cb;

  public:
    A2ATaskHalfHalf(A2AFieldView<vobj> &left, A2AFieldView<vobj> &right,
                    Vector<Integer> &gammaIndices, GridBase *grid,
                    int cb = Even, bool oddShifts = false)
        : A2ATaskBase(left, right, gammaIndices, grid), _cb(cb),
          _odd_shifts(oddShifts) {}
  };

  class A2ATaskHalfHalfLocalNoMom : public A2ATaskHalfHalf {
  protected:
    A2AFieldView<cobj> &_gamma;

  public:
    A2ATaskHalfHalfLocalNoMom(A2AFieldView<vobj> &left,
                              A2AFieldView<vobj> &right,
                              A2AFieldView<cobj> &gamma,
                              Vector<Integer> &gammaIndices, GridBase *grid,
                              int cb = Even)
        : A2ATaskHalfHalf(left, right, gammaIndices, grid, cb, false),
          _gamma(gamma) {}

    double getFlops();

    template <typename MatType> void execute(MatType &result, int orthogDir);
  };

  class A2ATaskHalfHalfOneLinkNoMom : public A2ATaskHalfHalf {
  protected:
    A2AFieldView<vColourMatrix> &_links_left, &_links_right;
    A2AStencilView<vobj> &_right_stencil;
    A2AStencilView<vColourMatrix> &_link_stencil;

  public:
    A2ATaskHalfHalfOneLinkNoMom(A2AFieldView<vobj> &left,
                                A2AFieldView<vobj> &right,
                                A2AStencilView<vobj> &rightStencil,
                                A2AFieldView<vColourMatrix> &linksLeft,
                                A2AFieldView<vColourMatrix> &linksRight,
                                A2AStencilView<vColourMatrix> &linkStencil,
                                Vector<Integer> &gammaIndices, GridBase *grid,
                                int cb = Even)
        : A2ATaskHalfHalf(left, right, gammaIndices, grid, cb, true),
          _links_left(linksLeft), _links_right(linksRight),
          _link_stencil(linkStencil), _right_stencil(rightStencil) {}

    double getFlops();

    template <typename MatType> void execute(MatType &result, int orthogDir);
  };

public:
  GridBase *_grid, *_cb_grid;

  LatticeGaugeField *_U;
  const std::vector<StagGamma::SpinTastePair> &_gammas;
  const std::vector<ComplexField> &_mom;

  std::vector<ComplexField> _stag_phase_E, _stag_phase_O;
  std::vector<LatticeColourMatrix> _Umu_E, _Umu_O;

  Vector<Integer> _gamma_indices_local, _gamma_indices_comm;
  std::vector<int> _shift_dirs, _shift_displacements;
  double _flops;

  A2AFieldView<vobj> _view_right_E, _view_right_O, _view_left_E, _view_left_O;
  A2AFieldView<cobj> _view_gamma_E, _view_gamma_O, _view_mom;
  A2AFieldView<vColourMatrix> _view_links_E, _view_links_O;
  A2AStencilView<vobj> _view_stencil_right_E, _view_stencil_right_O;
  A2AStencilView<vColourMatrix> _view_stencil_gauge_E, _view_stencil_gauge_O;

private:
  FermionField *_r_addr_E, *_l_addr_E;

public:
  A2AWorkerMILC() = delete;
  A2AWorkerMILC(GridBase *grid,
                const std::vector<StagGamma::SpinTastePair> &gammas,
                const std::vector<ComplexField> &mom,
                LatticeGaugeField *U = nullptr, GridBase *cbGrid = nullptr);

  ~A2AWorkerMILC();

  double getFlops() { return _flops; }
  void setFlops(double flops) { _flops = flops; }

  void buildLocalPhases();
  void buildGaugeLinks();

  template <typename TensorType> // output: rank 5 tensor, e.g.
                                 // Eigen::Tensor<ComplexD, 5>
                                 void StagMesonFieldNoGlobalSum(
                                     TensorType &mat,
                                     const FermionField *lhs_wi_E,
                                     const FermionField *lhs_wi_O,
                                     const FermionField *rhs_vj_E,
                                     const FermionField *rhs_vj_O,
                                     int orthog_dir,
                                     double *t_kernel = nullptr);
};

#ifndef MF_SUM_ARRAY_MAX
#define MF_SUM_ARRAY_MAX 16
#endif

#if 0
#ifndef GRID_SIMT
#define accelerator_for2dNB_shm(iter1, num1, iter2, num2, nsimd, shm, ...)     \
  accelerator_for2d(iter1, num1, iter2, num2, nsimd, {                         \
    cobj sum[shm];                                                             \
    __VA_ARGS__                                                                \
  });
#else

#define accelerator_for2dNB_shm(iter1, num1, iter2, num2, nsimd, shm, ...)     \
  {                                                                            \
    int nt = acceleratorThreads();                                             \
    typedef uint64_t Iterator;                                                 \
    auto lambda = [=] accelerator(Iterator iter1, Iterator iter2,              \
                                  Iterator lane,                               \
                                  cobj * sum) mutable { __VA_ARGS__; };        \
    dim3 cu_threads(nsimd, acceleratorThreads(), 1);                           \
    dim3 cu_blocks((num1 + nt - 1) / nt, num2, 1);                             \
    ShmLambdaApply<cobj, decltype(lambda)>                                     \
        <<<cu_blocks, cu_threads, shm * sizeof(cobj), computeStream>>>(        \
            num1, num2, nsimd, lambda);                                        \
  }
#endif
#endif

#define A2A_TASK_COMMON()                                                      \
  const int sizeLOut = result.dimension(3);                                    \
  const int sizeROut = result.dimension(4);                                    \
  int sizeR = sizeROut;                                                        \
  int sizeL = sizeLOut;                                                        \
                                                                               \
  const int nGamma = this->_gamma_indices.size();                              \
                                                                               \
  const int simdSize = this->_grid->Nsimd();                                   \
  const int reducedOrthogDimSize = this->_grid->_rdimensions[orthogDir];       \
                                                                               \
  const int nBlocks = this->_grid->_slice_nblock[orthogDir];                   \
                                                                               \
  const int localSpatialVolume = this->_grid->_ostride[orthogDir];             \
                                                                               \
  FermView *viewL_p = &this->_left[0];                                         \
  FermView *viewR_p = &this->_right[0];                                        \
  int localOrthogDimSize = this->_grid->_ldimensions[orthogDir];               \
                                                                               \
  int Nt = this->_grid->GlobalDimensions()[orthogDir];                         \
                                                                               \
  Integer *indexG_p = &this->_gamma_indices[0];                                \
                                                                               \
  int pd = this->_grid->_processors[orthogDir];                                \
  int pc = this->_grid->_processor_coor[orthogDir];                            \
                                                                               \
  auto result_p = result.data();                                               \
                                                                               \
  int localStride = 2;                                                         \
  int orthogSimdSize = this->_grid->_simd_layout[orthogDir];                   \
  Coordinate *icoor_p = &this->_i_coor_container[0];

#define A2A_TASK_HALF_COMMON()                                                 \
  A2A_TASK_COMMON();                                                           \
  sizeL = sizeL / 2;                                                           \
  sizeR = sizeR / 2;                                                           \
  int mult_L = 2;                                                              \
  int mult_R = 2;                                                              \
  bool even = this->_cb == Even;                                               \
  bool oddShifts = this->_odd_shifts;                                          \
  const int MFrvol =                                                           \
      localStride * reducedOrthogDimSize * sizeL * sizeR * nGamma;             \
                                                                               \
  Vector<cobj> shmem(MFrvol);                                                  \
  cobj *shm_p = &shmem[0];                                                     \
                                                                               \
  accelerator_for(r, MFrvol, 1, { shm_p[r] = Zero(); });

#define A2A_KERNEL_COMPUTE_ONELINK_TIMEDIR_NOMOM()                             \
  accelerator_for2d(ls_index, localStride *sizeL, r_index, sizeR, simdSize, {  \
    int site_offset = ls_index % localStride;                                  \
    int l_index = ls_index / localStride;                                      \
                                                                               \
    calcColourMatrix link_ahead, link_behind;                                  \
    calcSpinor left, shift_ahead, shift_behind;                                \
    calcScalar sum[MF_SUM_ARRAY_MAX];                                          \
                                                                               \
    StencilEntry *SE;                                                          \
    int ptype, ss;                                                             \
    int shmem_base = localStride * nGamma * reducedOrthogDimSize *             \
                         (l_index + sizeL * r_index) +                         \
                     site_offset;                                              \
                                                                               \
    for (int rt = 0; rt < reducedOrthogDimSize; rt++) {                        \
                                                                               \
      for (int mu = 0; mu < nGamma; mu++) {                                    \
        sum[mu] = Zero();                                                      \
      }                                                                        \
                                                                               \
      int shmem_idx = localStride * nGamma * rt + shmem_base;                  \
      for (int so = 0; so < localSpatialVolume; so += localStride) {           \
        ss = rt * localSpatialVolume + so + site_offset;                       \
                                                                               \
        left = coalescedRead(viewL_p[l_index][ss]);                            \
                                                                               \
        for (int mu = 0; mu < nGamma; mu++) {                                  \
          link_ahead = coalescedRead(viewGL_p[mu][ss]);                        \
                                                                               \
          SE = stencilR_p->GetEntry(ptype, 2 * mu, ss);                        \
          if (SE->_is_local) {                                                 \
            shift_ahead = coalescedReadPermute(viewR_p[r_index][SE->_offset],  \
                                               ptype, SE->_permute);           \
          } else {                                                             \
            shift_ahead = coalescedRead(                                       \
                bufRight_p[r_index * haloBuffRightSize + SE->_offset]);        \
          }                                                                    \
          acceleratorSynchronise();                                            \
                                                                               \
          SE = stencilR_p->GetEntry(ptype, 2 * mu + 1, ss);                    \
          if (SE->_is_local) {                                                 \
            shift_behind = coalescedReadPermute(viewR_p[r_index][SE->_offset], \
                                                ptype, SE->_permute);          \
          } else {                                                             \
            shift_behind = coalescedRead(                                      \
                bufRight_p[r_index * haloBuffRightSize + SE->_offset]);        \
          }                                                                    \
          acceleratorSynchronise();                                            \
                                                                               \
          SE = stencilG_p[mu].GetEntry(ptype, 0, ss);                          \
          if (SE->_is_local) {                                                 \
            link_behind = adj(coalescedReadPermute(viewGR_p[mu][SE->_offset],  \
                                                   ptype, SE->_permute));      \
          } else {                                                             \
            link_behind =                                                      \
                adj(coalescedRead(bufGauge_p[offsetG_p[mu] + SE->_offset]));   \
          }                                                                    \
          acceleratorSynchronise();                                            \
                                                                               \
          sum[mu] += innerProduct(left, link_ahead * shift_ahead +             \
                                            link_behind * shift_behind);       \
        }                                                                      \
        for (int mu = 0; mu < nGamma; mu++) {                                  \
          coalescedWrite(shm_p[shmem_idx + mu * localStride], sum[mu]);        \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  });

#define A2A_KERNEL_COMPUTE_LOCAL_TIMEDIR_NOMOM()                               \
  accelerator_for2d(ls_index, localStride *sizeL, r_index, sizeR, simdSize, {  \
    int site_offset = ls_index % localStride;                                  \
    int l_index = ls_index / localStride;                                      \
                                                                               \
    int ss, shmem_base = site_offset + localStride * nGamma *                  \
                                           reducedOrthogDimSize *              \
                                           (l_index + sizeL * r_index);        \
    calcScalar gamma_phase, temp_site, sum[MF_SUM_ARRAY_MAX];                  \
                                                                               \
    for (int rt = 0; rt < reducedOrthogDimSize; rt++) {                        \
                                                                               \
      for (int mu = 0; mu < nGamma; mu++) {                                    \
        sum[mu] = Zero();                                                      \
      }                                                                        \
                                                                               \
      int shmem_idx = localStride * nGamma * rt + shmem_base;                  \
      for (int so = 0; so < localSpatialVolume; so += localStride) {           \
        ss = rt * localSpatialVolume + so + site_offset;                       \
        temp_site = innerProduct(coalescedRead(viewL_p[l_index][ss]),          \
                                 coalescedRead(viewR_p[r_index][ss]));         \
                                                                               \
        for (int mu = 0; mu < nGamma; mu++) {                                  \
          gamma_phase = coalescedRead(viewG_p[mu][ss]);                        \
          sum[mu] += gamma_phase * temp_site;                                  \
        }                                                                      \
      }                                                                        \
      for (int mu = 0; mu < nGamma; mu++) {                                    \
        coalescedWrite(shm_p[shmem_idx + mu * localStride], sum[mu]);          \
      }                                                                        \
    }                                                                          \
  });

#define A2A_KERNEL_IO_NOMOM()                                                  \
  accelerator_for2d(li_index, sizeL *orthogSimdSize, r_index, sizeR, 1, {      \
    int simdOffset = li_index % orthogSimdSize;                                \
    int l_index = li_index / orthogSimdSize;                                   \
    ExtractBuffer<Scalar_s> extracted(simdSize);                               \
    cobj strideSum;                                                            \
    scalar_type temp;                                                          \
    int shmem_idx = localStride * reducedOrthogDimSize * nGamma *              \
                    (l_index + sizeL * r_index);                               \
    for (int i = 0; i < nGamma * reducedOrthogDimSize; i++) {                  \
                                                                               \
      int mu = i % nGamma;                                                     \
      int rt = i / nGamma;                                                     \
                                                                               \
      strideSum = Zero();                                                      \
      for (int j = 0; j < localStride; j++) {                                  \
        strideSum += shm_p[shmem_idx + i * localStride + j];                   \
      }                                                                        \
      extract(strideSum, extracted);                                           \
                                                                               \
      temp = scalar_type(0.0);                                                 \
      for (int idx = 0; idx < simdSize; idx++) {                               \
        if (icoor_p[idx][orthogDir] == simdOffset) {                           \
          temp += TensorRemove(extracted[idx]);                                \
        }                                                                      \
        acceleratorSynchronise();                                              \
      }                                                                        \
                                                                               \
      int lt = rt + simdOffset * reducedOrthogDimSize;                         \
      int t_out = lt + pc * localOrthogDimSize;                                \
                                                                               \
      int ij_dx =                                                              \
          mult_R *                                                             \
          (r_index +                                                           \
           sizeR * mult_L * (l_index + sizeL * (t_out + Nt * indexG_p[mu])));  \
                                                                               \
      result_p[ij_dx] += temp;                                                 \
                                                                               \
      if (even && oddShifts) {                                                 \
        result_p[ij_dx + 1] -= temp;                                           \
        result_p[ij_dx + sizeROut + 1] -= temp;                                \
        result_p[ij_dx + sizeROut] += temp;                                    \
      } else if (even) {                                                       \
        result_p[ij_dx + 1] += temp;                                           \
        result_p[ij_dx + sizeROut + 1] += temp;                                \
        result_p[ij_dx + sizeROut] += temp;                                    \
                                                                               \
      } else if (oddShifts) {                                                  \
        result_p[ij_dx + 1] += temp;                                           \
        result_p[ij_dx + sizeROut + 1] -= temp;                                \
        result_p[ij_dx + sizeROut] -= temp;                                    \
      } else {                                                                 \
        result_p[ij_dx + 1] -= temp;                                           \
        result_p[ij_dx + sizeROut + 1] += temp;                                \
        result_p[ij_dx + sizeROut] -= temp;                                    \
      }                                                                        \
      acceleratorSynchronise();                                                \
    }                                                                          \
  });

template <typename FImpl>
double A2AWorkerMILC<FImpl>::A2ATaskHalfHalfLocalNoMom::getFlops() {
  // One complex multiply takes 6 floating point ops (4 mult, 2 add)
  // --> complex inner product is 3 complex mult, 2 complex add = 3*6 + 2*2 = 22
  // double precision floating ops

  // For each vector and at each lattice site:
  //  - one inner product
  //  - For each gamma and momentum
  //    - multiply by gamma phase
  //    - multiply by momentum phase
  //    - sum
  return (22.0 + (6.0 + 2.0) * (this->_gamma_indices.size()));
}

template <typename FImpl>
template <typename MatType>
void A2AWorkerMILC<FImpl>::A2ATaskHalfHalfLocalNoMom::execute(MatType &result,
                                                              int orthogDir) {

  A2A_TASK_HALF_COMMON();

  ComplexView *viewG_p = &this->_gamma[0];

  assert(orthogDir ==
         Tdir); // This kernel assumes lattice is coalesced over time slices
  assert(nBlocks == 1);

  A2A_KERNEL_COMPUTE_LOCAL_TIMEDIR_NOMOM();

  A2A_KERNEL_IO_NOMOM();
}

template <typename FImpl>
double A2AWorkerMILC<FImpl>::A2ATaskHalfHalfOneLinkNoMom::getFlops() {
  // matrix*vector = 3 inner products
  // current code:
  // innerProduct(left,link_ahead*shift_ahead+adj(link_behind)*shift_behind)
  //  = matrix*vector + matrix* vector --> inner product = 7 inner products and
  //  1 complex sum

  // For each vector, each gamma, and at each lattice site:
  //  - one inner product
  //  - two su(3) matrix*vector ops
  //  - one complex sum
  return ((7 * 22.0 + 2.0) * this->_gamma_indices.size());
}

template <typename FImpl>
template <typename MatType>
void A2AWorkerMILC<FImpl>::A2ATaskHalfHalfOneLinkNoMom::execute(MatType &result,
                                                                int orthogDir) {

  A2A_TASK_HALF_COMMON();

  assert(orthogDir ==
         Tdir); // This kernel assumes lattice is coalesced over time slices
  assert(nBlocks == 1);

  // Pointers for accelerator indexing
  GaugeView *viewGL_p = &this->_links_left[0];
  GaugeView *viewGR_p = &this->_links_right[0];
  GaugeStencilView *stencilG_p =
      &this->_link_stencil[0]; // Gauge tencil for shifted links
  FermStencilView *stencilR_p =
      &this->_right_stencil[0]; // Gauge tencil for shifted links

  vobj *bufRight_p = &this->_right_stencil.getBuffer(
      0); // buffer for shifted kets in halo region
  vColourMatrix *bufGauge_p = &this->_link_stencil.getBuffer(
      0); // buffer for shifted links in halo region

  Integer *offsetG_p = &this->_link_stencil.getOffset(
      0); // buffer offsets for shifted links in halo region
  int haloBuffRightSize = this->_right_stencil.getOffset(1);

  A2A_KERNEL_COMPUTE_ONELINK_TIMEDIR_NOMOM();

  A2A_KERNEL_IO_NOMOM();
}

template <typename FImpl>
A2AWorkerMILC<FImpl>::A2AWorkerMILC(
    GridBase *grid, const std::vector<StagGamma::SpinTastePair> &gammas,
    const std::vector<ComplexField> &mom, LatticeGaugeField *U,
    GridBase *cbGrid)
    : _grid(grid), _gammas(gammas), _mom(mom), _U(U), _cb_grid(cbGrid) {

  _l_addr_E = nullptr;
  _r_addr_E = nullptr;

  StagGamma spinTaste;
  if (_U != nullptr) {
    spinTaste.setGaugeField(*_U);
  }

  // Organize gammas into local/non-local
  for (int i = 0; i < _gammas.size(); i++) {

    spinTaste.setSpinTaste(_gammas[i]);

    int shift = (spinTaste._spin ^ spinTaste._taste);
    if (shift != 0) {

      assert(_U != nullptr);

      _gamma_indices_comm.push_back(i);

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
    } else {
      _gamma_indices_local.push_back(i);
    }
  }

  if (_gamma_indices_local.size() > 0) {
    buildLocalPhases();
  }

  _view_mom.reserve(_mom.size());
  for (int i = 0; i < _mom.size(); ++i) {
    _view_mom.addView(_mom[i]);
  }

  if (_gamma_indices_comm.size() > 0) {
    double t0 = usecond();
    buildGaugeLinks();
    double t1 = usecond();
    std::cout << GridLogPerformance
              << " MesonField one link timings: build link fields:"
              << (t1 - t0) / 1000 << "ms" << std::endl;
  }
}

template <typename FImpl> A2AWorkerMILC<FImpl>::~A2AWorkerMILC() {
  _view_left_E.closeViews();
  _view_left_O.closeViews();
  _view_right_E.closeViews();
  _view_right_O.closeViews();
  _view_stencil_right_E.closeViews();
  _view_stencil_right_O.closeViews();
  _view_links_E.closeViews();
  _view_links_O.closeViews();
  _view_stencil_gauge_E.closeViews();
  _view_stencil_gauge_O.closeViews();
  _view_gamma_E.closeViews();
  _view_gamma_O.closeViews();
  _view_mom.closeViews();
}
template <typename FImpl> void A2AWorkerMILC<FImpl>::buildLocalPhases() {
  int nGamma_local = _gamma_indices_local.size();

  _stag_phase_E.resize(nGamma_local, _cb_grid);
  _stag_phase_O.resize(nGamma_local, _cb_grid);

  { // Set up staggered phases
    StagGamma spinTaste;
    ComplexField temp(_grid);
    int mu;

    _view_gamma_E.reserve(nGamma_local);
    _view_gamma_O.reserve(nGamma_local);
    for (int i = 0; i < nGamma_local; i++) {
      mu = _gamma_indices_local[i];

      temp = 1.0;

      spinTaste.setSpinTaste(_gammas[mu]);
      spinTaste.applyPhase(temp, temp); // store spin-taste phase

      pickCheckerboard(Even, _stag_phase_E[i], temp);
      pickCheckerboard(Odd, _stag_phase_O[i], temp);

      _view_gamma_E.addView(_stag_phase_E[i]);
      _view_gamma_O.addView(_stag_phase_O[i]);
    }
  }
}

template <typename FImpl> void A2AWorkerMILC<FImpl>::buildGaugeLinks() {
  int mu, nGamma_comm = _gamma_indices_comm.size();

  _Umu_E.resize(nGamma_comm, _cb_grid);
  _Umu_O.resize(nGamma_comm, _cb_grid);

  StagGamma spinTaste;
  LatticeColourMatrix Umu_temp(_U->Grid());

  for (int i = 0; i < nGamma_comm; i++) {
    mu = _gamma_indices_comm[i];
    spinTaste.setSpinTaste(_gammas[mu]);

    Umu_temp = PeekIndex<LorentzIndex>(
        *_U, _shift_dirs[2 * i]); // Store full lattice links in shift direction

    spinTaste.applyPhase(Umu_temp, Umu_temp); // store spin-taste phase

    pickCheckerboard(Even, _Umu_E[i], Umu_temp);
    pickCheckerboard(Odd, _Umu_O[i], Umu_temp);
  }

  double t0 = usecond();
  int size = _shift_dirs.size() / 2;

  for (int i = 0; i < size; ++i) {
    auto ptr = std::move(std::unique_ptr<GaugeStencil>(
        new GaugeStencil(_cb_grid, 1, Even, {_shift_dirs[2 * i]}, {-1})));
    _view_stencil_gauge_E.addStencil(ptr);
    _view_stencil_gauge_E.append(_Umu_E[i]);

    ptr = std::move(std::unique_ptr<GaugeStencil>(
        new GaugeStencil(_cb_grid, 1, Odd, {_shift_dirs[2 * i]}, {-1})));
    _view_stencil_gauge_O.addStencil(ptr);
    _view_stencil_gauge_O.append(_Umu_O[i]);
  }
  _view_stencil_gauge_E.openViews();
  _view_stencil_gauge_O.openViews();

  double t1 = usecond();

  std::cout << GridLogPerformance
            << " MesonField one link timings: gauge comms:" << (t1 - t0) / 1000
            << "ms" << std::endl;

  _view_links_E.reserve(nGamma_comm);
  _view_links_O.reserve(nGamma_comm);
  for (int i = 0; i < nGamma_comm; ++i) {
    _view_links_E.addView(_Umu_E[i]);
    _view_links_O.addView(_Umu_O[i]);
  }
}

template <class FImpl>
template <typename TensorType>
void A2AWorkerMILC<FImpl>::StagMesonFieldNoGlobalSum(
    TensorType &mat, const FermionField *lhs_wi_E, const FermionField *lhs_wi_O,
    const FermionField *rhs_vj_E, const FermionField *rhs_vj_O, int orthog_dir,
    double *t_kernel) {
  assert(_cb_grid->CheckerBoarded(orthog_dir) != 1);

  bool checkerL = lhs_wi_E[0].Grid()->_isCheckerBoarded;
  bool checkerR = rhs_vj_E[0].Grid()->_isCheckerBoarded;

  int sizeL = mat.dimension(3) / 2;
  int sizeR = mat.dimension(4) / 2;

  double t0 = usecond();

  if (_l_addr_E != lhs_wi_E) {
    _l_addr_E = const_cast<FermionField *>(lhs_wi_E);
    _view_left_E.closeViews();
    _view_left_O.closeViews();

    _view_left_E.reserve(sizeL);
    _view_left_O.reserve(sizeL);
    for (int i = 0; i < sizeL; ++i) {
      _view_left_E.addView(lhs_wi_E[i]);
      _view_left_O.addView(lhs_wi_O[i]);
    }
  }

  if (_r_addr_E != rhs_vj_E) {
    _r_addr_E = const_cast<FermionField *>(rhs_vj_E);
    _view_right_E.closeViews();
    _view_right_O.closeViews();
    _view_stencil_right_E.closeViews();
    _view_stencil_right_O.closeViews();

    if (_gamma_indices_comm.size() > 0) {
      auto ptr = std::move(std::unique_ptr<FermStencil>(
          new FermStencil(_cb_grid, _shift_dirs.size(), Even, _shift_dirs,
                          _shift_displacements)));
      _view_stencil_right_E.addStencil(ptr);

      ptr = std::move(std::unique_ptr<FermStencil>(
          new FermStencil(_cb_grid, _shift_dirs.size(), Odd, _shift_dirs,
                          _shift_displacements)));
      _view_stencil_right_O.addStencil(ptr);
      for (int i = 0; i < sizeR; ++i) {
        _view_stencil_right_E.append(rhs_vj_E[i]);
        _view_stencil_right_O.append(rhs_vj_O[i]);
      }
      _view_stencil_right_E.openViews();
      _view_stencil_right_O.openViews();
    }

    _view_right_E.reserve(sizeR);
    _view_right_O.reserve(sizeR);
    for (int i = 0; i < sizeR; ++i) {
      _view_right_E.addView(rhs_vj_E[i]);
      _view_right_O.addView(rhs_vj_O[i]);
    }
  }
  double t1 = usecond();
  std::cout << GridLogPerformance
            << " MesonField timings: left/right views+comms:"
            << (t1 - t0) / 1000 << "ms" << std::endl;

  if (t_kernel)
    *t_kernel = -usecond();

  setFlops(0.0);
  // Run any nonlocal gamma operators
  if (_gamma_indices_comm.size() > 0) {
    if (checkerL && checkerR) {
      A2ATaskHalfHalfOneLinkNoMom task_e(_view_left_E, _view_right_O,
                                         _view_stencil_right_O, _view_links_E,
                                         _view_links_O, _view_stencil_gauge_O,
                                         _gamma_indices_comm, _cb_grid, Even);
      task_e.execute(mat, orthog_dir);

      A2ATaskHalfHalfOneLinkNoMom task_o(_view_left_O, _view_right_E,
                                         _view_stencil_right_E, _view_links_O,
                                         _view_links_E, _view_stencil_gauge_E,
                                         _gamma_indices_comm, _cb_grid, Odd);
      task_o.execute(mat, orthog_dir);

      setFlops(task_e.getFlops());
    } else if (checkerL) {
      assert(0);
    } else if (checkerR) {
      assert(0);
    } else {
      assert(0);
    }
  }

  // Run any local gamma operators
  if (_gamma_indices_local.size() > 0) {

    if (checkerL && checkerR) {
      A2ATaskHalfHalfLocalNoMom task_e(_view_left_E, _view_right_E,
                                       _view_gamma_E, _gamma_indices_local,
                                       _cb_grid, Even);
      task_e.execute(mat, orthog_dir);

      A2ATaskHalfHalfLocalNoMom task_o(_view_left_O, _view_right_O,
                                       _view_gamma_O, _gamma_indices_local,
                                       _cb_grid, Odd);
      task_o.execute(mat, orthog_dir);

      setFlops(getFlops() + task_e.getFlops());
    } else if (checkerL) {
      assert(0);
    } else if (checkerR) {
      assert(0);
    } else {
      assert(0);
    }
  }
  if (t_kernel)
    *t_kernel += usecond();

  int resultStride = mat.dimension(2) * mat.dimension(3) * mat.dimension(4);
  int nGamma = mat.dimension(0) * mat.dimension(1);

  _grid->GlobalSumVector(mat.data(), nGamma * resultStride);
}

NAMESPACE_END(Grid);

#undef A2A_TASK_COMMON
#undef A2A_TASK_HALF_COMMON
#undef A2A_KERNEL_COMPUTE_ONELINK_TIMEDIR_NOMOM
#undef A2A_KERNEL_COMPUTE_LOCAL_TIMEDIR_NOMOM
#undef A2A_KERNEL_IO_NOMOM
