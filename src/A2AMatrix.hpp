/*
 * A2AMatrixMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
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
#ifndef A2A_Matrix__hpp_
#define A2A_Matrix__hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Grid/Eigen/unsupported/CXX11/Tensor>
#ifdef USE_MKL
#include "mkl.h"
#include "mkl_cblas.h"
#endif

#ifndef HADRONS_A2AM_NAME
#define HADRONS_A2AM_NAME "a2aMatrix"
#endif

#ifndef HADRONS_A2AM_IO_TYPE
#define HADRONS_A2AM_IO_TYPE ComplexF
#endif

#define HADRONS_A2AM_PARALLEL_IO

BEGIN_HADRONS_NAMESPACE

// general A2A matrix set based on Eigen tensors and Grid-allocated memory
// Dimensions:
//   0 - ext - external field (momentum, EM field, ...)
//   1 - str - spin-color structure
//   2 - t   - timeslice
//   3 - i   - left  A2A mode index
//   4 - j   - right A2A mode index
template <typename T>
using A2AMatrixSet = Eigen::TensorMap<Eigen::Tensor<T, 5, Eigen::RowMajor>>;

template <typename T>
using A2AMatrix = Eigen::Matrix<T, -1, -1, Eigen::RowMajor>;

template <typename T>
using A2AMatrixTr = Eigen::Matrix<T, -1, -1, Eigen::ColMajor>;

// Per-side staged-input content stamp: identifies WHICH block's staging a
// kernel call carries. A side's staged contents are a pure function of
// (execute call, block coordinate): the block loop is the sole writer of
// the staged buffers and staging is deterministic per coordinate -- pair
// sources refill their grow-only buffers per requested range, the SwapFn
// path refills _lowBuf and re-swaps evec slots per block (restoring them
// after), high modes pick a fixed vector element per coordinate, and the
// i==j alias hands over the right side's mode-i staging (the same content
// get() would produce). Addresses are NOT content keys: pair-source
// buffers and _lowBuf are deliberately stable across refills. epoch
// disambiguates repeated execute() calls on the same kernel.
struct A2ASideStamp {
  unsigned long long epoch = 0;
  int block = -1;
  bool operator==(const A2ASideStamp &o) const {
    return epoch == o.epoch && block == o.block;
  }
  bool operator!=(const A2ASideStamp &o) const { return !(*this == o); }
};

/******************************************************************************
 *                      Abstract class for A2A kernels                        *
 ******************************************************************************/
template <typename T, typename Field> class A2AKernelMILC {
public:
  A2AKernelMILC(void) = default;
  virtual ~A2AKernelMILC(void) = default;
  virtual void operator()(A2AMatrixSet<T> &m, const Field *left_e,
                          const Field *left_o, const Field *right_e,
                          const Field *right_o) = 0;
  // Stamped overload: kernels that track per-side staged-input validity
  // key on (epoch, block) content stamps; the default forwards to the
  // pointer-based entry so untracked kernels work unchanged.
  virtual void operator()(A2AMatrixSet<T> &m, const Field *left_e,
                          const Field *left_o, const Field *right_e,
                          const Field *right_o, const A2ASideStamp &lstamp,
                          const A2ASideStamp &rstamp) {
    (*this)(m, left_e, left_o, right_e, right_o);
  }
  virtual double flops(const unsigned int blockSizei,
                       const unsigned int blockSizej, int cbDiv = 1) = 0;
  virtual double bytes(const unsigned int blockSizei,
                       const unsigned int blockSizej) = 0;
  virtual double kernelTime() = 0;
  virtual double globalSumTime() = 0;
};

/******************************************************************************
 *              Abstract class for CB low-mode pair sources                   *
 ******************************************************************************/
// Supplier of checkerboarded low-mode (even, odd) fermion arrays for the
// A2AMatrix block loop (A2AMatrixBlockComputationMILC::execute). For a
// request [start, start+n) of eigenvectors, get() fills source-managed,
// grow-only buffers with the even-CB and odd-CB parts of each evec in the
// range and returns pointers to the two length-n arrays -- first always
// Checkerboard()==Even content, second always Checkerboard()==Odd. The
// eigenpack itself is never mutated (unlike the legacy SwapFn mechanism,
// which destructively swaps evec content in place and restores it after
// each block). Interface consumed through envCreateDerived/envGet, the
// same base-typed product pattern as FMat/noise modules; A2AKernelMILC
// above is the framework-interface precedent.
template <typename Field> class A2ALowModePairSourceMILC {
public:
  A2ALowModePairSourceMILC(void) = default;
  virtual ~A2ALowModePairSourceMILC(void) = default;
  // (even array, odd array) for evecs [start, start+n); the buffers are
  // owned by the source and valid until the next get() call.
  virtual std::pair<const Field *, const Field *> get(const int start,
                                                      const int n) = 0;
};

/******************************************************************************
 *                  Class to handle A2A matrix block HDF5 I/O                 *
 ******************************************************************************/
template <typename T> class A2AMatrixIoMILC {
public:
  // constructors
  A2AMatrixIoMILC(void) = default;
  A2AMatrixIoMILC(std::string filename, std::string dataname,
                  const unsigned int nt, const unsigned int ni = 0,
                  const unsigned int nj = 0, const unsigned int ni_start = 0,
                  const unsigned int nj_start = 0);
  // destructor
  ~A2AMatrixIoMILC(void) = default;
  // access
  unsigned int getNi(void) const;
  unsigned int getNj(void) const;
  unsigned int getNt(void) const;
  size_t getSize(void) const;
  // file allocation
  template <typename MetadataType>
  void initFile(const MetadataType &d, const unsigned int chunkSize_i,
                const unsigned int chunkSize_j);
  // block I/O
  void saveBlock(const T *data, const unsigned int i, const unsigned int j,
                 const unsigned int blockSizei, const unsigned int blockSizej);
  void saveBlock(const A2AMatrixSet<T> &m, const unsigned int ext,
                 const unsigned int str, const unsigned int i,
                 const unsigned int j);
  template <template <class> class Vec, typename VecT>
  void load(Vec<VecT> &v, double *tRead = nullptr, GridBase *grid = nullptr);

private:
  std::string _filename{""}, _dataname{""};
  unsigned int _nt{0}, _ni{0}, _nj{0}, _ni_start{0}, _nj_start{0};
};
/******************************************************************************
 *                  Wrapper for A2A matrix block computation                  *
 ******************************************************************************/
template <typename T, typename Field, typename MetadataType, typename TIo = T>
class A2AMatrixBlockComputationMILC {
private:
  struct IoHelper {
    A2AMatrixIoMILC<TIo> io;
    MetadataType md;
    unsigned int e, s, i, j;
  };
  typedef std::function<std::string(const unsigned int, const unsigned int)>
      FilenameFn;
  typedef std::function<MetadataType(const unsigned int, const unsigned int)>
      MetadataFn;
  typedef std::function<void(int)> SwapFn;

public:
  // constructor
  A2AMatrixBlockComputationMILC(GridBase *grid, const unsigned int orthogDim,
                                const unsigned int next,
                                const unsigned int nstr,
                                const unsigned int blockSize,
                                TimerArray *tArray = nullptr);
  // execution
  void execute(const std::vector<Field> &left, const std::vector<Field> &right,
               A2AKernelMILC<T, Field> &kernel, const FilenameFn &ionameFn,
               const FilenameFn &filenameFn, const MetadataFn &metadataFn,
               std::vector<Field> *evecs = nullptr,
               const std::vector<ComplexD> &evals = {},
               const SwapFn *swapEvecCheckerFn = nullptr,
               A2ALowModePairSourceMILC<Field> *pairSourceLeft = nullptr,
               A2ALowModePairSourceMILC<Field> *pairSourceRight = nullptr);

private:
  // I/O handler
  void saveBlock(const A2AMatrixSet<TIo> &m, IoHelper &h);

private:
  TimerArray *_tArray;
  GridBase *_grid;
  unsigned int _orthogDim, _nt, _next, _nstr, _blockSize, _min_i, _min_j;
  std::vector<IoHelper> _nodeIo;
  std::vector<Field> _lowBuf_i, _lowBuf_j;
  // Content-generation counter for per-side stamps (the Computation
  // object persists per module via envGetTmp, so the bump separates
  // repeated execute() calls on the same kernel).
  unsigned long long _stampEpoch = 0;
};

/******************************************************************************
 *                     A2AMatrixIoMILC template implementation *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T>
A2AMatrixIoMILC<T>::A2AMatrixIoMILC(std::string filename, std::string dataname,
                                    const unsigned int nt,
                                    const unsigned int ni,
                                    const unsigned int nj,
                                    const unsigned int ni_start,
                                    const unsigned int nj_start)
    : _filename(filename), _dataname(dataname), _nt(nt), _ni(ni), _nj(nj),
      _ni_start(ni_start), _nj_start(nj_start) {}

// access //////////////////////////////////////////////////////////////////////
template <typename T> unsigned int A2AMatrixIoMILC<T>::getNt(void) const {
  return _nt;
}

template <typename T> unsigned int A2AMatrixIoMILC<T>::getNi(void) const {
  return _ni;
}

template <typename T> unsigned int A2AMatrixIoMILC<T>::getNj(void) const {
  return _nj;
}

template <typename T> size_t A2AMatrixIoMILC<T>::getSize(void) const {
  return _nt * _ni * _nj * sizeof(T);
}

// file allocation /////////////////////////////////////////////////////////////
template <typename T>
template <typename MetadataType>
void A2AMatrixIoMILC<T>::initFile(const MetadataType &d,
                                  const unsigned int chunkSize_i,
                                  const unsigned int chunkSize_j) {
#ifdef HAVE_HDF5
  std::vector<hsize_t> dim = {static_cast<hsize_t>(_nt),
                              static_cast<hsize_t>(_ni),
                              static_cast<hsize_t>(_nj)},
                       chunk = {static_cast<hsize_t>(_nt),
                                static_cast<hsize_t>(chunkSize_i),
                                static_cast<hsize_t>(chunkSize_j)};
  H5NS::DataSpace dataspace(dim.size(), dim.data());
  H5NS::DataSet dataset;
  H5NS::DSetCreatPropList plist;

  // create empty file just with metadata
  {
    Hdf5Writer writer(_filename);
    write(writer, _dataname, d);
  }

  // create the dataset
  Hdf5Reader reader(_filename, false);

  push(reader, _dataname);
  auto &group = reader.getGroup();
  plist.setChunk(chunk.size(), chunk.data());
  plist.setFletcher32();
  dataset = group.createDataSet(HADRONS_A2AM_NAME, Hdf5Type<T>::type(),
                                dataspace, plist);
#else
  HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

// block I/O ///////////////////////////////////////////////////////////////////
template <typename T>
void A2AMatrixIoMILC<T>::saveBlock(const T *data, const unsigned int i,
                                   const unsigned int j,
                                   const unsigned int blockSizei,
                                   const unsigned int blockSizej) {
#ifdef HAVE_HDF5
  Hdf5Reader reader(_filename, false);
  std::vector<hsize_t> count = {_nt, blockSizei, blockSizej},
                       offset = {0, static_cast<hsize_t>(i),
                                 static_cast<hsize_t>(j)},
                       stride = {1, 1, 1}, block = {1, 1, 1};
  H5NS::DataSpace memspace(count.size(), count.data()), dataspace;
  H5NS::DataSet dataset;
  push(reader, _dataname);
  auto &group = reader.getGroup();
  dataset = group.openDataSet(HADRONS_A2AM_NAME);
  dataspace = dataset.getSpace();
  dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                            stride.data(), block.data());
  dataset.write(data, Hdf5Type<T>::type(), memspace, dataspace);
#else
  HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

template <typename T>
void A2AMatrixIoMILC<T>::saveBlock(const A2AMatrixSet<T> &m,
                                   const unsigned int ext,
                                   const unsigned int str, const unsigned int i,
                                   const unsigned int j) {
  unsigned int blockSizei = m.dimension(3);
  unsigned int blockSizej = m.dimension(4);
  unsigned int nstr = m.dimension(1);
  size_t offset = (ext * nstr + str) * _nt * blockSizei * blockSizej;

  saveBlock(m.data() + offset, i, j, blockSizei, blockSizej);
}

template <typename T>
template <template <class> class Vec, typename VecT>
void A2AMatrixIoMILC<T>::load(Vec<VecT> &v, double *tRead, GridBase *grid) {
#ifdef HAVE_HDF5
  std::vector<hsize_t> hdim;
  H5NS::DataSet dataset;
  H5NS::DataSpace dataspace;
  H5NS::CompType datatype;

  unsigned int myRank = 0, nRank = 1;
  if (grid) {
    myRank = grid->ThisRank(), nRank = grid->RankCount();
  }

  Hdf5Reader reader(_filename);
  push(reader, _dataname);
  auto &group = reader.getGroup();
  dataset = group.openDataSet(HADRONS_A2AM_NAME);
  datatype = dataset.getCompType();
  dataspace = dataset.getSpace();
  hdim.resize(dataspace.getSimpleExtentNdims());
  dataspace.getSimpleExtentDims(hdim.data());

  if ((_nt * _ni * _nj != 0) and
      ((hdim[0] < _nt) or (hdim[1] < (_ni + _ni_start)) or
       (hdim[2] < (_nj + _nj_start)))) {
    HADRONS_ERROR(
        Size, "all-to-all matrix size mismatch (got " +
                  std::to_string(hdim[0]) + "x" + std::to_string(hdim[1]) +
                  "x" + std::to_string(hdim[2]) + ", expected at least" +
                  std::to_string(_nt) + "x" + std::to_string(_ni + _ni_start) +
                  "x" + std::to_string(_nj + _nj_start));
  } else if (_ni * _nj == 0) {
    if (hdim[0] != _nt) {
      HADRONS_ERROR(Size, "all-to-all time size mismatch (got " +
                              std::to_string(hdim[0]) + ", expected " +
                              std::to_string(_nt) + ")");
    }
    _ni = hdim[1] - _ni_start;
    _nj = hdim[2] - _nj_start;
  }

  std::vector<hsize_t> count = {1, static_cast<hsize_t>(_ni),
                                static_cast<hsize_t>(_nj)},
                       stride = {1, 1, 1}, block = {1, 1, 1},
                       memCount = {static_cast<hsize_t>(_ni),
                                   static_cast<hsize_t>(_nj)};
  H5NS::DataSpace memspace(memCount.size(), memCount.data());

  std::cout << "Loading timeslice";
  std::cout.flush();
  *tRead = 0.;
  if (grid) {
    Vector<A2AMatrix<T>> buf(_nt, A2AMatrix<T>(_ni, _nj));

    int broadcastSize = sizeof(T) * buf[0].size();
    for (int t = myRank; t < _nt; t += nRank) {

      std::vector<hsize_t> offset = {static_cast<hsize_t>(t), _ni_start,
                                     _nj_start};

      std::cout << " " << t;
      std::cout.flush();

      dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                stride.data(), block.data());

      if (tRead)
        *tRead -= usecond();
      dataset.read(buf[t].data(), datatype, memspace, dataspace);
      if (tRead)
        *tRead += usecond();
    }
    grid->Barrier();
    for (int t = 0; t < _nt; t++) {
      int rank = t % nRank;
      grid->Broadcast(rank, buf[t].data(), broadcastSize);
      // grid->SendToRecvFrom(buf[idx].data(),grid->BossRank(),buf[idx].data(),rank,broadcastSize);
    }
    grid->Barrier();
    for (int t = 0; t < _nt; t++) {
      v[t] = buf[t].template cast<VecT>();
    }
  } else {

    A2AMatrix<T> buf(_ni, _nj);
    int broadcastSize = sizeof(T) * buf.size();
    for (unsigned int tp1 = _nt; tp1 > 0; --tp1) {
      unsigned int t = tp1 - 1;
      std::vector<hsize_t> offset = {static_cast<hsize_t>(t), _ni_start,
                                     _nj_start};

      if (t % 10 == 0) {
        std::cout << " " << t;
        std::cout.flush();
      }
      dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                stride.data(), block.data());
      if (tRead)
        *tRead -= usecond();
      dataset.read(buf.data(), datatype, memspace, dataspace);
      if (tRead)
        *tRead += usecond();
      v[t] = buf.template cast<VecT>();
    }
  }

  std::cout << std::endl;
#else
  HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

/******************************************************************************
 *               A2AMatrixBlockComputation template implementation            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
A2AMatrixBlockComputationMILC<T, Field, MetadataType, TIo>::
    A2AMatrixBlockComputationMILC(GridBase *grid, const unsigned int orthogDim,
                                  const unsigned int next,
                                  const unsigned int nstr,
                                  const unsigned int blockSize,
                                  TimerArray *tArray)
    : _grid(grid), _nt(grid->GlobalDimensions()[orthogDim]),
      _orthogDim(orthogDim), _next(next), _nstr(nstr), _blockSize(blockSize),
      _tArray(tArray) {}

#undef START_TIMER
#undef STOP_TIMER
#undef GET_TIMER

#define START_TIMER(name)                                                      \
  if (_tArray)                                                                 \
  _tArray->startTimer(name)
#define STOP_TIMER(name)                                                       \
  if (_tArray)                                                                 \
  _tArray->stopTimer(name)
#define GET_TIMER(name) ((_tArray != nullptr) ? _tArray->getDTimer(name) : 0.)

// execution ///////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputationMILC<T, Field, MetadataType, TIo>::execute(
    const std::vector<Field> &left, const std::vector<Field> &right,
    A2AKernelMILC<T, Field> &kernel, const FilenameFn &ionameFn,
    const FilenameFn &filenameFn, const MetadataFn &metadataFn,
    std::vector<Field> *evecs, const std::vector<ComplexD> &evals,
    const SwapFn *swapEvecCheckerFn,
    A2ALowModePairSourceMILC<Field> *pairSourceLeft,
    A2ALowModePairSourceMILC<Field> *pairSourceRight) {
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over _blockSize factors
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////

  _stampEpoch++; // this call's stamps are a new content generation

  std::vector<T> mCache;
  mCache.resize(_nt * _next * _nstr * _blockSize * _blockSize);
  // LOG(Message) << "cache bytes: " << mCache_bytes << std::endl;

  MemoryManager::PrintBytes();

  std::vector<TIo> mBuf;

  const int nPairSources =
      (pairSourceLeft != nullptr) + (pairSourceRight != nullptr);
  if (swapEvecCheckerFn != nullptr && nPairSources != 0) {
    HADRONS_ERROR(Argument,
                  "A2AMatrixBlockComputationMILC::execute: pass either a "
                  "SwapFn or low-mode pair sources, not both");
  }
  if (nPairSources == 1) {
    HADRONS_ERROR(
        Argument,
        "A2AMatrixBlockComputationMILC::execute: checkerboarded low modes "
        "need BOTH a left and a right pair source");
  }
  if (nPairSources == 2 && pairSourceLeft == pairSourceRight) {
    HADRONS_ERROR(
        Argument,
        "A2AMatrixBlockComputationMILC::execute: left and right pair "
        "sources must be distinct objects (the block loop holds the "
        "right-side result across the whole left-side loop, so the two "
        "sides cannot share one bank)");
  }
  bool checkerboarded_low =
      (swapEvecCheckerFn != nullptr || nPairSources == 2);
  int Ncb = checkerboarded_low
                ? 2
                : 1; // Ncb == 2 if the low modes are checkerboarded

  RealD norm = 1.0;

  int N_low = 0;
  if (evecs != nullptr) {
    norm = 1.0 / ::sqrt(norm2(evecs->at(0))); // Calculate norm of eigenvectors
    N_low =
        Ncb *
        evecs->size(); // N_low is the number of evecs for M + evecs for Mdag
  }

  int N_i = left.size();  // Total number of bra vectors to contract
  int N_j = right.size(); // Total number of ket vectors to contract

  // If eigs were passed in but either left or right are empty, just do cross
  // terms (don't calculate low-low)
  bool skip_low_left = false, skip_low_right = false;
  if (N_low != 0) {
    if (left.size() == 0 && right.size() != 0) {
      skip_low_right = true;
      N_i = N_low;
    } else if (left.size() != 0 && right.size() == 0) {
      skip_low_left = true;
      N_j = N_low;
    } else {
      N_i += N_low;
      N_j += N_low;
    }
  }

  _min_i = _min_j = _blockSize;

  if (N_i < _blockSize)
    _min_i = N_i;
  if (N_j < _blockSize)
    _min_j = N_j;

  double flops, bytes, t_kernel, t_gsum;
  double nodes = _grid->NodeCount();

  if (checkerboarded_low) {
    norm = norm / ::sqrt(2); // Reduce checkerboarded norm to 1/sqrt(2)

    if (_blockSize % 2 != 0) {
      HADRONS_ERROR(
          Implementation,
          "Blocksize must be divisible by 2 for checkerboarded low modes");
    }

    if (swapEvecCheckerFn != nullptr) {
      if (!skip_low_left)
        _lowBuf_i.resize(
            _blockSize / 2,
            evecs->at(0).Grid()); // storage for caching checkerboards
      if (!skip_low_right)
        _lowBuf_j.resize(_blockSize / 2, evecs->at(0).Grid());
    }
  }

  int NBlock_i = N_i / _blockSize +
                 (((N_i % _blockSize) != 0)
                      ? 1
                      : 0); // Round up on the number of blocks to compute
  int NBlock_j = N_j / _blockSize + (((N_j % _blockSize) != 0) ? 1 : 0);

  bool low_i, low_j;
  int i, j, evec_i, evec_j, N_ii, N_jj;

  // Per-block untimed-work attribution (A2A host-loops instrumentation).
  // TimerArray::getDTimer is CUMULATIVE (never reset-on-read), so every
  // per-block value below uses the negate-accumulate idiom of ioTime
  // (v = -GET_TIMER(name) before START_TIMER, v += GET_TIMER(name)
  // after STOP_TIMER). The stage timer is open across block boundaries:
  // a block's "stage" is everything from the previous block's untimed
  // line (or loop entry) to this block's mCache-zero region -- evec
  // restore, loop advance and L/R tile staging land in the NEXT block's
  // stage.
  double zeroUs, stageUs, wrapUs, normUs, barrierUs;
  stageUs = -GET_TIMER("untimed: stage");
  START_TIMER("untimed: stage");

  j = 0, evec_j = 0;
  while (j < N_j) { // While we still have bra vectors to contract

    low_j = j < N_low && !skip_low_right;

    const Field *r_temp_e, *r_temp_o;
    if (low_j) {
      N_jj = MIN(N_low - j, _blockSize);

      if (pairSourceRight != nullptr) {
        auto pairs = pairSourceRight->get(evec_j, N_jj / 2);
        r_temp_e = pairs.first;
        r_temp_o = pairs.second;
      } else if (checkerboarded_low) {
        for (int idxj = evec_j; idxj < (MIN(N_low, j + N_jj) / 2); idxj++) {
          _lowBuf_j[idxj - evec_j] = evecs->at(
              idxj); // Cache original evecs to avoid excessive Meooe ops.
          (*swapEvecCheckerFn)(idxj); // Swap original evec checkerboard to
                                      // complementary checkerboard.
        }
        if (_lowBuf_j[0].Checkerboard() == Even) {
          r_temp_e = &_lowBuf_j[0];
          r_temp_o = &evecs->at(evec_j);
        } else {
          r_temp_o = &_lowBuf_j[0];
          r_temp_e = &evecs->at(evec_j);
        }
      } else {
        r_temp_e = &evecs->at(evec_j);
        r_temp_o = nullptr;
      }
    } else {
      N_jj = MIN(N_j - j, _blockSize);

      if (skip_low_right)
        r_temp_e = &right[j];
      else
        r_temp_e = &right[j - N_low];

      r_temp_o = nullptr;
    }

    i = 0, evec_i = 0;
    while (i < N_i) { // While we still have ket vectors to contract

      low_i = i < N_low && !skip_low_left;

      const Field *l_temp_e, *l_temp_o;
      if (low_i) {
        N_ii = MIN(N_low - i, _blockSize);

        if (low_j && i == j) {
          l_temp_e = r_temp_e;
          if (checkerboarded_low)
            l_temp_o = r_temp_o;
          else
            l_temp_o = nullptr;

        } else if (pairSourceLeft != nullptr) {
          auto pairs = pairSourceLeft->get(evec_i, N_ii / 2);
          l_temp_e = pairs.first;
          l_temp_o = pairs.second;
        } else if (checkerboarded_low) {
          for (int idxi = evec_i; idxi < (MIN(N_low, i + N_ii) / 2); idxi++) {
            _lowBuf_i[idxi - evec_i] = evecs->at(idxi);
            (*swapEvecCheckerFn)(idxi);
          }
          if (_lowBuf_i[0].Checkerboard() == Even) {
            l_temp_e = &_lowBuf_i[0];
            l_temp_o = &evecs->at(evec_i);
          } else {
            l_temp_o = &_lowBuf_i[0];
            l_temp_e = &evecs->at(evec_i);
          }
        } else {
          l_temp_e = &evecs->at(evec_i);
          l_temp_o = nullptr;
        }
      } else {
        N_ii = MIN(N_i - i, _blockSize);
        if (skip_low_left)
          l_temp_e = &left[i];
        else
          l_temp_e = &left[i - N_low];

        l_temp_o = nullptr;
      }

      // close this block's stage interval (L/R staging done); the mBlock
      // view construction and the block log line below are sub-us gaps
      STOP_TIMER("untimed: stage");
      stageUs += GET_TIMER("untimed: stage");

      // A2AMatrixSet<T> mBlock(mCache_p, _next, _nstr, _nt, N_ii, N_jj);
      A2AMatrixSet<T> mBlock(mCache.data(), _next, _nstr, _nt, N_ii, N_jj);

      // The former ~1 GB host zeroing of mCache here (thread_for over
      // mBlock.size()) was dead code: every kernel entry path fully
      // overwrites the block before execute reads it again in the
      // conversion pass --
      //   - stencil Full: finishMatOutput copies the whole block D2H and
      //     global-sums it (Grid/GridMilc/a2a/A2AWorkerStencil.h:147-167);
      //   - stencil CB: assembleMat's fused stores write every slot of
      //     the (interleaved) output for all three ContractTypes
      //     (pre-fusion: reconstructLegacy wrote every element,
      //     MesonField.hpp BothHalf/LeftHalf/RightHalf tables);
      //   - legacy/ASlash kernels: StagMesonField zeroes the DEVICE mat
      //     then copies the full block back (Grid/GridMilc/a2a/A2AWorker.h).
      // Proven bit-identical by the gt/devfix whole-tree comparison
      // (test/compare_a2a_hostfix.py) across multi-block schedules --
      // stale-data reuse from a partially-overwritten block would show
      // as a bitwise deviation. a197df6 precedent (dead zero-init drop).
      // The zero timer bracket is kept (empty region -> ~0 in the
      // attribution line) so before/after benchmark logs stay
      // parser-comparable.
      zeroUs = -GET_TIMER("untimed: zero");
      START_TIMER("untimed: zero");
      STOP_TIMER("untimed: zero");
      zeroUs += GET_TIMER("untimed: zero");

      LOG(Message) << "All-to-all matrix block "
                   << i / _blockSize + NBlock_i * j / _blockSize + 1 << "/"
                   << NBlock_i * NBlock_j << " [" << i << " .. " << i + N_ii - 1
                   << ", " << j << " .. " << j + N_jj - 1 << "]" << std::endl;

      flops = 0.0;
      bytes = 0.0;

      wrapUs = -GET_TIMER("kernel");
      START_TIMER("kernel");
      kernel(mBlock, l_temp_e, l_temp_o, r_temp_e, r_temp_o,
             A2ASideStamp{_stampEpoch, i}, A2ASideStamp{_stampEpoch, j});
      STOP_TIMER("kernel");
      wrapUs += GET_TIMER("kernel");

      flops += kernel.flops(N_ii, N_jj, (low_j ? Ncb : 1) * (low_i ? Ncb : 1));
      bytes += kernel.bytes(N_ii, N_jj);

      t_kernel = kernel.kernelTime();
      t_gsum = kernel.globalSumTime();
      // kernel() residue: the "kernel" bracket minus the worker's timed
      // execute + global sum (device D2H, setL/R re-pack, host CB work).
      // GRID_TRACE (ZeroInit/SetLeft/SetRight/CopyFromDevice/GlobalSum)
      // decomposes the Grid-side parts of this residue.
      wrapUs -= t_kernel + t_gsum;

      normUs = -GET_TIMER("untimed: normalise");
      START_TIMER("untimed: normalise");

      // Build this rank's I/O task list BEFORE the conversion pass: the
      // conversion below only touches (e,s) slabs this rank owns (files
      // f = myRank, myRank + nRank, ... of the _next*_nstr pairs), so
      // file-less ranks convert nothing and skip mBuf entirely. The
      // filename/ioname/metadata builders are pure (no side effects, no
      // I/O), so hoisting is behavior-neutral; the list contents are
      // identical to the former in-IO-section build.
#ifdef HADRONS_A2AM_PARALLEL_IO
      unsigned int myRank = _grid->ThisRank(), nRank = _grid->RankCount();
      _nodeIo.clear();
      for (int f = myRank; f < _next * _nstr; f += nRank) {
        IoHelper h;

        h.i = i;
        h.j = j;
        h.e = f / _nstr;
        h.s = f % _nstr;
        h.io = A2AMatrixIoMILC<TIo>(filenameFn(h.e, h.s), ionameFn(h.e, h.s),
                                    _nt, N_i, N_j);
        h.md = metadataFn(h.e, h.s);
        _nodeIo.push_back(h);
      }
#endif

      if (!_nodeIo.empty()) {
        mBuf.resize(mBlock.size());

        // Per-ket-column coefficient table: coeff is a pure function of
        // jj and block-level scalars (low_i/low_j, norm, evec_j, evals)
        // -- character-identical arithmetic to the former per-jj
        // computation inside the conversion loop, evaluated once per
        // block (bit-identical values).
        ComplexD *evals_p = evals.empty() ? nullptr : (ComplexD *)&evals[0];
        std::vector<T> coeffTab(N_jj);
        for (int jj = 0; jj < N_jj; ++jj) {
          T coeff = 1.0;
          if (low_i || low_j) {
            coeff = norm; // Normalize low modes appropriately
            if (low_i && low_j) {
              coeff *= coeff;
            }

            if (low_j) {
              coeff = coeff / evals_p[evec_j + (jj / 2)];
              if (jj % 2 == 1)
                coeff = conjugate(coeff);
            }
          }

          coeffTab[jj] = coeff;
        }

        // Convert + precision-narrow ONLY this rank's I/O-owned (e,s)
        // slabs of the block into mBuf (TIo). The former loop
        // parallelised over jj with the innermost loop striding t at
        // N_ii*N_jj elements (4 MB at production blocks) -- every access
        // a cache/TLB miss, ~1 GB read + ~0.5 GB write per block per
        // rank on replicated data. jj is the RowMajor unit-stride axis
        // (idx = jj + N_jj*(ii + N_ii*t) + slabBase), so this streaming
        // form -- parallel over the slab's (t,ii) rows, jj contiguous
        // innermost -- is bandwidth-bound and computes the SAME
        // per-element multiply (coeffTab[jj] * cache_p[...]: one
        // ComplexD multiply + one narrowing store) exactly once per
        // element, by the slab's owner: bit-identical files (each file
        // receives only its owned slab), at a fraction of the collective
        // traffic.
        const int Lt = _nt;
        const T *cache_p = mBlock.data();
        TIo *result_p = mBuf.data();
        for (const auto &h : _nodeIo) {
          // one contiguous (e,s) slab: the same hyperslab saveBlock
          // writes below (offset (e*nstr + s)*Lt*N_ii*N_jj)
          const size_t slabBase =
              (size_t)(h.e * _nstr + h.s) * Lt * N_ii * N_jj;
          thread_for(ti, Lt * N_ii, {
            const size_t base = slabBase + (size_t)N_jj * ti;
            for (int jj = 0; jj < N_jj; ++jj) {
              result_p[base + jj] = coeffTab[jj] * cache_p[base + jj];
            }
          });
        }
      } // else: file-less rank -- nothing to convert, no mBuf allocation

      // View for the IO section: pointer+dims bookkeeping only (Eigen
      // TensorMap construction performs no data access); on file-less
      // ranks mBuf is empty and the view is never dereferenced (the
      // _nodeIo loop in the IO section does not run).
      A2AMatrixSet<TIo> mIOBlock(mBuf.data(), _next, _nstr, _nt, N_ii, N_jj);

      STOP_TIMER("untimed: normalise");
      normUs += GET_TIMER("untimed: normalise");

      // perf
      LOG(Message) << "Kernel perf " << flops / t_kernel / 1.0e3 / nodes
                   << " Gflop/s/node " << std::endl;
      LOG(Message) << "Kernel Time: " << t_kernel << " us." << std::endl;
      LOG(Message) << "Global Sum Time: " << t_gsum << " us." << std::endl;

      // IO
      double blockSize, ioTime;

      LOG(Message) << "Writing block to disk" << std::endl;
      ioTime = -GET_TIMER("IO: write block");
      barrierUs = -GET_TIMER("untimed: barrier");
      START_TIMER("IO: total");
      makeFileDir(filenameFn(0, 0), _grid);

#ifdef HADRONS_A2AM_PARALLEL_IO
      START_TIMER("untimed: barrier");
      _grid->Barrier();
      STOP_TIMER("untimed: barrier");
      // task list (_nodeIo) was built before the conversion pass above
      // parallel IO
      for (auto &h : _nodeIo) {
        saveBlock(mIOBlock, h);
      }
      START_TIMER("untimed: barrier");
      _grid->Barrier();
      STOP_TIMER("untimed: barrier");
#else
      assert(0);
#endif
      STOP_TIMER("IO: total");
      blockSize =
          static_cast<double>(_next * _nstr * _nt * N_ii * N_jj * sizeof(TIo));
      ioTime += GET_TIMER("IO: write block");
      barrierUs += GET_TIMER("untimed: barrier");
      LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
                   << ioTime << " us ("
                   << blockSize / ioTime * 1.0e6 / 1024 / 1024 << " MB/s)"
                   << std::endl;

      // Additive attribution line; the five contract lines above are
      // byte-identical to the pre-instrumentation build. Values in us.
      LOG(Message) << "A2A block untimed: zero " << zeroUs << " stage "
                   << stageUs << " wrap " << wrapUs << " normalise "
                   << normUs << " barrier " << barrierUs << std::endl;

      // reopen the stage interval: evec restore + loop advance + the
      // next block's L/R staging are the NEXT block's stage
      stageUs = -GET_TIMER("untimed: stage");
      START_TIMER("untimed: stage");

      if (swapEvecCheckerFn != nullptr && low_i &&
          (skip_low_right || i != j)) {
        for (int idxi = evec_i; idxi < (MIN(N_low, i + N_ii) / 2); idxi++)
          evecs->at(idxi) = _lowBuf_i[idxi - evec_i];
      }
      i += N_ii;
      evec_i += (N_ii / Ncb);
    } // End while (i < N_i) Loop

    if (swapEvecCheckerFn != nullptr && low_j) {
      for (int idxj = evec_j; idxj < (MIN(N_low, j + N_jj) / 2); idxj++)
        evecs->at(idxj) = _lowBuf_j[idxj - evec_j];
    }
    j += N_jj;
    evec_j += (N_jj / Ncb);
  } // End while (j < N_j) Loop
  // close the trailing open stage interval (its partial value is
  // discarded; all timers must be stopped for the module summary table)
  STOP_TIMER("untimed: stage");
  // acceleratorFreeDevice(mCache_p);
  MemoryManager::PrintBytes();
}

// I/O handler /////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputationMILC<T, Field, MetadataType, TIo>::saveBlock(
    const A2AMatrixSet<TIo> &m, IoHelper &h) {
  if ((h.i == 0) and (h.j == 0)) {
    START_TIMER("IO: file creation");
    h.io.initFile(h.md, _min_i, _min_j);
    STOP_TIMER("IO: file creation");
  }
  START_TIMER("IO: write block");
  h.io.saveBlock(m, h.e, h.s, h.i, h.j);
  STOP_TIMER("IO: write block");
}

#undef START_TIMER
#undef STOP_TIMER
#undef GET_TIMER

END_HADRONS_NAMESPACE

#endif // A2A_Matrix__hpp_
