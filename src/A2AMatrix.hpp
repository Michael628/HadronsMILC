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
  virtual double flops(const unsigned int blockSizei,
                       const unsigned int blockSizej, int cbDiv = 1) = 0;
  virtual double bytes(const unsigned int blockSizei,
                       const unsigned int blockSizej) = 0;
  virtual double kernelTime() = 0;
  virtual double globalSumTime() = 0;
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
               const SwapFn *swapEvecCheckerFn = nullptr);

private:
  // I/O handler
  void saveBlock(const A2AMatrixSet<TIo> &m, IoHelper &h);

private:
  TimerArray *_tArray;
  GridBase *_grid;
  unsigned int _orthogDim, _nt, _next, _nstr, _blockSize, _min_i, _min_j;
  std::vector<IoHelper> _nodeIo;
  std::vector<Field> _lowBuf_i, _lowBuf_j;
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
    const SwapFn *swapEvecCheckerFn) {
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over _blockSize factors
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////

  std::vector<T> mCache;
  mCache.resize(_nt * _next * _nstr * _blockSize * _blockSize);
  // LOG(Message) << "cache bytes: " << mCache_bytes << std::endl;

  MemoryManager::PrintBytes();

  std::vector<TIo> mBuf;

  bool checkerboarded_low = (swapEvecCheckerFn != nullptr);
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

    if (!skip_low_left)
      _lowBuf_i.resize(
          _blockSize / 2,
          evecs->at(0).Grid()); // storage for caching checkerboards
    if (!skip_low_right)
      _lowBuf_j.resize(_blockSize / 2, evecs->at(0).Grid());
  }

  int NBlock_i = N_i / _blockSize +
                 (((N_i % _blockSize) != 0)
                      ? 1
                      : 0); // Round up on the number of blocks to compute
  int NBlock_j = N_j / _blockSize + (((N_j % _blockSize) != 0) ? 1 : 0);

  bool low_i, low_j;
  int i, j, evec_i, evec_j, N_ii, N_jj;

  j = 0, evec_j = 0;
  while (j < N_j) { // While we still have bra vectors to contract

    low_j = j < N_low && !skip_low_right;

    const Field *r_temp_e, *r_temp_o;
    if (low_j) {
      N_jj = MIN(N_low - j, _blockSize);

      if (checkerboarded_low) {
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

      // A2AMatrixSet<T> mBlock(mCache_p, _next, _nstr, _nt, N_ii, N_jj);
      A2AMatrixSet<T> mBlock(mCache.data(), _next, _nstr, _nt, N_ii, N_jj);

      // accelerator_for(idx,mBlock.size(),1,{
      thread_for(idx, mBlock.size(), { mCache[idx] = 0.0; });

      LOG(Message) << "All-to-all matrix block "
                   << i / _blockSize + NBlock_i * j / _blockSize + 1 << "/"
                   << NBlock_i * NBlock_j << " [" << i << " .. " << i + N_ii - 1
                   << ", " << j << " .. " << j + N_jj - 1 << "]" << std::endl;

      flops = 0.0;
      bytes = 0.0;

      START_TIMER("kernel");
      kernel(mBlock, l_temp_e, l_temp_o, r_temp_e, r_temp_o);
      STOP_TIMER("kernel");

      flops += kernel.flops(N_ii, N_jj, (low_j ? Ncb : 1) * (low_i ? Ncb : 1));
      bytes += kernel.bytes(N_ii, N_jj);

      t_kernel = kernel.kernelTime();
      t_gsum = kernel.globalSumTime();

      mBuf.resize(mBlock.size());
      A2AMatrixSet<TIo> mIOBlock(mBuf.data(), _next, _nstr, _nt, N_ii, N_jj);

      {
        int next = _next, nstr = _nstr, Lt = _nt;
        ComplexD *evals_p = (ComplexD *)&evals[0];
        TIo *result_p = mBuf.data();
        T *cache_p = mBlock.data();
        // accelerator_for(jj,N_jj,1,{
        thread_for(jj, N_jj, {
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

          for (int ii = 0; ii < N_ii; ii++)
            for (int e = 0; e < next; e++)
              for (int s = 0; s < nstr; s++)
                for (int t = 0; t < Lt; t++) {
                  int idx = jj + N_jj * (ii + N_ii * (t + Lt * (s + nstr * e)));
                  // If the ket vectors (corresponding to the solves) are low
                  // modes, multiply by the eigenvals
                  result_p[idx] = coeff * cache_p[idx];
                }
        });
      }

      // perf
      LOG(Message) << "Kernel perf " << flops / t_kernel / 1.0e3 / nodes
                   << " Gflop/s/node " << std::endl;
      LOG(Message) << "Kernel Time: " << t_kernel << " us." << std::endl;
      LOG(Message) << "Global Sum Time: " << t_gsum << " us." << std::endl;

      // IO
      double blockSize, ioTime;
      unsigned int myRank = _grid->ThisRank(), nRank = _grid->RankCount();

      LOG(Message) << "Writing block to disk" << std::endl;
      ioTime = -GET_TIMER("IO: write block");
      START_TIMER("IO: total");
      makeFileDir(filenameFn(0, 0), _grid);

#ifdef HADRONS_A2AM_PARALLEL_IO
      _grid->Barrier();
      // make task list for current node
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
      // parallel IO
      for (auto &h : _nodeIo) {
        saveBlock(mIOBlock, h);
      }
      _grid->Barrier();
#else
      assert(0);
#endif
      STOP_TIMER("IO: total");
      blockSize =
          static_cast<double>(_next * _nstr * _nt * N_ii * N_jj * sizeof(TIo));
      ioTime += GET_TIMER("IO: write block");
      LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
                   << ioTime << " us ("
                   << blockSize / ioTime * 1.0e6 / 1024 / 1024 << " MB/s)"
                   << std::endl;

      if (checkerboarded_low && low_i && (skip_low_right || i != j)) {
        for (int idxi = evec_i; idxi < (MIN(N_low, i + N_ii) / 2); idxi++)
          evecs->at(idxi) = _lowBuf_i[idxi - evec_i];
      }
      i += N_ii;
      evec_i += (N_ii / Ncb);
    } // End while (i < N_i) Loop

    if (checkerboarded_low && low_j) {
      for (int idxj = evec_j; idxj < (MIN(N_low, j + N_jj) / 2); idxj++)
        evecs->at(idxj) = _lowBuf_j[idxj - evec_j];
    }
    j += N_jj;
    evec_j += (N_jj / Ncb);
  } // End while (j < N_j) Loop
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
