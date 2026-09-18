/*
 * LoadMesonField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
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
#ifndef HadronsMILC_MIO_LoadMesonField_hpp_
#define HadronsMILC_MIO_LoadMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <A2AMatrix.hpp>
#include <EigenPack.hpp>
#include <Modules/MContraction/MesonField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Load a meson-field HDF5 file                           *
 ******************************************************************************/
/*  Load a meson field [nt, N_i, N_j] previously written by TMesonFieldMILC
    (one HDF5 file per momentum x spin-taste at
    "<output>.<traj>/<dataname>.h5", dataset "a2aMatrix" inside a group named
    by the dataname, metadata group serialized from MesonFieldMILCMetadata).

    file        Name of the HDF5 file. The token "@traj@" (if present) is
                replaced by the current trajectory number, e.g.
                "mfout.@traj@/G1_G1_0_0_0.h5".
    dataset     Name of the HDF5 group holding the "a2aMatrix" dataset, i.e.
                the dataname used by the writer (spinGammaName_tasteGammaName_
                px_py_pz, e.g. "G1_G1_0_0_0").
    side        Optional eigenvalue un-weighting mode:
                  ""    (default) load as-is
                  "bra" rows are bra-side eigenvector rows: they carry a
                        uniform normalization only, nothing to un-weight
                  "ket" columns are ket-side eigenvector columns carrying a
                        norm/eval_k factor (conjugated on odd columns, see
                        A2AMatrix.hpp coeff folding); multiply the eigenvalues
                        back in, inverting the writer's folding
    lowModes    Name of the eigenpack module (MassShiftEigenPack) providing
                the eval_k = 2m + i*lambda_D values; required when
                side=ket (the eigenvalues are not stored in the HDF5 file)

    The loader refuses files that are not identity spin-taste (G1,G1) at zero
    momentum: only for that kernel are the stored elements plain per-timeslice
    overlaps <left_i|right_j>(t), which downstream consumers (the
    StagLMAMesonField solver) require.
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadMesonFieldPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LoadMesonFieldPar, std::string, file,
                                  std::string, dataset, std::string, side,
                                  std::string, lowModes);
  LoadMesonFieldPar(void) : side(""), lowModes("") {}
};

template <typename FImpl, typename Pack>
class TLoadMesonField : public Module<LoadMesonFieldPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  // constructor
  TLoadMesonField(const std::string name);
  // destructor
  virtual ~TLoadMesonField(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadMesonField,
                    ARG(TLoadMesonField<STAGIMPL, MassShiftEigenPack<STAGIMPL>>),
                    MIO);

/******************************************************************************
 *                       TLoadMesonField implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TLoadMesonField<FImpl, Pack>::TLoadMesonField(const std::string name)
    : Module<LoadMesonFieldPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TLoadMesonField<FImpl, Pack>::getInput(void) {
  std::vector<std::string> in = {};
  if (!par().lowModes.empty()) {
    in.push_back(par().lowModes);
  }

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLoadMesonField<FImpl, Pack>::getOutput(void) {
  std::vector<std::string> out = {getName()};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TLoadMesonField<FImpl, Pack>::setup(void) {
  int nt = env().getDim().back();

  if (par().file.empty() || par().dataset.empty()) {
    HADRONS_ERROR(Argument,
                  "both 'file' and 'dataset' parameters must be provided");
  }
  if ((par().side != "") && (par().side != "bra") && (par().side != "ket")) {
    HADRONS_ERROR(Argument, "side must be '', 'bra' or 'ket' (got '" +
                                par().side + "')");
  }
  if ((par().side == "ket") && par().lowModes.empty()) {
    HADRONS_ERROR(Argument,
                  "side=ket requires the 'lowModes' eigenpack reference for "
                  "the eigenvalues");
  }

  envCreate(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, getName(), 1, nt);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TLoadMesonField<FImpl, Pack>::execute(void) {
#ifdef HAVE_HDF5
  std::string fileName = par().file;
  int nt = env().getDim().back();

  tokenReplace(fileName, "traj", vm().getTrajectory());
  LOG(Message) << "Loading meson field '" << par().dataset << "' from '"
               << fileName << "'" << std::endl;

  // validate metadata: only identity spin-taste at zero momentum stores plain
  // per-timeslice overlaps, which is what downstream consumers reconstruct
  // from
  {
    MContraction::MesonFieldMILCMetadata md;
    Hdf5Reader reader(fileName);

    // metadata read is unexercised territory: a corrupt/missing metadata
    // group (schema drift) must surface clearly at load time, not silently
    // pass
    try
    {
      read(reader, par().dataset, md);
    }
    catch (const std::exception &e)
    {
      HADRONS_ERROR(Io, "failed reading meson-field metadata from '"
                            + fileName + "/" + par().dataset + "': " + e.what());
    }

    bool zeroMom = true;
    for (auto &pmu : md.momentum) {
      if (pmu != 0) {
        zeroMom = false;
      }
    }
    if (!zeroMom || (md.gamma_spin != StagGamma::StagAlgebra::G1) ||
        (md.gamma_taste != StagGamma::StagAlgebra::G1)) {
      HADRONS_ERROR(Argument,
                    "meson field '" + fileName + "' is not identity "
                    "spin-taste (G1,G1) at zero momentum: its elements carry "
                    "spin-taste/momentum phases that cannot be reconstructed "
                    "from by this module's consumers");
    }
  }

  // load the full [nt, N_i, N_j] table (rank-strided read + broadcast).
  // This read is inlined rather than going through A2AMatrixIoMILC::load():
  // load()'s VecT template parameter only admits scalar element containers
  // (EigenDiskVector-style, its sole historical caller), and cannot be
  // instantiated for the std::vector<A2AMatrix<...>> held in the environment
  // (its body casts each timeslice with .cast<VecT>(), which requires VecT
  // to be a scalar type). The rank-striding, hyperslab reads, broadcast and
  // size validation below mirror load()'s grid path.
  auto &mf =
      envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, getName());
  {
    GridBase *grid = envGetGrid(FermionField);
    unsigned int myRank = grid->ThisRank(), nRank = grid->RankCount();

    Hdf5Reader reader(fileName);
    push(reader, par().dataset);
    auto &group = reader.getGroup();
    H5NS::DataSet dataset;
    H5NS::DataSpace dataspace;
    H5NS::CompType datatype;

    dataset = group.openDataSet(HADRONS_A2AM_NAME);
    datatype = dataset.getCompType();
    dataspace = dataset.getSpace();
    std::vector<hsize_t> hdim(dataspace.getSimpleExtentNdims());
    dataspace.getSimpleExtentDims(hdim.data());

    if (hdim.size() != 3) {
      HADRONS_ERROR(Size, "all-to-all matrix '" + fileName + "/" +
                              par().dataset + "/" + HADRONS_A2AM_NAME +
                              "' is not 3-dimensional");
    }
    if (hdim[0] != static_cast<hsize_t>(nt)) {
      HADRONS_ERROR(Size, "all-to-all time size mismatch (got " +
                              std::to_string(hdim[0]) + ", expected " +
                              std::to_string(nt) + ")");
    }
    const hsize_t ni = hdim[1], nj = hdim[2];

    // size every timeslice on every rank before the broadcast loop
    for (int t = 0; t < nt; ++t) {
      mf[t].resize(ni, nj);
    }

    std::vector<hsize_t> count = {1, ni, nj}, stride = {1, 1, 1},
                         block = {1, 1, 1}, memCount = {ni, nj};
    H5NS::DataSpace memspace(memCount.size(), memCount.data());

    double tRead = 0.;
    int broadcastSize = sizeof(HADRONS_A2AM_IO_TYPE) * mf[0].size();

    std::cout << "Loading timeslice";
    std::cout.flush();
    for (int t = myRank; t < nt; t += nRank) {
      std::vector<hsize_t> offset = {static_cast<hsize_t>(t), 0, 0};

      std::cout << " " << t;
      std::cout.flush();

      dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                stride.data(), block.data());
      tRead -= usecond();
      dataset.read(mf[t].data(), datatype, memspace, dataspace);
      tRead += usecond();
    }
    grid->Barrier();
    for (int t = 0; t < nt; ++t) {
      grid->Broadcast(t % nRank, mf[t].data(), broadcastSize);
    }
    grid->Barrier();
    std::cout << std::endl;

    LOG(Message) << "Read "
                 << nt * ni * nj * sizeof(HADRONS_A2AM_IO_TYPE)
                 << " bytes in " << tRead << " us (" << nt << "x" << ni << "x"
                 << nj << " timeslices x rows x columns)" << std::endl;
  }

  // optional eigenvalue un-weighting of ket-side columns
  if (par().side == "ket") {
    auto &epack = envGet(Pack, par().lowModes);
    unsigned int nEvec = epack.evec.size();
    unsigned int nj = static_cast<unsigned int>(mf[0].cols());

    if (nj == nEvec) {
      HADRONS_ERROR(Implementation,
                    "ket-side un-weighting of full-volume eigenvector files "
                    "is unsupported: the writer's coefficient folding "
                    "mis-indexes eigenvalues for non-checkerboarded ket "
                    "columns (jj/2 indexing in A2AMatrix.hpp)");
    } else if (nj != 2 * nEvec) {
      HADRONS_ERROR(Size, "ket-side un-weighting: file has " +
                              std::to_string(nj) + " columns, expected " +
                              std::to_string(2 * nEvec) + " (paired columns "
                              "for " + std::to_string(nEvec) +
                              " eigenvectors)");
    }

    // invert the writer's folding: stored_{2k} = (norm/eval_k) * raw and
    // stored_{2k+1} = conjugate(norm/eval_k) * raw (A2AMatrix.hpp coeff
    // folding); the uniform norm is deliberately kept (the solver reconciles
    // it), only the eigenvalue factors are multiplied back in
    for (int t = 0; t < nt; ++t) {
      for (unsigned int k = 0; k < nEvec; ++k) {
        auto eval = static_cast<HADRONS_A2AM_IO_TYPE>(epack.eval[k]);

        mf[t].col(2 * k) *= eval;
        mf[t].col(2 * k + 1) *= conjugate(eval);
      }
    }
    LOG(Message) << "Un-weighted eigenvalue factors on " << 2 * nEvec
                 << " ket-side columns" << std::endl;
  } else if (par().side == "bra") {
    LOG(Message) << "bra-side rows carry a uniform normalization only; "
                    "nothing to un-weight"
                 << std::endl;
  }
#else
  HADRONS_ERROR(Implementation, "meson field I/O needs HDF5 library");
#endif
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MIO_LoadMesonField_hpp_
