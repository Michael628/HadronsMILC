/*
 * LoadMilc.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2024
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Michael Lynch <michaellynch628@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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
#ifndef HadronsMILC_MIO_LoadMilc_hpp_
#define HadronsMILC_MIO_LoadMilc_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

#include <cstring>
#include <fstream>
#include <string>
#include <vector>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Load a MILC v5 gauge configuration

 file                    Namestem of the MILC gauge file to read. The current
                         trajectory number is appended as "<file>.<traj>".
 exitOnChecksumMismatch  If true, a MILC sum29/sum31 checksum mismatch is a
                         fatal error; if false (default) it is only logged.
                         Mirrors MILC's own non-fatal checksum warning and
                         NerscIO::exitOnReadPlaquetteMismatch().

 The plain MILC v5 format (magic 0x4e87) is single-precision, full 3x3 on
 disk, in natural (x-fastest) lexicographic order. The on-disk fsu3_matrix[4]
 is byte-identical to Grid's LorentzColourMatrixF, so the read is a plain
 BINARYIO_LEXICOGRAPHIC load (no checkerboard remap) followed by an
 element-wise float->double munger. The MILC sum29/sum31 integrity checksum
 (word-wise rotate+XOR, distinct from Nersc/SciDAC) is reimplemented here and
 verified against the stored values.
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MIO)

//----------------------------------------------------------------------------
// MILC v5 on-disk constants (see milc_qcd/include/file_types.h,
// io_lat.h and generic/io_lat{4,_utils}.c)
//----------------------------------------------------------------------------
static constexpr uint32_t MILC_V5_MAGIC        = 0x4e87u; // GAUGE_VERSION_NUMBER
static constexpr int      MILC_NATURAL_ORDER   = 0;       // serial/natural site order
static constexpr int      MILC_HEADER_BYTES    = 88;      // magic+dims+timestamp+order
static constexpr int      MILC_CHECKSUM_BYTES  = 8;       // sum29 (4) + sum31 (4)
static constexpr int      MILC_PAYLOAD_OFFSET  = MILC_HEADER_BYTES + MILC_CHECKSUM_BYTES; // 96

//----------------------------------------------------------------------------
// Container for the parsed MILC v5 binary header + stored checksum block.
// The on-disk serialization order (authoritative, from swrite_gauge_hdr in
// io_lat_utils.c) is: magic(4), dims[4](16), time_stamp[64], order(4) = 88 B,
// followed by the checksum block sum29(4)@88, sum31(4)@92 (note: written in
// the OPPOSITE order to the in-memory gauge_check struct field order).
//----------------------------------------------------------------------------
struct MilcHeader
{
    uint32_t         magic{0};
    std::vector<int> dims{std::vector<int>(4, 0)};
    std::string      time_stamp;
    int              order{0};
    uint32_t         sum29{0};
    uint32_t         sum31{0};
    bool             byteReverse{false}; // MILC byterevflag: host must swap bytes
};

//----------------------------------------------------------------------------
// Self-contained, module-local MILC v5 reader. Reads the plain single-
// precision full-3x3 format into a double-precision GaugeField, detects
// endianness from the magic number, reimplements the MILC sum29/sum31
// checksum, and computes the plaquette. Modeled after Grid's OpenQcdIO.
//----------------------------------------------------------------------------
class MilcIO : public BinaryIO {
public:
    typedef Lattice<vLorentzColourMatrixD> GaugeField;

    //----- helper: is the host big-endian? ---------------------------------
    static inline bool hostIsBigEndian(void)
    {
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
        return true;
#elif defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
        return false;
#else
        uint32_t i = 0x01020304u;
        char     c;
        std::memcpy(&c, &i, sizeof(char));
        return (c == 0x01);
#endif
    }

    //----- helper: Grid format string for the detected file endianness -----
    // byteReverse is MILC's byterevflag (host-relative). The Grid format
    // string is the file's *absolute* endianness (Grid's be/le converters
    // are themselves host-aware), so we undo the host-relativity here.
    static inline std::string milcFileFormat(const bool byteReverse)
    {
        const bool fileBigEndian = byteReverse ? !hostIsBigEndian()
                                               : hostIsBigEndian();
        return fileBigEndian ? std::string("IEEE32BIG") : std::string("IEEE32");
    }

    //----- read + parse the 96-byte header/checksum prefix -----------------
    // Runs on every rank (96 bytes only); validates magic, order and the
    // lattice dimensions against the grid. Returns the payload byte offset.
    static inline int readHeader(std::string file, GridBase *grid,
                                 MilcHeader &header)
    {
        std::ifstream fin(file, std::ios::in | std::ios::binary);
        if (!fin.is_open())
        {
            HADRONS_ERROR(Io, "cannot open MILC gauge file '" + file + "'");
        }

        //--- magic number (4 bytes) + endianness detection ---
        // MILC read_gauge_hdr: raw magic == 0x4e87 => host order; its
        // byte-reversal == 0x4e87 => byte-swapped file.
        uint32_t rawMagic = 0;
        fin.read(reinterpret_cast<char *>(&rawMagic), sizeof(uint32_t));
        GRID_ASSERT(!fin.fail());
        bool byteReverse;
        if (rawMagic == MILC_V5_MAGIC)
        {
            byteReverse = false;
        }
        else if (byte_reverse32(rawMagic) == MILC_V5_MAGIC)
        {
            byteReverse = true;
        }
        else
        {
            HADRONS_ERROR(Definition,
                          "unrecognized magic number in MILC gauge file '"
                          + file + "' (expected 0x4e87); not a plain MILC v5 "
                          "configuration (archive/SciDAC/old formats are "
                          "unsupported)");
        }
        header.magic       = MILC_V5_MAGIC;
        header.byteReverse = byteReverse;

        //--- dims[4] (16 B) ---
        uint32_t dims[4];
        fin.read(reinterpret_cast<char *>(dims), sizeof(dims));
        GRID_ASSERT(!fin.fail());

        //--- time_stamp[64] ---
        char ts[64];
        fin.read(ts, sizeof(ts));
        GRID_ASSERT(!fin.fail());

        //--- order (4 B) ---
        uint32_t order = 0;
        fin.read(reinterpret_cast<char *>(&order), sizeof(order));
        GRID_ASSERT(!fin.fail());

        // byte-swap the integer header fields to host order if needed
        if (byteReverse)
        {
            for (int d = 0; d < 4; d++)
                dims[d] = byte_reverse32(dims[d]);
            order = byte_reverse32(order);
        }
        header.dims.assign(dims, dims + 4);
        header.time_stamp = std::string(ts, sizeof(ts));
        {
            const auto nul = header.time_stamp.find('\0');
            if (nul != std::string::npos)
                header.time_stamp.resize(nul);
        }
        header.order = static_cast<int>(order);

        //--- checksum block: sum29@88, sum31@92 (note: sum29 first on disk) ---
        uint32_t sum29 = 0, sum31 = 0;
        fin.read(reinterpret_cast<char *>(&sum29), sizeof(sum29));
        fin.read(reinterpret_cast<char *>(&sum31), sizeof(sum31));
        GRID_ASSERT(!fin.fail());
        if (byteReverse)
        {
            sum29 = byte_reverse32(sum29);
            sum31 = byte_reverse32(sum31);
        }
        header.sum29 = sum29;
        header.sum31 = sum31;

        const int data_start = static_cast<int>(fin.tellg()); // 96
        fin.close();

        //--- validate site ordering (coordinate-list/checkpoint unsupported) ---
        if (header.order != MILC_NATURAL_ORDER)
        {
            HADRONS_ERROR(Definition,
                          "unsupported site order "
                          + std::to_string(header.order) + " in MILC gauge file '"
                          + file + "' (only natural order / 0 is supported)");
        }

        //--- validate lattice dimensions against the grid ---
        GRID_ASSERT(grid->_ndimension == Nd);
        for (int d = 0; d < Nd; d++)
            GRID_ASSERT(static_cast<int>(grid->_fdimensions[d]) == header.dims[d]);

        return data_start;
    }

    //----- MILC sum29/sum31 checksum (local contribution) -----------------
    // Mirrors BinaryIO::ScidacChecksum structure but uses MILC's word-wise
    // rotate+XOR with a GLOBAL word index (rank = globalWord % {29,31}).
    // Computed over the endian-corrected (host-order) payload, matching
    // MILC which byte-reverses before accumulating. The rank==0 shift guard
    // avoids the 32-bit undefined-behaviour shift. Returns LOCAL per-rank
    // contributions; the caller reduces with GlobalXOR (as MILC's g_xor32
    // does, and as BinaryIO::IOobject does for SciDAC).
    template <class fobj>
    static inline void milcChecksum(GridBase *grid, std::vector<fobj> &fbuf,
                                    uint32_t &milc_csum29, uint32_t &milc_csum31)
    {
        const uint64_t size32 = sizeof(fobj) / sizeof(uint32_t); // 72 for fsu3_matrix[4]
        const int      nd     = grid->_ndimension;
        const uint64_t lsites = grid->lSites();

        Coordinate local_vol   = grid->LocalDimensions();
        Coordinate local_start = grid->LocalStarts();
        Coordinate global_vol  = grid->FullDimensions();

        milc_csum29 = 0;
        milc_csum31 = 0;

        thread_region
        {
            Coordinate coor(nd);
            uint32_t   csum29_thr = 0, csum31_thr = 0;

            thread_for_in_region(local_site, lsites,
            {
                // local lexicographic site -> global site (x-fastest, == MILC natural order)
                Lexicographic::CoorFromIndex(coor, local_site, local_vol);
                for (int d = 0; d < nd; d++)
                    coor[d] += local_start[d];
                int64_t global_site;
                Lexicographic::IndexFromCoor(coor, global_site, global_vol);

                const uint32_t *site_buf
                    = reinterpret_cast<const uint32_t *>(&fbuf[local_site]);
                const uint64_t  base
                    = static_cast<uint64_t>(global_site) * size32;

                for (uint64_t k = 0; k < size32; k++)
                {
                    const uint64_t  idx    = base + k;
                    const uint32_t  v      = site_buf[k];
                    const uint32_t  rank29 = idx % 29u;
                    const uint32_t  rank31 = idx % 31u;

                    csum29_thr ^= (v << rank29)
                                | (rank29 == 0u ? 0u : v >> (32u - rank29));
                    csum31_thr ^= (v << rank31)
                                | (rank31 == 0u ? 0u : v >> (32u - rank31));
                }
            });

            thread_critical
            {
                milc_csum29 ^= csum29_thr;
                milc_csum31 ^= csum31_thr;
            }
        }
    }

    //----- read a MILC v5 configuration into a double GaugeField ----------
    template <class GaugeStats = PeriodicGaugeStatistics>
    static inline void readConfiguration(GaugeField &Umu, MilcHeader &header,
                                         std::string file,
                                         const bool  exitOnMismatch = false,
                                         GaugeStats  GaugeStatisticsCalculator
                                             = GaugeStats())
    {
        GridBase *grid = Umu.Grid();

        uint64_t       offset = readHeader(file, grid, header);
        const std::string format = milcFileFormat(header.byteReverse);

        LOG(Message) << "MILC configuration '" << file << "' "
                     << (header.byteReverse ? "byte-swapped" : "host-order")
                     << " (" << format << "), time stamp '"
                     << header.time_stamp << "', site order " << header.order
                     << ", lattice " << header.dims[0] << "x" << header.dims[1]
                     << "x" << header.dims[2] << "x" << header.dims[3]
                     << std::endl;

        typedef LorentzColourMatrixF fobj; // on-disk single-precision fsu3_matrix[4]
        typedef typename GaugeField::vector_object::scalar_object sobj; // double

        const uint64_t       lsites = grid->lSites();
        std::vector<fobj>    iodata(lsites);     // endian-corrected on disk objects
        std::vector<sobj>    scalardata(lsites); // munged double objects

        // parallel MPI-IO read; IOobject swaps iodata to host order in place
        // (and computes Grid's own Nersc/SciDAC checksums, unused here).
        float    wf = 0.f;
        uint32_t nersc_csum, scidac_csuma, scidac_csumb;
        IOobject(wf, grid, iodata, file, offset, format,
                 BINARYIO_READ | BINARYIO_LEXICOGRAPHIC,
                 nersc_csum, scidac_csuma, scidac_csumb);

        // MILC sum29/sum31 over the endian-corrected payload.
        uint32_t milc_csum29, milc_csum31;
        milcChecksum(grid, iodata, milc_csum29, milc_csum31);
        grid->GlobalXOR(milc_csum29);
        grid->GlobalXOR(milc_csum31);

        const bool match29 = (milc_csum29 == header.sum29);
        const bool match31 = (milc_csum31 == header.sum31);

        if (grid->IsBoss())
        {
            LOG(Message) << "MILC checksum '" << file << "': "
                         << "sum29 computed 0x" << std::hex << milc_csum29
                         << " stored 0x" << header.sum29
                         << (match29 ? " OK" : " MISMATCH")
                         << " | sum31 computed 0x" << milc_csum31
                         << " stored 0x" << header.sum31
                         << (match31 ? " OK" : " MISMATCH")
                         << std::dec << std::endl;
        }
        if (!match29 || !match31)
        {
            if (exitOnMismatch)
            {
                HADRONS_ERROR(Io, "MILC gauge configuration checksum mismatch "
                                   "in '" + file + "'");
            }
            else if (grid->IsBoss())
            {
                LOG(Message) << "MILC checksum mismatch in '" << file
                             << "' (continuing; set exitOnChecksumMismatch=true "
                             "to make fatal)" << std::endl;
            }
        }

        // munge float -> double (element-wise, no transpose/conjugate) + vectorize
        GaugeSimpleMunger<fobj, sobj> munge;
        thread_for(x, lsites, { munge(iodata[x], scalardata[x]); });
        vectorizeFromLexOrdArray(scalardata, Umu);
        grid->Barrier();

        // plaquette / link trace (informational; MILC binary header carries none)
        FieldMetaData stats;
        GaugeStatisticsCalculator(Umu, stats);
        LOG(Message) << "MILC configuration '" << file << "' plaquette "
                     << stats.plaquette << " link_trace " << stats.link_trace
                     << std::endl;
    }
};

/******************************************************************************
 * LoadMilc module parameters
 ******************************************************************************/
class LoadMilcPar : Serializable {
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadMilcPar,
                                    std::string, file,
                                    bool,        exitOnChecksumMismatch);
    // warn by default on a sum29/sum31 mismatch (mirrors MILC's non-fatal
    // warning); opt-in to a fatal error via exitOnChecksumMismatch=true.
    LoadMilcPar(void)
        : exitOnChecksumMismatch(false) {}
};

/******************************************************************************
 * TLoadMilc module
 ******************************************************************************/
template <typename GImpl> class TLoadMilc : public Module<LoadMilcPar> {
public:
    GAUGE_TYPE_ALIASES(GImpl, );

public:
    // constructor
    TLoadMilc(const std::string name);
    // destructor
    virtual ~TLoadMilc(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadMilc, TLoadMilc<GIMPL>, MIO);

/******************************************************************************
 *                       TLoadMilc implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLoadMilc<GImpl>::TLoadMilc(const std::string name)
    : Module<LoadMilcPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TLoadMilc<GImpl>::getInput(void) {
    return {};
}

template <typename GImpl>
std::vector<std::string> TLoadMilc<GImpl>::getOutput(void) {
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl> void TLoadMilc<GImpl>::setup(void) {
    envCreateLat(GaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl> void TLoadMilc<GImpl>::execute(void) {
    MilcHeader  header;
    std::string fileName
        = par().file + "." + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading MILC configuration from file '" << fileName << "'"
                 << std::endl;

    auto &U = envGet(GaugeField, getName());
    MilcIO::readConfiguration(U, header, fileName, par().exitOnChecksumMismatch);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MIO_LoadMilc_hpp_
