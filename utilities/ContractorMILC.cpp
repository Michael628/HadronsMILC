/*
 * Contractor.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#include <Hadrons/Global.hpp>
#include <src/HadronsMILC/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/EigenPack.hpp>

using namespace Grid;
using namespace Hadrons;

#define TIME_MOD(t) (((t) + par.global.nt) % par.global.nt)

namespace ContractorMILC
{
    class TrajRange: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step);
    };
    
    class GlobalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        TrajRange, trajCounter,
                                        unsigned int, nt,
                                        std::string, diskVectorDir,
                                        std::string, output);
    };

    class EpackPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(EpackPar,
                                        std::string, fileStem,
                                        bool, evalFile,
                                        bool, multiFile,
                                        int, nEigs,
                                        RealD, massOld,
                                        RealD, massNew);
        EpackPar(void): 
        fileStem{"N/A"}, evalFile{false}, multiFile{false},
        nEigs{0}, massOld{0.0},massNew{0.0} {}
    };

    class EigenValueFile: Serializable {
    public:
      GRID_SERIALIZABLE_CLASS_MEMBERS(EigenValueFile,
         std::vector<RealD>, evals);
    };

    class A2AMatrixPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                        std::string, file,
                                        std::string, dataset,
                                        unsigned int, cacheSize,
                                        unsigned int, ni,
                                        unsigned int, nj,
                                        unsigned int, niOffset,
                                        unsigned int, njOffset,
                                        std::string, name,
                                        EpackPar, epack);
    };

    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string, terms,
                                        std::vector<std::string>, times,
                                        std::string, translations,
                                        bool, translationAverage,
                                        bool, spaceNormalize);
    };

    class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<ContractorMILC::A2AMatrixPar>,  a2aMatrix,
                                        ProductPar, contraction,
                                        std::vector<unsigned int>, times,
                                        std::vector<ComplexD>, correlator);
    };
}

struct ContractorMILCPar
{
    ContractorMILC::GlobalPar                  global;
    std::vector<ContractorMILC::A2AMatrixPar>  a2aMatrix;
    std::vector<ContractorMILC::ProductPar>    product;
};

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times,
                 std::vector<unsigned int> &current,
                 const unsigned int depth)
{
    if (depth > 0)
    {
        for (auto t: times[times.size() - depth])
        {
            current[times.size() - depth] = t;
            makeTimeSeq(timeSeq, times, current, depth - 1);
        }
    }
    else
    {
        timeSeq.push_back(current);
    }
}

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times)
{
    std::vector<unsigned int> current(times.size());

    makeTimeSeq(timeSeq, times, current, times.size());
}

void saveCorrelator(const ContractorMILC::CorrelatorResult &result, const std::string dir, 
                    const unsigned int dt, const unsigned int traj)
{
    std::string              fileStem = "", filename;
    std::vector<std::string> terms = strToVec<std::string>(result.contraction.terms);

    for (unsigned int i = 0; i < terms.size() - 1; i++)
    {
        fileStem += terms[i] + "_" + std::to_string(result.times[i]) + "_";
    }
    fileStem += terms.back();
    if (!result.contraction.translationAverage)
    {
        fileStem += "_dt_" + std::to_string(dt);
    }
    filename = dir + "/" + ModuleBase::resultFilename(fileStem, traj);
    LOG(Message) << "Saving correlator to '" << filename << "'" << std::endl;
    makeFileDir(dir);
    ResultWriter writer(filename);
    write(writer, fileStem, result);
}

std::set<unsigned int> parseTimeRange(const std::string str, const unsigned int nt)
{
    std::regex               rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
    std::smatch              sm;
    std::vector<std::string> rstr = strToVec<std::string>(str);
    std::set<unsigned int>   tSet;

    for (auto &s: rstr)
    {
        std::regex_match(s, sm, rex);
        if (sm[1].matched)
        {
            unsigned int t;
            
            t = std::stoi(sm[1].str());
            if (t >= nt)
            {
                HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
            }
            tSet.insert(t);
        }
        else if (sm[2].matched)
        {
            unsigned int ta, tb;

            ta = std::stoi(sm[3].str());
            tb = std::stoi(sm[4].str());
            if ((ta >= nt) or (tb >= nt))
            {
                HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
            }
            for (unsigned int ti = ta; ti <= tb; ++ti)
            {
                tSet.insert(ti);
            }
        }
    }

    return tSet;
}

struct Sec
{
    Sec(const double usec)
    {
        seconds = usec/1.0e6;
    }
    
    double seconds;
};

inline std::ostream & operator<< (std::ostream& s, const Sec &&sec)
{
    s << std::setw(10) << sec.seconds << " sec";

    return s;
}

struct Flops
{
    Flops(const double flops, const double fusec)
    {
        gFlopsPerSec = flops/fusec/1.0e3;
    }
    
    double gFlopsPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Flops &&f)
{
    s << std::setw(10) << f.gFlopsPerSec << " GFlop/s";

    return s;
}

struct Bytes
{
    Bytes(const double bytes, const double busec)
    {
        gBytesPerSec = bytes/busec*1.0e6/1024/1024/1024;
    }
    
    double gBytesPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Bytes &&b)
{
    s << std::setw(10) << b.gBytesPerSec << " GB/s";

    return s;
}

int main(int argc, char* argv[])
{
    // parse command line
    std::string   parFilename;

    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> <grid flags>";
        std::cerr << std::endl;
        
        return EXIT_FAILURE;
    }
    parFilename = argv[1];

     Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

   // parse parameter file
    ContractorMILCPar par;
    //    unsigned int  nMat,nCont;
    {
    XmlReader     reader(parFilename);

    read(reader, "global",    par.global);
    read(reader, "a2aMatrix", par.a2aMatrix);
    read(reader, "product",   par.product);
    }
    //    nMat  = par.a2aMatrix.size();
    //    nCont = par.product.size();

    // create diskvectors
    std::map<std::string, std::vector<ComplexD> > evalMap;
    std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;
    //    unsigned int                                     cacheSize;

    GridCartesian*    grid     = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                     GridDefaultSimd(Nd,vComplex::Nsimd()),
                                     GridDefaultMpi());

    auto &dims = grid->_gdimensions;
    RealD vol = 1.;
    for (int idx=0;idx < dims.size()-1;idx++) {
        vol *= dims[idx];
    }

    for (auto &p: par.a2aMatrix)
    {
        std::string dirName = par.global.diskVectorDir + "/" + p.name;
        if (!grid->IsBoss()) {
            int rank = grid->ThisRank();
            dirName += "-tmp-";
            dirName += std::to_string(rank);
        }

        a2aMat.emplace(p.name, EigenDiskVector<ComplexD>(dirName, par.global.nt, p.cacheSize));
    }

    // trajectory loop
    for (unsigned int traj = par.global.trajCounter.start; 
         traj < par.global.trajCounter.end; traj += par.global.trajCounter.step)
    {
        LOG(Message) << ":::::::: Trajectory " << traj << std::endl;

        // load data
        for (auto &p: par.a2aMatrix)
        {
            std::string filename = p.file;
            double      t;

            tokenReplace(filename, "traj", traj);
            LOG(Message) << "======== Loading '" << filename << "'" << std::endl;

            A2AMatrixIoMILC<HADRONS_A2AM_IO_TYPE> a2aIo(filename, p.dataset, par.global.nt,p.ni,p.nj,p.niOffset,p.njOffset);

            a2aIo.load(a2aMat.at(p.name), &t,grid);
            LOG(Message) << "Read " << a2aIo.getSize() << " bytes in " << t/1.0e6 
                    << " sec, " << a2aIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;

            if (!p.epack.fileStem.empty()) {

                std::vector<RealD> evals;
                int nEvals;
                RealD massOld = p.epack.massOld;
                RealD massNew = p.epack.massNew;
                if (p.epack.evalFile) {
                    ResultReader reader(ModuleBase::resultFilename(p.epack.fileStem,traj));
                    read(reader,"evals",evals);
                } else {
                    PackRecord record;
                    std::string filename = p.epack.fileStem + "." + std::to_string(traj) + (p.epack.multiFile ? "" : ".bin");
                    EigenPackIo::readEvals(evals,record,0,p.epack.nEigs,filename,
                       p.epack.multiFile);

                }
                nEvals = evals.size();

                evalMap.emplace(p.name,std::vector<ComplexD>(2 * nEvals));

                ComplexD newEval;
                int newIndex = 2*nEvals-p.njOffset-1;
                for (int i = nEvals; i > p.njOffset/2; i--) {
                    if (massOld == massNew) {
                        if (i == nEvals) {
                            LOG(Message) << "Removing  Eigenvalue weight from " << p.name << " with mass=" << massNew << std::endl;
                        }
                        newEval = ComplexD(2.0*p.epack.massOld,sqrt(evals[i-1]));
                    } else {
                        if (i == nEvals) {
                            LOG(Message) << "Shifting Eigenvalues of " << p.name << " from mass=" << massOld << " to mass=" << massNew << std::endl;
                        }
                        newEval = ComplexD(2.0*p.epack.massOld,sqrt(evals[i-1]))/ComplexD(2.0*p.epack.massNew,sqrt(evals[i-1]));
                    }

                    evalMap.at(p.name)[newIndex-1] = newEval;
                    evalMap.at(p.name)[newIndex] = conjugate(newEval);
                    newIndex -= 2;
                }
            }
        }

        if (grid->IsBoss()) {
            for (auto &p: par.product)
            {
                std::vector<std::string>               term = strToVec<std::string>(p.terms);
                std::vector<std::set<unsigned int>>    times;
                std::vector<std::vector<unsigned int>> timeSeq;
                std::set<unsigned int>                 translations;
                std::vector<A2AMatrixTr<ComplexD>>     lastTerm(par.global.nt);
                A2AMatrix<ComplexD>                    prod, tmp, ref;
                TimerArray                             tAr;
                double                                 fusec, busec, flops, bytes;
    	    //	    double  tusec;
                ContractorMILC::CorrelatorResult           result;             

                tAr.startTimer("Total");
                LOG(Message) << "======== Contraction tr(";
                for (unsigned int g = 0; g < term.size(); ++g)
                {
                    std::cout << term[g] << ((g == term.size() - 1) ? ')' : '*');
                }
                std::cout << std::endl;
                if (term.size() != p.times.size() + 1)
                {
                    HADRONS_ERROR(Size, "number of terms (" + std::to_string(term.size()) 
                                + ") different from number of times (" 
                                + std::to_string(p.times.size() + 1) + ")");
                }
                for (auto &s: p.times)
                {
                    times.push_back(parseTimeRange(s, par.global.nt));
                }
                for (auto &m: par.a2aMatrix)
                {
                    if (std::find(result.a2aMatrix.begin(), result.a2aMatrix.end(), m) == result.a2aMatrix.end())
                    {
                        result.a2aMatrix.push_back(m);
                        tokenReplace(result.a2aMatrix.back().file, "traj", traj);
                    }
                }
                result.contraction = p;
                result.correlator.resize(par.global.nt, 0.);

                translations = parseTimeRange(p.translations, par.global.nt);
                makeTimeSeq(timeSeq, times);
                LOG(Message) << timeSeq.size()*translations.size()*(term.size() - 2) << " A*B, "
                        << timeSeq.size()*translations.size()*par.global.nt << " tr(A*B)"
                        << std::endl;

                LOG(Message) << "* Caching transposed last term" << std::endl;
                for (unsigned int t = 0; t < par.global.nt; ++t)
                {
                    tAr.startTimer("Disk vector overhead");
                    ref = a2aMat.at(term.back())[t];
                    tAr.stopTimer("Disk vector overhead");

                    tAr.startTimer("Transpose caching");
                    lastTerm[t].resize(ref.rows(), ref.cols());
                    thread_for( j,ref.cols(),{
                      for (unsigned int i = 0; i < ref.rows(); ++i)
                      {
                        if (evalMap.count(term.back()) > 0 && j < evalMap.at(term.back()).size() ) {
                            lastTerm[t](i, j) = ref(i, j)*evalMap.at(term.back())[j];
                        } else {
                            lastTerm[t](i, j) = ref(i, j);
                        }
                      }
            		});
                    tAr.stopTimer("Transpose caching");
                }
                bytes = par.global.nt*lastTerm[0].rows()*lastTerm[0].cols()*sizeof(ComplexD);
                LOG(Message) << Sec(tAr.getDTimer("Transpose caching")) << " " 
                          << Bytes(bytes, tAr.getDTimer("Transpose caching")) << std::endl;
                for (unsigned int i = 0; i < timeSeq.size(); ++i)
                {
                    unsigned int dti = 0;
                    auto         &t = timeSeq[i];

                    result.times = t;
                    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                    {
                        result.correlator[tLast] = 0.;
                    }
                    for (auto &dt: translations)
                    {
                        LOG(Message) << "* Step " << i*translations.size() + dti + 1
                                << "/" << timeSeq.size()*translations.size()
                                << " -- positions= " << t << ", dt= " << dt << std::endl;
                        if (term.size() > 2)
                        {
                            std::cout << std::setw(8) << "products";
                        }
                        flops  = 0.;
                        bytes  = 0.;
                        fusec  = tAr.getDTimer("A*B algebra");
                        busec  = tAr.getDTimer("A*B total");
                        tAr.startTimer("Linear algebra");
                        tAr.startTimer("Disk vector overhead");
                        prod = a2aMat.at(term[0])[TIME_MOD(t[0] + dt)];
                        tAr.stopTimer("Disk vector overhead");

                        if (evalMap.count(term[0]) > 0) {
                            thread_for_collapse2( j,prod.rows(),{
                                for(uint64_t k = 0; k < evalMap.at(term[0]).size(); k++) {
                                    prod(j, k) = prod(j, k)*evalMap.at(term[0])[k];
                                }
                            });
                        }

                        for (unsigned int j = 1; j < term.size() - 1; ++j)
                        {
                            tAr.startTimer("Disk vector overhead");
                            ref = a2aMat.at(term[j])[TIME_MOD(t[j] + dt)];
                            tAr.stopTimer("Disk vector overhead");
                            
                            if (evalMap.count(term[j]) > 0) {
                                thread_for(k,ref.rows(),{
                                    for(int l = 0; l < evalMap.at(term[j]).size(); l++){
                                        ref(k, l) = ref(k, l)*evalMap.at(term[j])[l];
                                    }
                                });
                            }

                            tAr.startTimer("A*B total");
                            tAr.startTimer("A*B algebra");
                            A2AContraction::mul(tmp, prod, ref);
                            tAr.stopTimer("A*B algebra");
                            flops += A2AContraction::mulFlops(prod, ref);
                            prod   = tmp;
                            tAr.stopTimer("A*B total");
                            bytes += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
                        }
                        if (term.size() > 2)
                        {
                            std::cout << Sec(tAr.getDTimer("A*B total") - busec) << " "
                                    << Flops(flops, tAr.getDTimer("A*B algebra") - fusec) << " " 
                                    << Bytes(bytes, tAr.getDTimer("A*B total") - busec) << std::endl;
                        }
                        LOG(Message) << std::setw(8) << "traces";
                        flops  = 0.;
                        bytes  = 0.;
                        fusec  = tAr.getDTimer("tr(A*B)");
                        busec  = tAr.getDTimer("tr(A*B)");
                        for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                        {
                            tAr.startTimer("tr(A*B)");
                            A2AContraction::accTrMul(result.correlator[TIME_MOD(tLast - dt)], prod, lastTerm[tLast]);
                            tAr.stopTimer("tr(A*B)");
                            flops += A2AContraction::accTrMulFlops(prod, lastTerm[tLast]);
                            bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
                        }
                        tAr.stopTimer("Linear algebra");
                        std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
                                << Flops(flops, tAr.getDTimer("tr(A*B)") - fusec) << " " 
                                << Bytes(bytes, tAr.getDTimer("tr(A*B)") - busec) << std::endl;

                        if (!p.translationAverage)
                        {
                            if (p.spaceNormalize) {
                                for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                                {
                                    result.correlator[tLast] /= vol;
                                }
                            }

                            saveCorrelator(result, par.global.output, dt, traj);
                            for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                            {
                                result.correlator[tLast] = 0.;
                            }
                        }
                        dti++;
                    }
                    if (p.translationAverage)
                    {
                        for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                        {
                            if (p.spaceNormalize) {
                                result.correlator[tLast] /= vol*translations.size();
                            } else {
                                result.correlator[tLast] /= translations.size();
                            }
                        }
                        saveCorrelator(result, par.global.output, 0, traj);
                    }
                }
                tAr.stopTimer("Total");
                printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
            }
        }
    }
    grid->Barrier();
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    
    return EXIT_SUCCESS;
}
