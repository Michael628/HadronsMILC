/*
 * EigenPackExtractEvalsMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef HadronsMILC_MUtilities_EigenPackExtractEvals_hpp_
#define HadronsMILC_MUtilities_EigenPackExtractEvals_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Load eigen vectors/values package                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EigenPackExtractEvalsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EigenPackExtractEvalsPar,
                                    std::string,  eigenPack,
                                    std::string,  output);
};

template <typename FImpl, typename Pack>
class TEigenPackExtractEvals: public Module<EigenPackExtractEvalsPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef BaseEigenPack<Field>   BasePack;

    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TEigenPackExtractEvals(const std::string name);
    // destructor
    virtual ~TEigenPackExtractEvals(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(EigenPackExtractEvals, ARG(TEigenPackExtractEvals<STAGIMPL,BaseFermionEigenPack<STAGIMPL> >), MUtilities);

/******************************************************************************
 *                    TEigenPackExtractEvals implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TEigenPackExtractEvals<FImpl,Pack>::TEigenPackExtractEvals(const std::string name)
: Module<EigenPackExtractEvalsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TEigenPackExtractEvals<FImpl,Pack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TEigenPackExtractEvals<FImpl,Pack>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename FImpl, typename Pack>
DependencyMap TEigenPackExtractEvals<FImpl, Pack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TEigenPackExtractEvals<FImpl,Pack>::setup(void)
{
    if (par().output.empty()) {
        HADRONS_ERROR(Parsing, "Must provide output file.");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TEigenPackExtractEvals<FImpl,Pack>::execute(void)
{
    auto &epack = envGet(BasePack, par().eigenPack);

    saveResult(par().output,"evals",epack.eval);

    // if (env().getGrid()->IsBoss())
    // {
    //     makeFileDir(par().output, env().getGrid());
    //     {
    //         ResultWriter WR(resultFilename(par().output));
    //         write(WR,"evals",EF);
    //     }
    // }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_EigenPackExtractEvals_hpp_
