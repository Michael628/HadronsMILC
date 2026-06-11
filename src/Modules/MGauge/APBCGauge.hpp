/*
 * APBCGauge.hpp, part of HadronsMILC
 */
#ifndef HadronsMILC_MGauge_APBCGauge_hpp_
#define HadronsMILC_MGauge_APBCGauge_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MGauge)

class APBCGaugePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(APBCGaugePar,
                                    std::string, gauge,
                                    std::string, boundary);
};

template <typename GField>
class TAPBCGauge: public Module<APBCGaugePar>
{
public:
    // constructor
    TAPBCGauge(const std::string name);
    // destructor
    virtual ~TAPBCGauge(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(APBCGauge, TAPBCGauge<LatticeGaugeField>, MGauge);

/******************************************************************************
 *                  TAPBCGauge implementation                                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GField>
TAPBCGauge<GField>::TAPBCGauge(const std::string name)
: Module<APBCGaugePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GField>
std::vector<std::string> TAPBCGauge<GField>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GField>
std::vector<std::string> TAPBCGauge<GField>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GField>
void TAPBCGauge<GField>::setup(void)
{
    envCreateLat(GField, getName());
    envTmpLat(LatticeColourMatrix, "link");
    envTmp(Lattice<iScalar<vInteger>>, "coord", 1, envGetGrid(GField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GField>
void TAPBCGauge<GField>::execute(void)
{
    const unsigned int nd = env().getNd();
    std::vector<int> boundary(nd, 1);

    boundary[nd - 1] = -1;
    if (!par().boundary.empty())
    {
        boundary = strToVec<int>(par().boundary);
    }
    if (boundary.size() != nd)
    {
        HADRONS_ERROR(Argument, "boundary has " + std::to_string(boundary.size())
                      + " entries, expected " + std::to_string(nd));
    }

    LOG(Message) << "Creating APBC gauge field '" << getName()
                 << "' from '" << par().gauge << "'" << std::endl;

    auto &in  = envGet(GField, par().gauge);
    auto &out = envGet(GField, getName());

    out = in;

    envGetTmp(LatticeColourMatrix, link);
    envGetTmp(Lattice<iScalar<vInteger>>, coord);

    for (int mu = 0; mu < Nd; ++mu)
    {
        LOG(Message) << "Boundary phase[" << mu << "] = "
                     << boundary[mu] << std::endl;

        if (boundary[mu] != 1)
        {
            LatticeCoordinate(coord, mu);
            link = PeekIndex<LorentzIndex>(out, mu);
            int dimSize = out.Grid()->GlobalDimensions()[mu] - 1;
            link = where((coord == dimSize),
                         static_cast<double>(boundary[mu])*link,
                         link);
            PokeIndex<LorentzIndex>(out, link, mu);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MGauge_APBCGauge_hpp_
