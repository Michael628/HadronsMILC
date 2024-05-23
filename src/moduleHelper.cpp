#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;


int myCreateModule(Application &app, std::string name, std::string type, XmlReader& reader) {

  int status = 1;

  LOG(Message) << "Building " << name << std::endl;

  if (type == "MSolver::StagLocalCoherenceLanczos300") {
    MSolver::TLocalCoherenceLanczos<STAGIMPL,300> module(name);

    module.parseParameters(reader,"options");

    app.createModule<MSolver::TLocalCoherenceLanczos<STAGIMPL,300> >(name, module.par());
    /*  } else if (type == "MContraction::StagA2AMesonField") {
    MContraction::TNewMesonField<STAGIMPL,MassShiftEigenPack<STAGIMPL> > module(name);

    module.parseParameters(reader,"options");

    app.createModule<MContraction::TNewMesonField<STAGIMPL,MassShiftEigenPack<STAGIMPL> > >(name, module.par());
    */
  // } else if (type == "MSolver::StagLMA") {
  //   MSolver::TMyLowModeProjMILC<STAGIMPL,MassShiftEigenPack<STAGIMPL> > module(name);

  //   module.parseParameters(reader,"options");

  //   app.createModule<MSolver::TMyLowModeProjMILC<STAGIMPL,MassShiftEigenPack<STAGIMPL> > >(name, module.par());

    //  } else if (type == "MIO::LoadIldg") {
    //    MIO::TLoadIldg<GIMPL> module(name);

    //    module.parseParameters(reader,"options");

    //    app.createModule<MIO::TLoadIldg<GIMPL> >(name,module.par());

    //  } else if (type == "MSolver::StagA2AVectors") {
    //    MSolver::TProbe<STAGIMPL, BaseFermionEigenPack<STAGIMPL>> module(name);

    //    module.parseParameters(reader,"options");

    //    app.createModule<MSolver::TProbe<STAGIMPL, BaseFermionEigenPack<STAGIMPL> > >(name, module.par());


    //    } else if (type == "MContraction::A2AMesonFieldMILC") {
    //      MContraction::TMyMesonField<STAGIMPL> module(name);

    //      module.parseParameters(reader,"options");

    //      app.createModule<MContraction::TMyMesonField<STAGIMPL> >(name, module.par());
  // } else if (type == "MSolver::StagMixedPrecisionRBPrecCG") {
  //   MSolver::TMixedPrecisionRBPrecCG<STAGIMPLF, STAGIMPLD, 380> module(name);

  //   module.parseParameters(reader,"options");

  //   app.createModule<MSolver::TMixedPrecisionRBPrecCG<STAGIMPLF, STAGIMPLD, 380> >(name, module.par());

  // } else if (type == "MNoise::StagTimeDilutedSpinColorDiagonal") {
  //   MNoise::TTimeDilutedSpinColorDiagonal<STAGIMPL> module(name);

  //   module.parseParameters(reader,"options");

  //   app.createModule<MNoise::TTimeDilutedSpinColorDiagonal<STAGIMPL> >(name);

  // } else if (type == "MIO::StagLoadFermionEigenPack") {
  //   MIO::TLoadEigenPack<FermionEigenPack<STAGIMPL>, GIMPL> module(name);

  //   module.parseParameters(reader,"options");

  //   app.createModule<MIO::TLoadEigenPack<FermionEigenPack<STAGIMPL>, GIMPL> >(name,module.par());

  // } else if (type == "MNoise::StagSparseSpinColorDiagonal") {
  //   MNoise::TSparseSpinColorDiagonal<STAGIMPL> module(name);

  //   module.parseParameters(reader,"options");

  //   app.createModule<MNoise::TSparseSpinColorDiagonal<STAGIMPL> >(name, module.par());

  // } else if (type == "MNoise::StagFullVolumeSpinColorDiagonal") {
  //   MNoise::TFullVolumeSpinColorDiagonal<STAGIMPL> module(name);

  //   module.parseParameters(reader,"options");

  //   app.createModule<MNoise::TFullVolumeSpinColorDiagonal<STAGIMPL> >(name, module.par());

  // } else if (type == "MSolver::StagA2AVectors") {
  //   MSolver::TA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>> module(name);

  //   module.parseParameters(reader,"options");

  //   app.createModule<MSolver::TA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL> > >(name, module.par());

  // } else if (type == "MIO::LoadGauge") {
  //   MIO::TLoadField<GIMPL::GaugeField> module(name);

  //   module.parseParameters(reader,"options");

  //   app.createModule<MIO::TLoadField<GIMPL::GaugeField> >(name, module.par());

  } else {
    status = 0;
  }

  return status;
}
