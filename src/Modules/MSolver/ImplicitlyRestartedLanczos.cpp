#include "ImplicitlyRestartedLanczos.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TImplicitlyRestartedLanczosMILC<STAGIMPL::FermionField>;
template class Grid::Hadrons::MSolver::TImplicitlyRestartedLanczosMILC<STAGIMPL::FermionField, STAGIMPLF::FermionField>;
