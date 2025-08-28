#include "ImprovedStaggered.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MAction;

template class Grid::Hadrons::MAction::TImprovedStaggeredMILC<STAGIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MAction::TImprovedStaggeredMILC<STAGIMPLF>;
#endif
