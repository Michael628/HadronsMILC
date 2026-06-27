#include "HighlyImprovedStaggered.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MAction;

template class Grid::Hadrons::MAction::THighlyImprovedStaggeredMILC<STAGIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MAction::THighlyImprovedStaggeredMILC<STAGIMPLF>;
#endif
