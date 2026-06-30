#include "HISQSmear.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

template class Grid::Hadrons::MGauge::THISQSmear<STAGIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MGauge::THISQSmear<STAGIMPLF>;
#endif
