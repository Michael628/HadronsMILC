#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include <Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " [parameter file] [Grid options]"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());
  LOG(Message) << "Grid initialized" << std::endl;

  // run setup ///////////////////////////////////////////////////////////////
  std::string paramFile = argv[1];
  Application app(paramFile);

  // execution — run() auto-calls parseParameterFile() because the filename
  // constructor set parameterFileName_ and no modules exist yet. This reads
  // <parameters>, creates all modules from <modules> (including optional
  // <subgrid> tags), then executes the trajectory loop.
  app.run();

  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();

  return EXIT_SUCCESS;
}
