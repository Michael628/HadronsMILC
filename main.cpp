#include "config.h"

// Grid pulls in lime_config.h, which #undefs PACKAGE_VERSION (it #undefs
// before and after using its own LIME version). Capture the literal value
// under a project-specific name before any Grid/Hadrons include.
static const char *const HADRONSMILC_VERSION = PACKAGE_VERSION;

#include <Grid/Version.h>
#include <GridMilc/util/Version.h>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include <Modules.hpp>

#include <string>
#include <vector>

#ifndef GITHASH
#define GITHASH "unknown"
#endif

using namespace Grid;
using namespace Hadrons;

namespace {

// Stable key=value identity lines. pyfm's audit parses these from run logs, so
// change the keys only together with the parser.
std::vector<std::string> versionLines(void) {
  return {
      std::string("HadronsMILC version=") + HADRONSMILC_VERSION +
          " git=" + HADRONSMILC_GIT,
      std::string("Grid git=") + GITHASH,
      std::string("Hadrons git=") + HADRONSMILC_HADRONS_GIT +
          " (configure-time)",
      std::string("GridMilc version=") + gridMilcVersion(),
      std::string("Dependency pins=") + HADRONSMILC_PINS,
  };
}

// Optional pyfm provenance record, <grid><provenance>. Every leaf is optional.
class ProvenancePar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ProvenancePar,
                                  std::string, pyfmVersion,
                                  std::string, pyfmSha,
                                  std::string, hadronsMilcCompat,
                                  std::string, generated);
};

// Guarded read (Application.cpp:274): push() returns false without a warning
// when the leaf is absent.
std::string readOptional(XmlReader &reader, const std::string &name) {
  std::string value;
  if (push(reader, name)) {
    pop(reader);
    read(reader, name, value);
  }
  return value;
}

// "v0.2.3" -> "0.2", "0.2" -> "0.2", "1" -> "1".
std::string majorMinor(std::string version) {
  if (!version.empty() && (version[0] == 'v' || version[0] == 'V')) {
    version.erase(0, 1);
  }
  const auto first = version.find('.');
  if (first == std::string::npos) {
    return version;
  }
  return version.substr(0, version.find('.', first + 1));
}

void checkProvenance(const std::string &paramFile) {
  XmlReader reader(paramFile, false, HADRONS_XML_TOPLEV);
  if (!push(reader, "provenance")) {
    return;
  }
  ProvenancePar par;
  par.pyfmVersion       = readOptional(reader, "pyfmVersion");
  par.pyfmSha           = readOptional(reader, "pyfmSha");
  par.hadronsMilcCompat = readOptional(reader, "hadronsMilcCompat");
  par.generated         = readOptional(reader, "generated");
  pop(reader);

  LOG(Message) << "Provenance pyfmVersion=" << par.pyfmVersion
               << " pyfmSha=" << par.pyfmSha
               << " hadronsMilcCompat=" << par.hadronsMilcCompat
               << " generated=" << par.generated << std::endl;
  if (par.hadronsMilcCompat.empty()) {
    return;
  }
  const std::string ours   = majorMinor(HADRONSMILC_VERSION);
  const std::string theirs = majorMinor(par.hadronsMilcCompat);
  if (ours != theirs) {
    LOG(Warning) << "Parameter file targets HadronsMILC " << theirs
                 << " but this binary is " << HADRONSMILC_VERSION
                 << " (MAJOR.MINOR mismatch)" << std::endl;
  }
}

} // namespace

int main(int argc, char *argv[]) {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--version") {
      for (const auto &line : versionLines()) {
        std::cout << line << std::endl;
      }
      return EXIT_SUCCESS;
    }
  }

  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " [parameter file] [Grid options]"
              << std::endl;
    std::cerr << "       " << argv[0] << " --version" << std::endl;
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
  for (const auto &line : versionLines()) {
    LOG(Message) << line << std::endl;
  }

  // run setup ///////////////////////////////////////////////////////////////
  std::string paramFile = argv[1];
  checkProvenance(paramFile);
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
