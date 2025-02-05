#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " [parameter file] [Grid options]" << std::endl;
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
    Application app;
    Application::GlobalPar par;
    Application::ObjectId id;
    
    std::string vecDesc, paramFile;

    paramFile = argv[1];
    XmlReader reader(paramFile, false, HADRONS_XML_TOPLEV);

    read(reader, "parameters", par);
    app.setPar(par);

    if (!push(reader, "modules"))
    {
        HADRONS_ERROR(Parsing, "Cannot open node 'modules' in parameter file '" 
                            + paramFile + "'");
    }
    if (!push(reader, "module"))
    {
        HADRONS_ERROR(Parsing, "Cannot open node 'modules/module' in parameter file '" 
                            + paramFile + "'");
    }
    do
    {
        read(reader, "id", id);
        // if (!myCreateModule(app, id.name, id.type, reader)) {
            app.createModule(id.name, id.type, reader);
        // }
    } while (reader.nextElement("module"));
    pop(reader);
    pop(reader);

    // execution
    app.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
