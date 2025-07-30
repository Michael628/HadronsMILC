# AGENTS.md - HadronsMILC Development Guide

## Build Commands
- **Build**: `./configure.sh && make` (uses autotools/automake)
- **Clean**: `make clean`
- **Python utilities**: `cd utilities/pyfm && pip install -e .`

## Code Style - C++
- Use Grid/Hadrons framework conventions with `BEGIN_HADRONS_NAMESPACE`/`END_HADRONS_NAMESPACE`
- Template classes prefixed with `T` (e.g., `TImprovedStaggeredMILC`, `TMesonMILC`)
- Module registration with `MODULE_REGISTER_TMP` macro
- Include guards: `#ifndef HadronsMILC_ModuleName_hpp_`
- Standard includes: `<Hadrons/Global.hpp>`, `<Hadrons/Module.hpp>`, `<Hadrons/ModuleFactory.hpp>`
- Use `GRID_SERIALIZABLE_CLASS_MEMBERS` for parameter classes
- Logging with `LOG(Message)` for output
- Error handling with `HADRONS_ERROR(Type, message)`

## Code Style - Python
- Use type hints and Pydantic models for configuration
- Import style: standard library first, then third-party, then local imports
- Use `setup_logging()` for debug output
- Configuration via YAML files with `utils.load_param()`

## Project Structure
- `src/HadronsMILC/Modules/`: C++ physics modules organized by category (MAction, MContraction, etc.)
- `utilities/pyfm/`: Python utilities for job management and data processing
- Main executable built from `main.cpp` using Hadrons Application framework