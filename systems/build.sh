#! /bin/bash

# Expected directory structure:
# TOPDIR (root directory where this script is run from)
# ├── HadronsMILC/           - Main application (git repo)
# │   ├── systems/           - System-specific configurations
# │   │   ├── ${CONFIG_SYSTEM}/
# │   │   │   ├── env.sh     - System environment configuration
# │   │   │   └── configure-params.sh - System build parameters
# │   │   └── configure-params_default.sh - Default configuration parameters
# │   ├── build${BUILD_EXT}/ - Build directory
# │   └── install${BUILD_EXT}/ - Installation directory
# ├── Grid/                  - Grid library (git repo, cloned from GitHub)
# │   ├── build${BUILD_EXT}/ - Build directory
# │   └── install${BUILD_EXT}/ - Installation directory
# ├── Hadrons/               - Hadrons library (git repo, cloned from GitHub)
# │   ├── build${BUILD_EXT}/ - Build directory
# │   └── install${BUILD_EXT}/ - Installation directory
# ├── deps/                  - Dependencies directory
# │   ├── gmp-6.2.1/         - GMP source (downloaded)
# │   ├── gmp/               - Symlink to gmp-6.2.1
# │   ├── mpfr-4.1.0/        - MPFR source (downloaded)
# │   ├── mpfr/              - Symlink to mpfr-4.1.0
# │   ├── lime-1.3.2/        - LIME source (downloaded)
# │   ├── lime/              - Symlink to lime-1.3.2
# │   ├── hdf5-1.10.10/      - HDF5 source (downloaded)
# │   ├── hdf5/              - Symlink to hdf5-1.10.10
# │   ├── openssl-3.3.1/     - OpenSSL source (downloaded)
# │   ├── openssl/           - Symlink to openssl-3.3.1
# │   └── install-${CONFIG_SYSTEM}/ - Installation directory for dependencies
# ├── env${BUILD_EXT}.sh     - Optional local environment configuration
# └── configure-params${BUILD_EXT}.sh - Optional local configuration parameters

TOPDIR=$(pwd)

function print_help() {
  cat << 'EOF'
Usage: ./build.sh [OPTIONS]

Dependency Options (build external libraries):
  --gmp              Build GMP library
  --mpfr             Build MPFR library
  --lime             Build LIME library
  --hdf5             Build HDF5 library
  --ssl              Build OpenSSL library

Build Configuration Options:
  --system <name>    Configure system (default: scalar)
  --ext <name>       Add extension suffix to build directories
  --threads <num>    Number of parallel threads for build (default: 4)
  --debug            Enable debug build
  --mpi-reduction    Enable MPI reduction
  --force            Force rebuild (clean build directories)
  --skip-make        Skip make step, only configure

Component Selection:
  --grid             Build Grid library
  --hadrons          Build Hadrons library
  --app              Build HadronsMILC application
  --all              Build all components (grid, hadrons, app)

Help:
  --help             Show this help message

Examples:
  ./build.sh --system perlmutter --app --threads 8
  ./build.sh --gmp --mpfr --system scalar --grid --hadrons
EOF
}

function parse_flags() {
  # Initialize all flags to false
  BUILD_DEBUG='false'
  BUILD_MPI_REDUCTION='false'
  BUILD_GMP='false'
  BUILD_MPFR='false'
  BUILD_LIME='false'
  BUILD_HDF5='false'
  BUILD_OPENSSL='false'
  FORCE_REBUILD='false'
  MAKE_COMPONENT='true'
  BUILD_EXT=''
  CONFIG_SYSTEM='scalar'
  THREADS=4
  BUILD_COMPONENTS=''

  # Parse all arguments
  while test $# -gt 0; do
    echo "param: $1"
    case "$1" in
      # Dependency flags
      --gmp)
        BUILD_GMP='true'
        shift
      ;;
      --threads)
        shift
        THREADS=$1
        shift
      ;;
      --ssl)
        BUILD_OPENSSL='true'
        shift
      ;;
      --mpfr)
        BUILD_MPFR='true'
        shift
      ;;
      --lime)
        BUILD_LIME='true'
        shift
      ;;
      --ext)
        shift
        BUILD_EXT="-$1"
        shift
      ;;
      --hdf5)
        BUILD_HDF5='true'
        shift
      ;;
      --all)
        BUILD_COMPONENTS="grid hadrons app"
        shift
      ;;

      # Build flags
      --debug)
        BUILD_DEBUG='true'
        shift
      ;;
      --mpi-reduction)
        BUILD_MPI_REDUCTION='true'
        shift
      ;;
      --force)
        FORCE_REBUILD='true'
        shift
      ;;
      --system)
        shift
        CONFIG_SYSTEM="$1"
        shift
      ;;
      --skip-make)
        MAKE_COMPONENT='false'
        shift
      ;;

      # Components
      --grid|--hadrons|--app)
          BUILD_COMPONENTS="${BUILD_COMPONENTS} ${1:2}"
        shift
      ;;

      --help)
        print_help
        exit 0
      ;;

      *)
        echo "Error: Unknown argument: $1"
        echo ""
        print_help
        exit 1
      ;;
    esac
  done

  BUILD_EXT="-${CONFIG_SYSTEM}${BUILD_EXT}"
}

function source_env() {
  # Load local environment if it exists, otherwise load system
  if [ -f "${TOPDIR}/env${BUILD_EXT}.sh" ]; then
    source ${TOPDIR}/env${BUILD_EXT}.sh
  elif [ -f "${TOPDIR}/HadronsMILC/systems/${CONFIG_SYSTEM}/env.sh" ]; then
    source ${TOPDIR}/HadronsMILC/systems/${CONFIG_SYSTEM}/env.sh
  fi
}

function source_config_params() {
  # Load default params, then overwrite with system and local params
  source ${TOPDIR}/HadronsMILC/systems/configure-params_default.sh
  # Overwrite with system params if they exist
  [ -f "${TOPDIR}/HadronsMILC/systems/${CONFIG_SYSTEM}/configure-params.sh" ] && source ${TOPDIR}/HadronsMILC/systems/${CONFIG_SYSTEM}/configure-params.sh
  # Overwrite with local params if they exist
  [ -f "${TOPDIR}/configure-params${BUILD_EXT}.sh" ] && source ${TOPDIR}/configure-params${BUILD_EXT}.sh
}

function dependencies() {
  # Installs gmp, mpfr, lime, and (sometimes) hdf5
  # Uses global BUILD_* variables set by parse_flags()

  WORKDIR=${TOPDIR}/deps
  INSTALLDIR=${WORKDIR}/install-${CONFIG_SYSTEM}

  mkdir -p ${WORKDIR}
  pushd ${WORKDIR}

  source_env
  module list

  source_config_params

  if [ $BUILD_OPENSSL = 'true' ]; then
          if [ ! -d openssl-3.3.1 ]; then
            wget https://www.openssl.org/source/openssl-3.3.1.tar.gz
            tar -xf openssl-3.3.1.tar.gz
            ln -s openssl-3.3.1/ openssl
            rm openssl-3.3.1.tar.gz
          fi

          pushd openssl-3.3.1
          rm -rf build${BUILD_EXT}
          mkdir -p build${BUILD_EXT}
          pushd build${BUILD_EXT}
          dependency_configure "openssl" "${INSTALLDIR}"
          make all install
          status=$?

          popd
          if [[ $status -ne 0 ]]; then
                  echo "openssl compile failed."
                  exit 1
          fi
          popd
  fi

  if [ $BUILD_GMP = 'true' ]; then
          if [ ! -d gmp-6.2.1 ]; then
            wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz
            tar -xf gmp-6.2.1.tar.xz
            ln -s gmp-6.2.1/ gmp
            rm gmp-6.2.1.tar.xz
          fi

          pushd gmp-6.2.1
          rm -rf build${BUILD_EXT}
          mkdir -p build${BUILD_EXT}
          pushd build${BUILD_EXT}
          dependency_configure "gmp" "${INSTALLDIR}"
          make all install
          status=$?

          popd
          if [[ $status -ne 0 ]]; then
                  echo "gmp compile failed."
                  exit 1
          fi
          popd
  fi

  if [ $BUILD_MPFR == 'true' ]; then
          if [ ! -d mpfr-4.1.0 ]; then
            wget https://www.mpfr.org/mpfr-4.1.0/mpfr-4.1.0.tar.gz
            tar -xvzf mpfr-4.1.0.tar.gz
            ln -s mpfr-4.1.0/ mpfr
            rm mpfr-4.1.0.tar.gz
          fi

          pushd mpfr-4.1.0
          rm -rf build${BUILD_EXT}
          mkdir -p build${BUILD_EXT}
          pushd build${BUILD_EXT}
          dependency_configure "mpfr" "${INSTALLDIR}"
          make all install
          status=$?
          popd

          if [ $status -ne 0 ]; then
                  echo "mpfr compile failed."
                  exit 1
          fi
          popd
  fi

  if [ $BUILD_LIME == 'true' ]; then
          if [ ! -d lime-1.3.2 ]; then
            wget http://usqcd-software.github.io/downloads/c-lime/lime-1.3.2.tar.gz
            tar -xvzf lime-1.3.2.tar.gz
            ln -s lime-1.3.2/ lime
            rm lime-1.3.2.tar.gz
          fi

          pushd lime-1.3.2
          rm -rf build${BUILD_EXT}
          mkdir -p build${BUILD_EXT}
          pushd build${BUILD_EXT}
          dependency_configure "lime" "${INSTALLDIR}"
          make all install
          status=$?
          popd

          if [ $status -ne 0 ]; then
                  echo "lime compile failed."
                  exit 1
          fi
          popd
  fi

  if [ $BUILD_HDF5 == 'true' ]; then
          if [ ! -d hdf5-1.10.10 ]; then
            wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.10/src/hdf5-1.10.10.tar.gz
            tar -xvzf hdf5-1.10.10.tar.gz
            ln -s hdf5-1.10.10/ hdf5
            rm hdf5-1.10.10.tar.gz
          fi

          pushd hdf5-1.10.10
          rm -rf build${BUILD_EXT}
          mkdir -p build${BUILD_EXT}
          pushd build${BUILD_EXT}
          dependency_configure "hdf5" "${INSTALLDIR}"
          make all install
          status=$?
          popd
          if [ $status -ne 0 ]; then
                  echo "hdf5 compile failed."
                  exit 1
          fi
          popd
  fi

  popd

}

function build-component() {
  local SOURCE=$1

  # Uses global BUILD_* variables set by parse_flags()

	case ${SOURCE} in
    grid)
      GIT_REPO=https://github.com/milc-qcd/Grid
      GIT_BRANCH="feature/LMI-develop"
      SRCDIR=${TOPDIR}/Grid
      ;;
    hadrons)
      GIT_REPO=https://github.com/milc-qcd/Hadrons
      GIT_BRANCH="feature/LMI-develop"
      SRCDIR=${TOPDIR}/Hadrons
      ;;
    app)
      GIT_REPO=https://github.com/Michael628/HadronsMILC
      GIT_BRANCH="develop"
      SRCDIR=${TOPDIR}/HadronsMILC
      ;;
    *)
      echo "Unsupported build type"
      echo "Usage $0 <grid|hadrons|app>"
      exit 1
    esac

	BUILDDIR=${SRCDIR}/build${BUILD_EXT}
	INSTALLDIR=${SRCDIR}/install${BUILD_EXT}

	if [ ! -d ${SRCDIR} ]
	then
	  mkdir -p ${SRCDIR}
	  echo "Fetching ${GIT_BRANCH} branch of ${SOURCE} package from github"
	  pushd ${TOPDIR}
	  git clone ${GIT_REPO} -b ${GIT_BRANCH}
    git submodule update --init
	  popd
	fi

	# Fetch Eigen package, set up Make.inc files and create Grid configure
	pushd ${SRCDIR}
	./bootstrap.sh
	popd

	# Configure only if not already configured
	mkdir -p ${BUILDDIR}
  if [ "$FORCE_REBUILD" = 'true' ]; then
    rm -rf ${BUILDDIR}/*
  fi
	pushd ${BUILDDIR}
	if [ ! -f Makefile ]
	then
	  echo "Configuring ${SOURCE} in ${BUILDDIR}"

    source_env
    module list > compile.out 2>&1

    source_config_params

    case ${SOURCE} in
      grid)
        grid_configure "${INSTALLDIR}" "${TOPDIR}" >> compile.out 2>&1
        status=$?
        echo "Configure exit status $status"
      ;;
      hadrons)
        hadrons_configure "${INSTALLDIR}" "${TOPDIR}" >> compile.out 2>&1
        status=$?
        echo "Configure exit status $status"
      ;;
      app)
        app_configure "${INSTALLDIR}" "${TOPDIR}" >> compile.out 2>&1
        status=$?
        echo "Configure exit status $status"
      ;;
    esac

  if [ "$MAKE_COMPONENT" = 'true' ]; then

	  if [ $status -ne 0 ]
	  then
	      echo "Quitting because of configure errors"
	  else
	    echo "Building in ${BUILDDIR}"
	    make V=1 -k -j$THREADS  >> compile.out 2>&1

	    echo "Installing in ${INSTALLDIR}"
	    make install >> compile.out 2>&1
	  fi

	fi
	fi
	popd
}

# Show help if no arguments provided
if [ $# -eq 0 ]; then
  print_help
  exit 0
fi

# Parse all command-line flags first
parse_flags "$@"

# Build dependencies if any were requested
if [ "$BUILD_GMP" = 'true' ] || [ "$BUILD_MPFR" = 'true' ] || \
   [ "$BUILD_LIME" = 'true' ] || [ "$BUILD_HDF5" = 'true' ] || \
   [ "$BUILD_OPENSSL" = 'true' ]; then
  dependencies
fi

# Build each requested component
for component in $BUILD_COMPONENTS; do
  build-component "$component"
done
