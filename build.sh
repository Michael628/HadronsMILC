#! /bin/bash

echo "CC=$CC"
echo "CPP=$CPP"
echo "CXX=$CXX"

force_flag='false'

print_usage() {
	printf "Usage: $0 [OPTION] ARCH SOURCE"
}

while test $# -gt 0; do
	case "$1" in
	-f | --force)
		force_flag='true'
		shift
		;;
	deps | grid | hadrons | app)
		SOURCE=$1
		shift
		;;
	scalar | mpi)
		ARCH=$1
		shift
		;;
	*)
		break
		;;
	esac
done

case ${ARCH} in
scalar | mpi) ;;
*)
	echo "Unsupported ARCH"
	echo "Usage $0 <scalar|mpi> <deps|grid|hadrons|app>"
	exit 1
	;;
esac

TOPDIR=$(pwd)

case ${SOURCE} in
grid)
	GIT_REPO=https://github.com/milc-qcd/Grid
	case ${ARCH} in
	old)
		GIT_BRANCH="feature/staggered-a2a-ml"
		;;
	*)
		GIT_BRANCH="feature/LMI-develop"
		;;
	esac
	SRCDIR=${TOPDIR}/Grid
	;;
hadrons)
	GIT_REPO=https://github.com/milc-qcd/Hadrons
	case ${ARCH} in
	old)
		GIT_BRANCH="feature/staggered-a2a-ml"
		;;
	*)
		GIT_BRANCH="feature/LMI-develop"
		;;
	esac
	SRCDIR=${TOPDIR}/Hadrons
	;;
app)
	GIT_REPO=https://github.com/Michael628/HadronsMILC
	GIT_BRANCH="main"
	SRCDIR=${TOPDIR}/HadronsMILC
	;;
deps)
	if [ ! -f deps.sh ]; then
		echo "Missing deps.sh file in directory ${TOPDIR}."
		exit 1
	else
		echo "Running ${TOPDIR}/deps.sh ${ARCH}."
		${TOPDIR}/deps.sh --all ${ARCH}
		status=?
		if [ $status -ne 0 ]; then
			echo "Dependency script failed. See output."
		fi
		exit 0
	fi
	;;
*)
	echo "Unsupported build type"
	echo "Usage $0 <scalar|mpi> <grid|hadrons|app>"
	exit 1
	;;
esac

BUILDDIR=${SRCDIR}/build-${ARCH}
INSTALLDIR=${SRCDIR}/install-${ARCH}

MAKE="make V=1"

if [ ! -d ${SRCDIR} ]; then
	mkdir -p ${SRCDIR}
	echo "Fetching ${GIT_BRANCH} branch of ${SOURCE} package from github"
	pushd ${TOPDIR}
	git clone ${GIT_REPO} -b ${GIT_BRANCH}
	popd
fi

# Fetch Eigen package, set up Make.inc files and create Grid configure
pushd ${SRCDIR}
git checkout ${GIT_BRANCH}
./bootstrap.sh
popd

# Configure only if not already configured
rm -rf ${BUILDDIR}
mkdir -p ${BUILDDIR}
pushd ${BUILDDIR}
if [ ! -f Makefile ] || [ "${force_flag}" == 'true' ]; then
	echo "Configuring ${SOURCE} for ${ARCH} in ${BUILDDIR}"

	case ${ARCH} in

	mpi)
		case ${SOURCE} in
		grid)
			${SRCDIR}/configure \
				--prefix=${INSTALLDIR} \
				--enable-debug \
				--enable-simd=GEN \
				--enable-comms=mpi \
				--with-lime=${TOPDIR}/deps/install/scalar \
				MPICXX=$MPICXX \
				CXXFLAGS="-std=c++17 -Wno-psabi"
			status=$?
			echo "Configure exit status $status"
			;;
		esac
		;;

	scalar)
		case ${SOURCE} in
		grid)
			${SRCDIR}/configure \
				--prefix=${INSTALLDIR} \
				--enable-debug \
				--enable-simd=GEN \
				--enable-comms=none \
				--enable-unified=no \
				--enable-shm=none \
				--enable-reduction=mpi \
				--with-lime=${TOPDIR}/deps/install/scalar \
				CXXFLAGS="-std=c++17 -Wno-psabi"
			status=$?
			echo "Configure exit status $status"
			;;
		hadrons)
			../configure \
				--prefix ${INSTALLDIR} \
				--with-grid=${TOPDIR}/Grid/install-${ARCH}
			status=$?
			echo "Configure exit status $status"
			;;
		app)
			../configure \
				--prefix ${INSTALLDIR} \
				--with-grid=${TOPDIR}/Grid/install-${ARCH} \
				--with-hadrons=${TOPDIR}/Hadrons/install-${ARCH}
			status=$?
			echo "Configure exit status $status"
			;;
		esac
		;;
	esac

	if [ $status -ne 0 ]; then
		echo "Quitting because of configure errors"
	else
		echo "Building in ${BUILDDIR}"
		${MAKE} -k -j4

		echo "Installing in ${INSTALLDIR}"
		${MAKE} install
	fi

fi
