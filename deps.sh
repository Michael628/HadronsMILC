#!/bin/bash
# Installs gmp, mpfr, lime, and (sometimes) hdf5
gmp_flag='false'
mpfr_flag='false'
lime_flag='false'
hdf5_flag='false'

ld_flag='-nostartfiles'

while test $# -gt 0; do
	echo "param: $1"
	case "$1" in
	--gmp)
		gmp_flag='true'
		shift
		;;
	--mpfr)
		mpfr_flag='true'
		shift
		;;
	--lime)
		lime_flag='true'
		shift
		;;
	--hdf5)
		hdf5_flag='true'
		shift
		;;
	--all)
		gmp_flag='true'
		mpfr_flag='true'
		lime_flag='true'
		hdf5_flag='true'
		shift
		;;
	scalar)
		ARCH=$1
		shift
		;;
	*)
		break
		;;
	esac
done

TOPDIR=$(pwd)
WORKDIR=${TOPDIR}/deps
INSTALLDIR=${WORKDIR}/install/${ARCH}

mkdir -p ${WORKDIR}
pushd ${WORKDIR}

if [ $gmp_flag = 'true' ]; then
	if [ ! -d gmp-6.2.1 ]; then
		wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz
		tar -xf gmp-6.2.1.tar.xz
	fi

	pushd gmp-6.2.1
	rm -rf build-${ARCH}
	mkdir -p build-${ARCH}
	pushd build-${ARCH}
	../configure \
		--prefix=${INSTALLDIR} \
		CXXFLAGS=${cxx_flags} CFLAGS=${cxx_flags} \
		LDFLAGS=${ld_flags}
	make all install
	status=$?

	popd
	if [[ $status -ne 0 ]]; then
		echo "gmp compile failed."
		exit 1
	fi
	popd
fi

if [ $mpfr_flag = 'true' ]; then
	if [ ! -d mpfr-4.1.0 ]; then
		wget https://www.mpfr.org/mpfr-4.1.0/mpfr-4.1.0.tar.gz
		tar -xvzf mpfr-4.1.0.tar.gz
	fi

	pushd mpfr-4.1.0
	rm -rf build-${ARCH}
	mkdir -p build-${ARCH}
	pushd build-${ARCH}
	../configure \
		--prefix=${INSTALLDIR} \
		--with-gmp=${INSTALLDIR} \
		CXXFLAGS=${cxx_flags} CFLAGS=${cxx_flags}
	make all install
	status=$?
	popd

	if [ $status -ne 0 ]; then
		echo "mpfr compile failed."
		exit 1
	fi
	popd
fi

if [ $lime_flag = 'true' ]; then
	if [ ! -d lime-1.3.2 ]; then
		wget http://usqcd-software.github.io/downloads/c-lime/lime-1.3.2.tar.gz
		tar -xvzf lime-1.3.2.tar.gz
	fi

	pushd lime-1.3.2
	rm -rf build-${ARCH}
	mkdir -p build-${ARCH}
	pushd build-${ARCH}
	../configure \
		--prefix=${INSTALLDIR} \
		CXXFLAGS=${cxx_flags} CFLAGS=${cxx_flags}
	make all install
	status=$?
	popd

	if [ $status -ne 0 ]; then
		echo "lime compile failed."
		exit 1
	fi
	popd
fi

if [ $hdf5_flag = 'true' ]; then
	if [ ! -d hdf5-1.10.10 ]; then
		#    wget https://www.hdfgroup.org/package/hdf5-1-10-7-tar-gz/?wpdmdl=15050&refresh=6329defe5b0d91663688446
		wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.10/src/hdf5-1.10.10.tar.gz
		tar -xvzf hdf5-1.10.10.tar.gz
	fi

	pushd hdf5-1.10.10
	rm -rf build-${ARCH}
	mkdir -p build-${ARCH}
	pushd build-${ARCH}
	../configure \
		--prefix=${INSTALLDIR} \
		--enable-cxx \
		CXXFLAGS=${cxx_flags} CFLAGS=${cxx_flags}
	make all install
	popd
	if [ $status -ne 0 ]; then
		echo "hdf5 compile failed."
		exit 1
	fi
	popd
fi
