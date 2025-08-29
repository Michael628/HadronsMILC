AM_CXXFLAGS += -I$(top_srcdir)

lib_LIBRARIES = libHadronsMILC.a

include modules.inc

libHadronsMILC_a_SOURCES = \
	$(modules_cpp)

libHadronsMILC_adir = $(includedir)
nobase_libHadronsMILC_a_HEADERS = \
	EigenPack.hpp            \
	A2AVectors.hpp            \
	A2AMatrix.hpp             \
	DilutedNoise.hpp          \
	Modules.hpp               \
	a2a/A2ATask.h             \
	a2a/A2AView.h             \
	a2a/A2AWorker.h           \
	a2a/A2AWorker_old.h       \
	spin/StagGamma.h          \
	$(modules_hpp)
