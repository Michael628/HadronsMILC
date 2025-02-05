AM_CXXFLAGS += -I$(top_srcdir)

lib_LIBRARIES = libHadronsMILC.a

include modules.inc

libHadronsMILC_a_SOURCES = \
	$(modules_cpp)

libHadronsMILC_adir = $(includedir)/HadronsMILC
nobase_libHadronsMILC_a_HEADERS = \
	A2AVectors.hpp            \
	A2AMatrix.hpp             \
	DilutedNoise.hpp          \
	Modules.hpp               \
	$(modules_hpp)
