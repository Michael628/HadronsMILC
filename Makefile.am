SUBDIRS = src/HadronsMILC utilities

AM_LDFLAGS += -L$(top_builddir)/src/HadronsMILC
AM_CXXFLAGS += -I$(top_srcdir)/src/HadronsMILC

bin_PROGRAMS = HadronsMILC

HadronsMILC_SOURCES = main.cpp
HadronsMILC_LDADD = -lHadronsMILC
