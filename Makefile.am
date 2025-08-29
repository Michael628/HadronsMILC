SUBDIRS = src utilities

AM_LDFLAGS += -L$(top_builddir)/src
AM_CXXFLAGS += -I$(top_srcdir)/src

bin_PROGRAMS = HadronsMILC

HadronsMILC_SOURCES = main.cpp
HadronsMILC_LDADD = -lHadronsMILC
