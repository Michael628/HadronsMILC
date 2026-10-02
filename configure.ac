AC_PREREQ([2.69])
AC_INIT([HadronsMILC], [0.2.1], [michaellynch628@gmail.com])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_ARG_WITH([grid],
    [AS_HELP_STRING([--with-grid=<prefix>],
    [try this for a non-standard install prefix of Grid])],
    [PATH="$with_grid/bin$PATH_SEPARATOR$PATH"]
    [CXXFLAGS="$CXXFLAGS -I$with_grid/include"]
    [LDFLAGS="$LDFLAGS -L$with_grid/lib"])
AC_ARG_WITH([hadrons],
    [AS_HELP_STRING([--with-hadrons=<prefix>],
    [try this for a non-standard install prefix of Hadrons])],
    [PATH="$with_hadrons/bin$PATH_SEPARATOR$PATH"]
    [CXXFLAGS="$CXXFLAGS -I$with_hadrons/include"]
    [LDFLAGS="$LDFLAGS -L$with_hadrons/lib"])
AC_ARG_ENABLE([strict-pins],
    [AS_HELP_STRING([--enable-strict-pins],
    [fail when Grid/Hadrons commits differ from the DEPENDENCIES pins])],
    [ac_STRICT_PINS=${enable_strict_pins}], [ac_STRICT_PINS=no])
AC_CHECK_PROG([HADCONF],[hadrons-config],[yes])
if test x"$HADCONF" != x"yes" ; then
    AC_MSG_ERROR([hadrons-config not found])
fi
if test x"$CXX" == x ; then
    CXX="`hadrons-config --cxx`"
elif test "$CXX" != "`hadrons-config --cxx`" ; then
    AC_MSG_WARN([CXX differs from that reported by hadrons-config])
fi
if test x"$CXXLD" == x ; then
    CXXLD="`hadrons-config --cxxld`"
elif test "$CXXLD" != "`hadrons-config --cxxld`" ; then
    AC_MSG_WARN([CXXLD differs from that reported by hadrons-config])
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_AR
AC_LANG([C++])

CXXFLAGS="$CXXFLAGS `hadrons-config --cxxflags`"
LDFLAGS="$LDFLAGS `hadrons-config --ldflags`"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
LIBS=" -lHadrons -lGridMilc $LIBS `hadrons-config --libs`"

AC_MSG_CHECKING([that a minimal Grid/Hadrons/GridMilc program compiles]);
AC_LINK_IFELSE(
        [AC_LANG_SOURCE([[
    #include <Grid/Grid.h>
    #include <Hadrons/Application.hpp>
    #include <GridMilc/GridMilc.h>

    using namespace Grid;
    using namespace Hadrons;

    int main(int argc, char *argv[])
    {
        Grid_init(&argc, &argv);
        Grid_finalize();
        Application application;
        application.run();

        return 0;
    }

    ]])],
        [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([Could not compile a minimal Hadrons/GridMilc program])])

# Version identities and dependency pins ######################################
HMILC_GIT="`git -C "$srcdir" describe --tags --always --dirty 2>/dev/null`"
if test x"$HMILC_GIT" = x ; then
    HMILC_GIT=unknown
fi
AC_DEFINE_UNQUOTED([HADRONSMILC_GIT], ["$HMILC_GIT"],
    [HadronsMILC git describe at configure time])

HADRONS_GIT="`hadrons-config --git 2>/dev/null`"
if test x"$HADRONS_GIT" = x ; then
    HADRONS_GIT=unknown
fi
HADRONS_SHA_ACTUAL="`echo "$HADRONS_GIT" | awk '{ print $2 }'`"
AC_DEFINE_UNQUOTED([HADRONSMILC_HADRONS_GIT], ["$HADRONS_GIT"],
    [Hadrons branch and commit from hadrons-config --git at configure time])

AC_CHECK_PROG([GRIDCONF], [grid-config], [yes])
GRID_VERSION_H=""
if test x"$GRIDCONF" = x"yes" ; then
    GRID_VERSION_H="`grid-config --prefix`/include/Grid/Version.h"
fi
if test ! -f "$GRID_VERSION_H" && test x"$with_grid" != x ; then
    GRID_VERSION_H="$with_grid/include/Grid/Version.h"
fi
GRID_SHA_ACTUAL=""
if test -f "$GRID_VERSION_H" ; then
    GRID_SHA_ACTUAL="`sed -n 's/^#define GITHASH "\(@<:@0-9a-f@:>@\{40\}\).*/\1/p' "$GRID_VERSION_H"`"
    if grep -q 'uncommited changes' "$GRID_VERSION_H" ; then
        AC_MSG_WARN([Grid was built from a tree with uncommitted changes ($GRID_VERSION_H)])
    fi
else
    AC_MSG_WARN([Grid/Version.h not found, cannot determine the Grid commit])
fi

PIN_FILE="$srcdir/DEPENDENCIES"
PIN_STATUS=unchecked
AC_MSG_CHECKING([for dependency pins in DEPENDENCIES])
if test -f "$PIN_FILE" ; then
    AC_MSG_RESULT([yes])
    PIN_STATUS=match
    for dep in grid hadrons ; do
        case $dep in
            grid)    actual="$GRID_SHA_ACTUAL" ;;
            hadrons) actual="$HADRONS_SHA_ACTUAL" ;;
        esac
        pin="`awk -v n="$dep" '$1 == n { print $3 }' "$PIN_FILE"`"
        if test x"$pin" = x ; then
            AC_MSG_WARN([no $dep pin in DEPENDENCIES])
        elif test x"$actual" = x ; then
            AC_MSG_WARN([cannot determine the $dep commit (DEPENDENCIES pin $pin)])
            PIN_STATUS=mismatch
        elif test "$actual" != "$pin" ; then
            AC_MSG_WARN([$dep commit $actual differs from DEPENDENCIES pin $pin])
            PIN_STATUS=mismatch
        else
            AC_MSG_NOTICE([$dep commit matches DEPENDENCIES pin $pin])
        fi
    done
    if test x"$PIN_STATUS" = x"mismatch" && test x"$ac_STRICT_PINS" = x"yes" ; then
        AC_MSG_ERROR([dependency commits differ from DEPENDENCIES pins (--enable-strict-pins)])
    fi
else
    AC_MSG_RESULT([no, skipping pin check])
fi
AC_DEFINE_UNQUOTED([HADRONSMILC_PINS], ["$PIN_STATUS"],
    [Dependency pin check result: match, mismatch or unchecked])

AC_SUBST([CXXLD])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([utilities/Makefile])
AC_OUTPUT
