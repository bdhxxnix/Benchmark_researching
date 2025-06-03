# B+ Tree Test Cases
# $Id: Makefile.am 130 2011-05-18 08:24:25Z tb $

noinst_PROGRAMS = testsuite

TESTS = testsuite

if HAVE_CPPUNIT

testsuite_SOURCES = TestRunner.cc

testsuite_SOURCES += InstantiationTest.cc SimpleTest.cc LargeTest.cc BoundTest.cc IteratorTest.cc StructureTest.cc DumpRestoreTest.cc RelationTest.cc

else

testsuite_SOURCES = TestTrue.cc

endif

AM_CXXFLAGS = -W -Wall -Wold-style-cast -Wshadow -DBTREE_DEBUG -I$(top_srcdir)/include @CPPUNIT_CFLAGS@
LDADD = @CPPUNIT_LIBS@
