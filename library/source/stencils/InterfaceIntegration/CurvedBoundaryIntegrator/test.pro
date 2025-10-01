######################################################################
# Von C.Pflaum
######################################################################

CONFIG-=windows
CONFIG+=console
CONFIG -= moc

QMAKE_CXX = ccache g++

DEFINES += _MYASSERT




unix:QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-parameter

#CONFIG += openmp

TARGET = run

CONFIG += debug

# Input
SOURCES += mainPoissonCurved.cc      # for Parikshit Upadhyaya 


# headers
HEADERS += ../interfaceMatrices.h  ../constantIntegrators.h  ../interatorBasisFunction.h
           
           





