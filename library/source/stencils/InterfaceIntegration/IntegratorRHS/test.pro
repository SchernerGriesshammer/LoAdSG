######################################################################
# Von C.Pflaum
######################################################################

CONFIG -= windows
CONFIG += console
CONFIG -= moc
CONFIG -= gui
CONFIG += c++11

QMAKE_CXX = ccache g++

QMAKE_CXXFLAGS_RELEASE *= -O3

DEFINES += _MYASSERT




unix:QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-parameter

TARGET = test

#CONFIG += release
CONFIG += debug

# Input
SOURCES += mainVariableHelm.cc        # for Tim Reinfels
HEADERS += integral.hpp util.hpp

# headers
HEADERS += ../interfaceMatrices.h  ../constantIntegrators.h  ../interatorBasisFunction.h
           
           





