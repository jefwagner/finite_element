###########################################################
# General Purpose makefile:
###########################################################
# Author: Jef Wagner
# Date: 2016-06-20
#
###########################################################

# The program name
PROGRAM = finite_element

# Get a list of source files, filtering out unit test
SRCS = $(filter-out test_%.cpp, $(wildcard *.cpp))
# Get a list of object files from the source files
OBJS = $(patsubst %.cpp, %.o, $(SRCS))
# Get a list of dependencies from the source files
DEPS = $(patsubst %.cpp, %.dep, $(SRCS))

# Choose our C and C++ compilers
CC = gcc
CXX = g++
PY = python

# Choose the compiler flags to set
CFlAGS = -Wall -g
CFLAGS += -O0 -ggdb
CFLAGS += -fopenmp

CXXFLAGS = -Wall -g
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -fopenmp

# Set the inclusion and library path

# To make this compatible between OSX and Linux,

OS := $(Shell uname)

ifeq (OS, Darwin)
INCFLAGS = -I/opt/local/include
INCFLAGS += -I/opt/local/include/eigen3

LDFLAGS = -L/opt/local/lib

else
INCFLAGS = -I/usr/include
INCFLAGS += -I/usr/include/eigen3
INCFLAGS += -I/usr/local/include

LDFLAGS = -L/usr/lib
LDFLAGS += -L/usr/local/lib

endif

# Set the libraries to link to
LIBS = -ltriangle
LIBS += -lm
LIBS += -fopenmp

# Declare the test and clean as phoney
# This lets Make know that there should
# not be any files called test or clean
.PHONEY: test clean

# Build the program.
# This is the default behavior if no targets are declared
$(PROGRAM): $(OBJS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

# Include the dependencies
# This has all of the additional dependencies
# for building the source files to object files.
# See the comment below for more detail.
ifneq ($(strip $(DEPS)),)
include $(DEPS)
endif

# Implicit rule for automatic dependency generation
# the -MM command has the compiler simply output
# a list of header files in the Make format. For example
# If the source file 'mesh.cpp' starts:
# > #include 'mesh.hpp'
# > #include 'utils.hpp'
# then this will create a file called 'mesh.dep' which
# reads:
# > mesh.o: mesh.cpp mesh.hpp utils.hpp
%.dep: %.cpp
	@set -e; rm -f $@; \
	 $(CXX) -MM $(CXXFLAGS) $< > $@

# Implicit rule generating object files
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCFLAGS) $< -o $@

# If we call the test variable, run the unit test python file
test:
	$(PY) unit_test_cpp.py

clean:
	rm -f *.o *.dep $(PROGRAM)
