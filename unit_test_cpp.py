#######################################################################
# Unit Testing for c++ projects
#######################################################################
# This file is used to implement a simple unit testing framework
# for c++ projects.
#
# The main Makefile for the project compiles the source files, but
# should ignore files that start with test_.
#
# Any file that starts with test_ should contain at least one test
# function, which has the signature:
#  > void test_{NAME}();
# The test file should call at the end the function print_status,
# which takes a boolean (or int) and a string {FUNC_NAME}, and either
# prints:
#  - {FUNC_NAME} failed.
#  + {FUNC_NAME} passed!
# The print_status function has the following signature:
#  > void print_status( int status, string func_name);
#
# This file will:
# 1. Generate all neccessary dependencies to compile all test_*.cpp
#    files
# 2. Compile all test_*.cpp and other neccessary files to object files
# 3. Create a test_main.cpp that calls all test_*() functions, then
#    Compile and link the program.
#######################################################################
from __future__ import print_function

import sys
from glob import glob
import os
from os.path import isfile, getmtime
import re
import shlex
from subprocess import Popen, PIPE

#
# Clean all the files generated for testing
#
def clean_test_files():
    test_files = glob('test_*.o')
    for file in test_files:
        os.remove(file)
    if isfile('test__main.cpp'):
        os.remove('test__main.cpp')
    if isfile('test.test'):
        os.remove('test.test')

#
# Load, and read the makefile to grab the compiler flags
#
def read_compiler_flags():
    """Read the compiler flags from the makefile"""
    f = open('Makefile','r')
    flags = {
    'CC' : '',
    'CXX' : '',
    'CFLAGS' : '',
    'CXXFLAGS' : '',
    'INCFLAGS' : '',
    'LDFLAGS' : '',
    'LIBS' : ''}
    for line in f:
        words = [ word.strip() for word in line.split('+=')]
        if len(words) == 1:
            words = [ word.strip() for word in line.split('=')]
        if len(words) > 1:
            if words[0] in flags:
                flags[words[0]] += (' ' + words[1])
    f.close()
    return flags

#
# Get the build dependencies for a file
#
def get_build_deps(file):
    """Use the -MM compiler flag to get a list of dependencies"""
    cmd_string = '{} -MM {}'.format(
        flags['CXX'],file)
    proc = Popen(shlex.split(cmd_string),stdout=PIPE, stderr=PIPE)
    (stdoutdata,stderrdata) = proc.communicate()
    rule = stdoutdata.decode('utf-8')
    # Process the rule to get target and dependencies
    target, deps = rule.strip().split(':')
    deps = deps.split()
    return (target, deps)

#
# Compile a source file to an object file
#
def compile_file(flags, file, target):
    print()
    print('Compiling {}'.format(file))
    cmd_string = '{} -c {} {} {} -o {}'.format(
        flags['CXX'],flags['CXXFLAGS'],flags['INCFLAGS'],file,target)
    print(cmd_string.strip())
    proc = Popen(shlex.split(cmd_string), stdout=PIPE, stderr=PIPE)
    (stdoutdata,stderrdata) = proc.communicate()
    if proc.returncode != 0:
        print('Error in compiling {}:'.format(file))
        print()
        print(stdoutdata.decode('utf-8'))
        print(stderrdata.decode('utf-8'))
        clean_test_files()
        sys.exit()
    else:
        print('Succesfully complied {}:'.format(file))
        print()

#
# Link all the files togeter
#
def link_files(flags, objects, target):
    print()
    print('Linking {}'.format(target))
    cmd_string = '{} {} {} {} -o {}'.format(
        flags['CXX'],flags['LDFLAGS'],flags['LIBS'],' '.join(objects),target)
    print(cmd_string.strip())
    proc = Popen(shlex.split(cmd_string), stdout=PIPE, stderr=PIPE)
    (stdoutdata,stderrdata) = proc.communicate()
    if proc.returncode != 0:
        print('Error in linking {}:'.format(target))
        print()
        print(stdoutdata.decode('utf-8'))
        print(stderrdata.decode('utf-8'))
        clean_test_files()
        sys.exit()
    else:
        print('Succesfully linked {}:'.format(target))
        print()

#
# Link all the files togeter
#
def run_test():
    cmd_string = './test.test'
    proc = Popen(shlex.split(cmd_string), stdout=PIPE, stderr=PIPE)
    (stdoutdata,stderrdata) = proc.communicate()
    if proc.returncode != 0:
        print('Error running ./test.test')
        print()
        print(stdoutdata.decode('utf-8'))
        print(stderrdata.decode('utf-8'))
#        clean_test_files()
        sys.exit()
    else:
        print(stdoutdata.decode('utf-8'))

#
# The source code that will give the main__test.cpp file
#
main_source = """
#include <iostream>

// prototypes
{}

int main(){}
    std::cout << std::endl;
    std::cout << "Testing utilility functions in utils.hpp" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
{}
    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::endl;
    return 0;
{}
"""

#######################################################################
# Start of the unit testing logic
#######################################################################
# Clean the system to start
clean_test_files()
# Get the compiler flags from the makefile
flags = read_compiler_flags()
# Start an empty list of non-test source files
extra_srcs = []
# Get a list of all test files
test_files =  glob('test_*.cpp')
# Start and empty list of test functions
test_functions = []
for file in test_files:
    # For each file get a list of dependecies in a Makefile rule using
    # the -MM compiler flag
    target, deps = get_build_deps(file)
    # Go through the dependecies to see if this will depend on any
    # non-test source and object files
    for dep in deps:
        if dep.split('.')[-1] == 'hpp':
            src_name = dep.split('.')[0]+'.cpp'
            if isfile(src_name):
                extra_srcs.append(src_name)
    # Go through the files and find the test_functions
    f = open(file,'r')
    for line in f:
        for word in line.split():
            if word.split('_')[0] == 'test':
                test_functions.append(re.sub('[(){};]','',word))
    f.close()
    # Compile the test file
    compile_file( flags, file, target)
# Optionally compile the other source files if the dependencies
# are newer than the object file
for src in extra_srcs:
    recompile = True
    obj, deps = get_build_deps(src)
#    print( 'src : {}'.format(src))
    dep_times = [getmtime(dep) for dep in deps]
    if isfile(obj):
        obj_time = getmtime(obj)
        if obj_time > max(dep_times):
            recompile = False
    if recompile:
        compile_file(flags, src, obj)
# Create the test__main.cpp file
prototypes = ['void '+func+'();' for func in test_functions]
prototypes = '\n'.join(prototypes)
function_calls = ['    '+func+'();' for func in test_functions]
function_calls = '\n'.join(function_calls)
f = open('test__main.cpp','w')
f.write(main_source.format(prototypes, '{', function_calls, '}'))
f.close()
# Compile the main test file!
compile_file(flags, 'test__main.cpp', 'test__main.o')
# Link the test files
test_objs = [re.sub('.cpp','.o',src) for src in test_files]
extra_objs = [re.sub('.cpp','.o',src) for src in extra_srcs]
all_objs = test_objs+extra_objs
all_objs.append('test__main.o')
link_files(flags, all_objs, 'test.test')
# Run the test
run_test()
print()
print('Successfully compiled and ran test.')
clean_test_files()
