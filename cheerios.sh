#!/bin/bash

#Set File Directory

SCRIPT_PATH = $(dirname $(readlink -f $0))

#Set Directory

cd
cd SCRIPT_PATH

#Run Make, not Make test

make > log 2>&1
if ["$?" -ne 0 ]; then
  echo "Build Failed"
  cat log
fi

#Run file, Note: This does not run make test so DO NOT use this script for
# debugging of any sort.

./finite_element

#This is the mathematica analysis part, this might have issues since this
# assumes that mathematica is already in your path

math -noprompt -run "<<PlotGenerator.nb"
