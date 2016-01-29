#!/bin/bash
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`

clear; clear
rm -rf MakeFiles/ CMakeCache.txt Makefile cmake_install.cmake icarus CTestTestfile.cmake Testing/ CMakeFiles
rm -rf makelog testlog
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "> all files and directories that have to do with cmake/make have been deleted"
echo "> the following files and directories are left:"
ls -a
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ">>> running cmake <<<"
cmake CMakeLists.txt
echo ">>> done <<<"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ">>> running make <<<"
echo "(in brackets: number of [errors, warnings] caused by compiling this file)"
echo -ne "[ 0%]"
make 1> makelog 2>&1 &
pid=$!
oldstr="old"; numwarnings=0; numerrors=0
while kill -0 $pid 2> /dev/null; do
  newstr=$(grep -s "%]" makelog | tail -1)
  if [[ ("$oldstr" != "$newstr") && (${#newstr} -gt 10) ]] ; then
    echo " ["$(($(grep -c "error:" makelog) - numerrors))", "$(($(grep -c "warning:" makelog) - numwarnings))"]"
    echo -ne $newstr
    oldstr=$newstr
    numerrors=$(grep -c "error:" makelog)
    numwarnings=$(grep -c "warning:" makelog)
  fi
done
echo " ["$(($(grep -c "error:" makelog) - numerrors))", "$(($(grep -c "warning:" makelog) - numwarnings))"]"
echo ">>> done <<<"
if [[ $(grep -c "error:" makelog) != 0 ]] ; then
  echo -ne "${red}"
else
  echo -ne "${green}"
fi
echo "> total errors: "$(grep -c "error:" makelog)
if [[ $(grep -c "warning:" makelog) != 0 ]] ; then
  echo -ne "${yellow}"
else
  echo -ne "${green}"
fi
echo "> total warnings: "$(grep -c "warning:" makelog)"${reset}"
echo "> errors and warnings from runnig make have been put into the file 'makelog'"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ">>> running ctest <<<"
ctest 1> testlog 2>&1 &
pid=$!
oldstr="old"
while kill -0 $pid 2> /dev/null; do
  newstr=$(grep "Test #" testlog | tail -1)
  if [[ ("$oldstr" != "$newstr") && ($newstr == *"sec"*) ]] ; then
    echo $newstr
    oldstr=$newstr
  fi
done
echo ">>> done <<<"
if [[ $(grep "tests passed" testlog) != *"100%"* ]] ; then
  echo -ne "${red}"
else
  echo -ne "${green}"
fi
echo "> "$(grep "tests passed" testlog)${reset}
echo "> "$(grep "Total Test time" testlog)
echo "> errors and warnings from runnig ctest have been put into the file 'testlog'"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
