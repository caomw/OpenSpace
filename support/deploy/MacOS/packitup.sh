#!/bin/sh
# This script packs up the executables, libraries, and subsets of the
# data to make a transportable package for OpenSpace.
#
# <insert LICENSE info here>
######################################################################

echo "* MacOS/packitup.sh invoked as $0"

##
# Start in $WORKSPACE

if [ -z "$WORKSPACE" ]; then
  HERE=`dirname $0`
  WORKSPACE=`cd $HERE/../../.. ; echo $PWD`
fi

if [ ! -d $WORKSPACE ]; then
 echo "WORKSPACE does not exist: $WORKSPACE"
 exit 1
fi

##
# Get version numbers for archive file name

cd $WORKSPACE/include/openspace
Major=`grep MAJOR openspace.h | awk '{print $5}' | tr -d ';'`
minor=`grep MINOR openspace.h | awk '{print $5}' | tr -d ';'`
patch=`grep PATCH openspace.h | awk '{print $5}' | tr -d ';'`
if [ -z "$Major$minor$patch" ]; then
  echo "Cannot determine version number. "
  exit 2
fi  

tarfile="OpenSpace-${Major}-${minor}-${patch}.tgz"

##
# Create deploy/OpenSpace to hold everything:

cd $WORKSPACE
mkdir -p deploy/OpenSpace

cd deploy/OpenSpace

##
# Copy top level files:

cp -p $WORKSPACE/LICENSE .
cp -p $WORKSPACE/CREDITS .
cp -p $WORKSPACE/openspace.cfg .
cp -p $WORKSPACE/Doxyfile  .

mkdir -p bin lib data

##
# Copy Debug executables, libraries, and scripts:

cp -p  $WORKSPACE/bin/openspace/Debug/openspace      bin
cp -p  $WORKSPACE/bin/openspace/Debug/openspacetest  bin
cp -Rp $WORKSPACE/bin/openspace/Debug/Launcher.app   bin

cp -p  $WORKSPACE/bin/lib/Debug/lib*  lib

cp -Rp $WORKSPACE/scripts .


##
# Copy basic data:

cd data
mkdir -p spice scene
cp -Rp $WORKSPACE/data/fonts  .

cd spice
cp $WORKSPACE/data/spice/*.bsp .
cd ..

cd scene
cp -p $WORKSPACE/data/scene/*.scene .
cp -Rp $WORKSPACE/data/scene/common .
cp -Rp $WORKSPACE/data/scene/earth .
cd ..


##
# Now pack it up in a .gz archive

cd $WORKSPACE/deploy
echo "* Creating compressed tar archive..."
tar czvf $tarfile  OpenSpace
