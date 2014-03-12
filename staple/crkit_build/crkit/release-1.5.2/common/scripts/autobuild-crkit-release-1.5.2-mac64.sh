#!/bin/bash

# *************************************************************
# CRKit UNIX/LINUX/Mac OS X BUILD SCRIPT - READ CAREFULLY
# *************************************************************
# This script will download and build CRKit and all the
# supporting software. 
#
# In order to run this script, you must have installed:
# gcc
# gmake
# subversion (svn)
# curl
#
# Building this software will require around 3GB of disk space.
#
# This script downloads CMake, ITK, VTK, Qt and CRKit. 
#
# *************************************************************

# Uncomment the line below to see the commands as they are executed.
##### set -x -v

# -------------------------------------------------------------
# Options that the user can supply
# -------------------------------------------------------------

# The directory where the build will take place
basedir="`pwd`/crkit_build"

# The installation directory, defaults to $basedir/install; 
# The executables and libraries are built into an `install' directory,
# and then packaged into a package (.dmg, .tar.gz or NSIS installer)
# ready for suitable deployment on this system or any other.
#
# To test the build, check the CRKit executabless in $instdir/crkit/bin
instdir=$basedir/install

echo instdir is $instdir

# The location of the build log file
logfile=$basedir/buildlog

# To use already installed software, set:
# cmake in PATH
# ITK_DIR
# qmake in PATH
# QT_DIR
# VTK_DIR
#

unset VTK_DIR
unset QT_DIR
unset ITK_DIR
PATH=/bin:/usr/bin

# Versions for the supporting software we are expecting. 
CMAKE_VERSION=2.8.3
ITK_VERSION=3.20.0
VTK_MAJOR_VERSION=5.6
VTK_MINOR_VERSION=1
VTK_VERSION=${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}
# For Max OS X < 10.7
# QT_VERSION=4.7.4
# For Mac OS X 10.7 :
# QT_VERSION=4.8.0
QT_VERSION=4.7.4
CRKIT_VERSION=release-1.5.2
CRKIT_FORCE_REBUILD=true

if [ "x`uname`" == "xDarwin" ] ; then
  QTPLATFORM="mac-opensource"
  QTPLATFORM="everywhere-opensource"
else 
  QTPLATFORM="x11-opensource"
  QTPLATFORM="everywhere-opensource"
fi

# -------------------------------------------------------------
# Check for necessary software
# -------------------------------------------------------------
for prog in "g++" "gcc" "make" "svn" "curl" "gunzip"
do
  PROGBIN=`which $prog`
  if [ "x${PROGBIN}" == "x" ]
  then
    echo "Cannot find $prog on your system"
    echo "Please ask the administrator to install $prog and to "
    echo "make sure that $prog is in your path."
    exit
  fi
done

# -------------------------------------------------------------
# Initialization
# -------------------------------------------------------------
echo "Preparing to build CRKit supporting software."
echo " This will take 4 hours and 3000MB of disk space"
echo "Log messages are placed into $logfile"

function detect_cmake {
  if [ ! `which cmake` ] ; then
    echo "You may set a default cmake."
  else
    cmakebin=`which cmake`
  fi
  if [ ! -z ${cmakebin} -a -x ${cmakebin} ]
  then
    echo "Using system cmake : ${cmakebin}"
  else
    cmakebin=$instdir/cmake/bin/cmake
  fi
}

# -------------------------------------------------------------
# Check out and build cmake
# -------------------------------------------------------------
function get_cmake {
  if [ -x ${cmakebin} ] 
  then
    echo "Skipping rebuilding of cmake."
# Place the newly built CMake first in our PATH
    PATH=$instdir/cmake/bin:$PATH
    export PATH
    return 0
  fi

  mkdir -p $basedir/cmake
  cd $basedir/cmake

  # Use curl to download the right version
  curl -o cmake-${CMAKE_VERSION}.tar.gz http://crkit.s3.amazonaws.com/cmake-${CMAKE_VERSION}.tar.gz
  gunzip -c cmake-${CMAKE_VERSION}.tar.gz | tar xvf -
  if [ $? -ne 0 ]; then echo "Downloading CMake failed." ; exit; fi
  mv cmake-${CMAKE_VERSION} CMake

  # Build cmake
  echo "Building CMake"
  cd $basedir/cmake/CMake
  ./configure --prefix=$instdir/cmake >> $logfile
  make -j 2 >> $logfile
  make install >> $logfile

  # Test to see if everything is ok
  if [ ! -x ${cmakebin} ]
  then
    echo "Failed to create CMake executable"
    exit -1
  fi

# Place the newly built CMake first in our PATH
  PATH=$instdir/cmake/bin:$PATH
  export PATH
}


# We need to use a VTK that has QVTK configured and installed.
function detect_vtk {
  if [ -f ${VTK_DIR}/UseVTK.cmake ] 
  then
    qvtk=`ls ${VTK_DIR}/libQVTK*.*`
    # Check if QVTK is compiled with this VTK.
    if [ "x${qvtk}" == "x" ] ; then
      echo "Did not detect QVTK widget with ${VTK_DIR}."
      echo "Building a new VTK to use QVTK."
      VTK_DIR=$instdir/vtk/lib/vtk-${VTK_MAJOR_VERSION}
    fi
    echo "Using system VTK : ${VTK_DIR}"
  else 
    VTK_DIR=$instdir/vtk/lib/vtk-${VTK_MAJOR_VERSION}
    echo "Set VTK_DIR to : ${VTK_DIR}"
  fi
}

# -------------------------------------------------------------
# Check out and build VTK
# -------------------------------------------------------------
function get_vtk {
  if [ -f ${VTK_DIR}/UseVTK.cmake ] 
  then
    echo "Skipping rebuilding of VTK."
    return 0
  fi

  mkdir -p $basedir/vtk/bingcc
  cd $basedir/vtk

  # Use curl to download VTK.
  curl -o vtk-${VTK_VERSION}.tar.gz http://crkit.s3.amazonaws.com/vtk-${VTK_VERSION}.tar.gz
  gunzip -c vtk-${VTK_VERSION}.tar.gz | tar xvf -
  if [ $? -ne 0 ]; then echo "Downloading VTK failed." ; exit; fi

  # Configure VTK.
  echo "Building VTK"
  cd $basedir/vtk/bingcc
  ${cmakebin} \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DDESIRED_QT_VERSION:STRING=4 \
    -DVTK_USE_CARBON:BOOL=OFF \
    -DVTK_USE_COCOA:BOOL=ON \
    -DVTK_USE_GUISUPPORT:BOOL=ON \
    -DVTK_USE_64BIT_IDS:BOOL=ON \
    -DVTK_USE_QT:BOOL=ON \
    -DVTK_USE_QVTK:BOOL=ON \
    -DQT_QMAKE_EXECUTABLE:PATH="${qmakebin}" \
    -DVTK_USE_HYBRID:BOOL=ON \
    -DVTK_USE_ANSI_STDLIB:BOOL=ON \
    -DVTK_USE_PARALLEL:BOOL=ON \
    -DVTK_USE_RENDERING:BOOL=ON \
    -DVTK_USE_PATENTED:BOOL=ON \
    -DCMAKE_CXX_FLAGS_RELEASE:STRING="-O3 -DNDEBUG" \
    -DCMAKE_INSTALL_PREFIX:PATH=$instdir/vtk \
    $basedir/vtk/VTK >> $logfile
  make -j 2 >> $logfile
  make install >> $logfile

  if [ ! -f ${VTK_DIR}/UseVTK.cmake ] ; then
      echo "VTK library ${VTK_DIR}/UseVTK.cmake failed to build!"
      exit -1
  fi
}

function detect_itk {
  if [ -f ${ITK_DIR}/UseITK.cmake ] ; then
    echo "Using ITK from ${ITK_DIR}/UseITK.cmake"
  else
    ITK_DIR=${instdir}/itk/lib/InsightToolkit
  fi
}

# -------------------------------------------------------------
# Check out and build ITK
# -------------------------------------------------------------
function get_itk {
  if [ -f ${ITK_DIR}/UseITK.cmake ] 
  then
    echo "Skipping rebuilding of ITK."
    return 0
  fi

  echo "Now building ITK for ${ITK_DIR}"
  mkdir -p $basedir/itk/bingcc
  cd $basedir/itk

  # Use curl to download ITK.
  curl -o InsightToolkit-${ITK_VERSION}.tar.gz http://crkit.s3.amazonaws.com/InsightToolkit-${ITK_VERSION}.tar.gz
  gunzip -c InsightToolkit-${ITK_VERSION}.tar.gz | tar xvf -
  if [ $? -ne 0 ]; then echo "Downloading InsightToolkit failed." ; exit; fi
  mv InsightToolkit-${ITK_VERSION} Insight

  # Configure ITK using CMake
  echo "Building ITK"
  cd $basedir/itk/bingcc
  ${cmakebin} \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DITK_USE_PORTABLE_ROUND:BOOL=ON \
    -DITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_CXX_FLAGS_RELEASE:STRING="-O3 -DNDEBUG" \
    -DCMAKE_INSTALL_PREFIX:PATH=$instdir/itk \
    $basedir/itk/Insight >> $logfile
  make -j 2 >> $logfile
  make install >> $logfile

  # Check whether the necessary libraries built
  if [ ! -f ${ITK_DIR}/UseITK.cmake ] 
  then
    echo "ITK library $ITK_DIR/UseITK.cmake failed to build."
      exit -1
  fi
}

function detect_qt {
  qmakebin=`which qmake`
  if [ ! -z ${qmakebin} -a -x $qmakebin ]
  then
    echo "Using system qmake : $qmakebin"
  else
    qmakebin=$instdir/qt-${QT_VERSION}/bin/qmake
  fi
}

# -------------------------------------------------------------
# Check out and build Qt
# -------------------------------------------------------------
function get_qt {
  if [ -x $qmakebin ]
  then
    echo "Skipping rebuilding of Qt."
    QTBIN=`dirname $qmakebin`
    PATH=${QTBIN}:${PATH}
    QTDIR=`dirname ${QTBIN}`
    export PATH QTDIR
    return 0
  fi
  mkdir -p $basedir/qt
  cd $basedir/qt

  if [ ! -f qt-${QTPLATFORM}-src-${QT_VERSION}.tar.gz ] ; then
    curl http://crkit.s3.amazonaws.com/qt-${QTPLATFORM}-src-${QT_VERSION}.tar.gz -o qt-${QTPLATFORM}-src-${QT_VERSION}.tar.gz
  fi
  if [ $? -ne 0 ]; then echo "Downloading Qt failed." ; exit; fi
  if [ ! -d qt-${QTPLATFORM}-src-${QT_VERSION} ] ; then
    echo "There is no directory : qt-${QTPLATFORM}-src-${QT_VERSION}"
    gunzip -c qt-${QTPLATFORM}-src-${QT_VERSION}.tar.gz | tar xvf -
  fi

  # Configure Qt using CMake
  echo "Building Qt"
  cd $basedir/qt/qt-${QTPLATFORM}-src-${QT_VERSION}

  QTCONFARGS=""
#   Qt 4.5.2 or higher requires a new option: -opensource
  QTCONFARGS="${QTCONFARGS} -opensource"

  echo Now running ./configure $QTCONFARGS -confirm-license -release \
    -no-phonon -no-sql-odbc -prefix $instdir/qt-${QT_VERSION}
  ./configure $QTCONFARGS -confirm-license -release \
    -no-phonon -no-sql-odbc -prefix $instdir/qt-${QT_VERSION} >> $logfile
  make >> $logfile
  make install >> $logfile

  # Check whether the necessary libraries built
  for bin in "moc" "qmake" "uic"
  do
    if [ ! -e $instdir/qt-${QT_VERSION}/bin/${bin} ]
    then
      echo "Qt binary $instdir/qt-${QT_VERSION}/bin/${bin} failed to build!"
      exit -1
    fi
  done

# Ensure the newly built Qt is first in the PATH
  QTDIR=$instdir/qt-${QT_VERSION}
  PATH=${QTDIR}/bin:$PATH
  export PATH QTDIR
}

# Find the source code for CRKit
function detect_crkit {
  if [ -z ${CRKIT_SOURCE} ] ; then
    CRKIT_SOURCE=${basedir}/crkit/${CRKIT_VERSION}
    echo "Setting CRKIT_SOURCE to be ${CRKIT_SOURCE}"
  elif [ -d ${CRKIT_SOURCE} ] ; then
    echo "Using CRKit source code : ${CRKIT_SOURCE}"
  else
    echo "Did not find source in ${CRKIT_SOURCE}".
    CRKIT_SOURCE=${basedir}/crkit/${CRKIT_VERSION}
    echo "Setting CRKIT_SOURCE to be ${CRKIT_SOURCE}"
  fi
}

# -------------------------------------------------------------
# Check out and build CRKit
# -------------------------------------------------------------
function get_crkit {
  if [ "x${CRKIT_FORCE_REBUILD}" != "xtrue" -a -f $instdir/bin/crlSTAPLE ] ; then
    echo "Skipping rebuilding of CRKit."
    return 0
  fi

  mkdir -p $basedir/crkit/bingcc
  cd $basedir/crkit

  if [ ! -d ${CRKIT_SOURCE} ] ; then
    # Use curl to download CRKit.
    curl -o crkit-${CRKIT_VERSION}.tar.gz http://crkit.s3.amazonaws.com/crkit-${CRKIT_VERSION}.tar.gz
    gunzip -c crkit-${CRKIT_VERSION}.tar.gz | tar xvf -
    if [ $? -ne 0 ]; then echo "Downloading CRKit failed." ; exit; fi
    mv crkit-${CRKIT_VERSION} ${CRKIT_VERSION}
    echo "Building from downloaded source code in ${CRKIT_SOURCE}."
  else
    echo "Going to build from CRKit code in ${CRKIT_SOURCE}."
  fi

  # Configure CRKit using CMake
  echo "Building CRKit"
  cd $basedir/crkit/bingcc
  ${cmakebin} \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_CXX_FLAGS_RELEASE:STRING="-O3 -DNDEBUG" \
    -DUSE_ITK:BOOL=ON \
    -DUSE_VTK:BOOL=ON \
    -DUSE_QT:BOOL=ON \
    -DDESIRED_QT_VERSION:STRING=4 \
    -DQT_QMAKE_EXECUTABLE:PATH="${qmakebin}" \
    -DITK_DIR:PATH="${ITK_DIR}" \
    -DVTK_DIR:PATH="${VTK_DIR}" \
    -DCRKIT_INSTALL_PREFIX:PATH=$instdir \
    -DCMAKE_INSTALL_PREFIX:PATH=$instdir \
    ${CRKIT_SOURCE} >> $logfile

  # Make only in the CRKit directory
  make 2>&1 >> $logfile
  make install >> $logfile

  if [ ! -f $instdir/bin/crlSTAPLE ]
  then
    echo "CRKit failed to build in $instdir/bin/ "
    exit -1;
  fi
}

# -------------------------------------------------------------
# Perform the actual build tasks
# -------------------------------------------------------------

detect_cmake
get_cmake

detect_itk
get_itk

detect_qt
get_qt

detect_vtk
get_vtk

detect_crkit
get_crkit

echo "CRKit build and installation completed."
echo "CRKit package is located in $instdir/bin/"

cd $basedir/crkit/bingcc
make package

# If we are running on a Mac, then we need to take the following steps to make
# an installable DMG file.
#
echo "Complete the DMG building process."
# First copy the fixed up MacOS content - dylibs.
cp $instdir/MacOS/* _CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin/CRKit.app/Contents/MacOS/
# Second, copy the fixed up frameworks and libraries for Qt:
mkdir _CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin/CRKit.app/Contents/Frameworks
cp -r $instdir/FrameWorks/Qt*.framework _CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin/CRKit.app/Contents/Frameworks

# Third, copy the fixed up binaries:
# The fixed up binaries expect to see ../MacOS and ../Frameworks.
# As a result, we copy them into MacOS
cp $instdir/bin/* _CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin/CRKit.app/Contents/MacOS

# Fourth: Remove the binaries that have been put into the bundle Resource,
# as these are not configured correctly.
rm _CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin/CRKit.app/Contents/Resources/bin/*
 
# Fifth, remake the DMG file:
rm _CPack_Packages/Darwin/Bundle/*.dmg
hdiutil create -srcdir _CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin _CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin.dmg
echo "DMG file is in $basedir/crkit/bingcc/_CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin called CRKIT-1.5.2-Darwin.dmg"
cp $basedir/crkit/bingcc/_CPack_Packages/Darwin/Bundle/CRKIT-1.5.2-Darwin.dmg $basedir/crkit/bingcc

# Sixth, remove the install libraries from the build tree.
mv $instdir ${instdir}.old

exit 0

