#!/bin/sh
echo "Install Gauss on OSX"
echo "Download Qt"
brew update
brew install cmake
brew install eigen
brew upgrade cmake
brew upgrade eigen
wget http://download.qt.io/official_releases/qt/5.9/5.9.0/qt-opensource-mac-x64-5.9.0.dmg

echo "Enter password to install Qt"
sudo hdiutil attach -noverify ./qt-opensource-mac-x64-5.9.0.dmg
/Volumes/qt-opensource-mac-x64-5.9.0/qt-opensource-mac-x64-5.9.0.app/Contents/MacOS/qt-opensource-mac-x64-5.9.0 --platform minimal --script qt_script_osx.qs

mkdir ./build
cd ./build


echo "Building make files (this can fail if homebrew CMake version is wrong)"
cmake .. -DCMAKE_PREFIX_PATH=~/Qt5.9.0/5.9/clang_64/lib/cmake -DEigen4_DIR=/usr/local/Cellar/eigen/3.3.4/include/eigen3 -DCMAKE_BUILD_TYPE=Release

make -j 2 all
