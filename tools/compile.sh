#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

mkdir out
sh -f configure.sh Emscripten-clang
cd build/Emscripten-clang-Release
cmake --build .
cd bin
cp ../../../data/pig.obj . 
sh ../../../tools/gen_emscripten_html.sh contours_viewer.js pig.obj
rm pig.obj
