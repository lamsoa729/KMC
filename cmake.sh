#!/bin/bash

do_build() {
  mkdir build
  cd build
  cmake ..
  make -j8
  cd ..
}

do_graph_build() {
  mkdir build
  cd build
  cmake -DGRAPH=1 ..
  make -j8
  cd ..
}

do_test_build() {
  mkdir build
  cd build
  cmake -DTESTS=1 ..
  make -j8
  if (make test)
  then
    cd ..
  else 
    ./tests/test_simcore
    cd ..
  fi
}

do_debug_build() {
  mkdir build
  cd build
  cmake -DDEBUG=1 ..
  make -j8
  cd ..
}

do_debug_graph_build() {
  mkdir build
  cd build
  cmake -DDEBUG=1 -DGRAPH=1 ..
  make -j8
  cd ..
}

do_docs_build() {
  mkdir build
  cd build
  cmake ..
  make docs
  cd ..
}

do_clean() {
  rm -rf build/CMake*
  rm -rf build/Makefile
  rm -rf build/cmake_install.cmake
  rm -rf build/src
  rm -rf build/Doxyfile
  rm -rf build/tests
}

do_usage() {
  echo "Usage: $0 [arg]"
  echo "arg must be one of:"
  echo "  clean   - remove temporary installation files"
  echo "  build   - build simcore without graphics"
  echo "  gbuild   - build simcore with graphics"
  echo "  debug   - build simcore in debug mode without graphics"
  echo "  gdebug  - build simcore in debug mode with graphics"
  echo "  test    - build simcore and run unit tests"
  echo "  docs    - build Doxygen documentation"
}

case $1 in
  clean)
    do_clean
    ;;
  build)
    do_build
    ;;
  debug)
    do_debug_build
    ;;
  gbuild)
    do_graph_build
    ;;
  gdebug)
    do_debug_graph_build
    ;;
  test)
    do_test_build
    ;;
  docs)
    do_docs_build
    ;;
  *)
    do_usage
    ;;
esac
