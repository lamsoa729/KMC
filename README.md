# KMCX - Kinetic Monte Carlo for Crosslinkers

## About

KMCX is code for calculating proper binding and unbinding behavior for crosslinker motor proteins on rod-like cytoskeletal filaments.

## Compiling

Requirements are cmake, a C++11 compiler, and the GSL library installed somewhere where cmake can find it (like */usr/local/lib*, etc).

To compile using cmake, do the usual 

```
mkdir build
cd build
cmake ..
make
```

## Documentation

By default, documentation will be built into the *build/documentation* folder in both html and latex using Doxygen (if it is available).

## Testing

Unit tests are written with catch2 and located in the *tests* folder. They are compiled by default. After building with cmake, run

```
make test
```

If any of the tests fail, the offending unit test binary should be run with `./tests/test_<rng_name>` to see which exact unit test failed.

