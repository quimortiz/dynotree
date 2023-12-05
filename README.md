
<p align="center">
  <img width="50%" height="auto" src="https://github.com/quimortiz/dynotree/blob/main/logo.svg">
</p>

# Dynotree

[![pre-commit](https://github.com/quimortiz/dynotree/actions/workflows/pre-commit.yml/badge.svg)](https://github.com/quimortiz/dynotree/actions/workflows/pre-commit.yml)
[![C/C++ CI](https://github.com/quimortiz/dynotree/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/quimortiz/dynotree/actions/workflows/c-cpp.yml)
[![PyPI version](https://badge.fury.io/py/pydynotree.svg)](https://pypi.org/project/pydynotree/)
[![codecov](https://codecov.io/gh/quimortiz/dynotree/graph/badge.svg?token=AFQEOG7QKW)](https://codecov.io/gh/quimortiz/dynotree)

A Dynamic Kd-Tree written in C++ with Python Bindings, supporting Euclidean, SO(2), SO(3), and more!

Dynotree supports both Euclidean and non-Euclidean spaces, as well as compound spaces. These spaces can be defined either at compile time for efficiency or at runtime for flexibility.

It is a header-only C++ library that depends solely on Eigen and offers Python bindings for rapid development and prototyping.

The primary application is motion planning. Check the `examples/python` directory where we implement
RRT, RRT star, and RRT connect in few lines code, utilizing dynotree for efficient nearest neighbor lookup.

This code has been employed in various research projects and is continuously tested using Github Actions. We have benchmarked it against nearest neighbor implementations from Nigh, OMPL, and bucket-pr-kdtree.

We welcome Bug Issues, Pull Requests, and feature suggestions.

The base KD-tree implementation derives from [bucket-pr-kdtree](https://github.com/jkflying/bucket-pr-kdtree), but our version has been enhanced to handle non-Euclidean spaces and support both Compile and Runtime spaces based on the Eigen API.

# Try it out!

## PYTHON

### PIP Package

```bash
pip3 install pydynotree
```

First example:
[rrt_free.py](test/python/rrt_free.py)

Second example:
[rrt.py](https://github.com/quimortiz/dyn_kdtree/blob/main/test/python/rrt.py)

### Python from source

Refer to the section on `Creating a Python package and installing from source` below.

## C++

Tests and examples can be found in [main.cpp](https://github.com/quimortiz/dyn_kdtree/blob/main/src/main.cpp).

TODO
```bash
git clone --recurse-submodules https://github.com/quimortiz/dyn_kdtree
cd dyn_kdtree
mkdir build
cd build
cmake ..
make
./main
```
TODO
To run a specific example or test, use:
```
./main --run_test=NAME_OF_TEST
```

# Development

## Python

### Create a Python Package and Install from Source

Create wheels with
```
CMAKE_ARGS="-DBUILD_PYTHON_BINDINGS=1" python3 -m build
```

The wheels will be located in the `dist/` directory. You can install the wheels with (adjust the package and python version as needed)
```
pip3 install dist/dynotree-0.0.4-cp38-cp38-linux_x86_64.whl
```

### Create a Python Package and Upload to PYPI

Create valid wheels for PYPI (using special docker container)

```bash
docker pull quay.io/pypa/manylinux2014_x86_64
docker run -it -v $(pwd):/io quay.io/pypa/manylinux2014_x86_64 /io/build-wheels.sh
python3 -m twine upload wheelhouse/*
```

For uploading to testpypi:
```bash
python3 -m twine upload --repository testpypi wheelhouse/*
```

### Stubs


Generate stubs from compiled file with pybind
```
pybind12-stubgen pydynotree
```

Copy them in the python package.

From
```
/home/quim/stg/dynotree/src/python/pydynotree
```
run
```
cp ../../../build_release_clang/bindings/python/stubs/pydynotree.pyi
```



## C++

Run tests with
```
TODO
```

# Documentation

[C++ Documentation](https://quimortiz.github.io/dynotree/index.html)

[Python Documentation](https://quimortiz.github.io/dynotree_pydoc/index.html)


The documentation is quite minimal. The goal is to provide a searchable index that exposes all available classes and functions.
You can instead check examples and tests to see how to use the code.
C++ documentation is built with Doxygen and GitHub CI.
Python documentation is built and deployed in another repository [dynotree_pydoc](https://github.com/quimortiz/dynotree_pydoc) (To Do -- should be built automatically from the PIP package).

To build the documentation locally:

```
doxygen
firefox doxydoc/html/index.html
```




# Dependencies

PIP Package: No dependencies (but example require `numpy` and `matplotlib`)

C++ header-only Libray - Eigen3.

For C++ test and benchmark: Eigen and Boost Testing


# Benchmark

For benchmarking against other platforms, refer to this repository:
* [bucket-pr-kdtree](https://github.com/jkflying/bucket-pr-kdtree)
* OMPL
* Nigh
* ...

Run benchmark with:

```
TODO
```


# Code Coverage

(work in progress)
```
cmake  -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=1 -DENABLE_TEST_COVERAGE=1 ..
make -j8
ctest -C Release -VV
gcov ./test/cpp/CMakeFiles/main.dir/main.cpp.gcno
lcov --capture --directory . --output-file coverage.info
lcov --remove coverage.info '/usr/include/*' 'src/*' '/usr/lib/*' -o filtered_coverage.info
genhtml filtered_coverage.info --output-directory out
find . -type f -name '*.html'
firefox out/index.html
```

(TODO: get this into github CI...)



# Roadmap

The code is stable and currently in used in Dynoplan. No API breaking changes are expected but ther are not guarantees.
