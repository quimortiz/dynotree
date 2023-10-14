
<p align="center">
  <img width="50%" height="auto" src="https://github.com/quimortiz/dynotree/blob/main/logo.svg">
</p>

# Dynotree

[![pre-commit](https://github.com/quimortiz/dynotree/actions/workflows/pre-commit.yml/badge.svg)](https://github.com/quimortiz/dynotree/actions/workflows/pre-commit.yml)
[![C/C++ CI](https://github.com/quimortiz/dynotree/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/quimortiz/dynotree/actions/workflows/c-cpp.yml)
[![PyPI version](https://badge.fury.io/py/dynotree.svg)](https://badge.fury.io/py/dynotree)

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
pip3 install dynotree
```

First example:
[rrt_free.py](test/python/rrt_free.py)

Second example:
[rrt.py](https://github.com/quimortiz/dyn_kdtree/blob/main/test/python/rrt.py)

### Python from source

Refer to the section on `Creating a Python package and installing from source` below.

## C++

Tests and examples can be found in [main.cpp](https://github.com/quimortiz/dyn_kdtree/blob/main/src/main.cpp).

```bash
git clone --recurse-submodules https://github.com/quimortiz/dyn_kdtree
cd dyn_kdtree
mkdir build
cd build
cmake ..
make
./main
```

To run a specific example or test, use:
```
./main --run_test=NAME_OF_TEST
```

# Development

## Create a Python Package and Install from Source

To create wheels:
```
CMAKE_ARGS="-DBUILD_PYTHON_BINDINGS=1" python3 -m build
```
The wheels will be located in the `dist/` directory. Install the wheels with the following (adjust the package and python version as needed):

```
pip3 install dist/dynotree-0.0.4-cp38-cp38-linux_x86_64.whl
```

## Create a Python Package and Upload to PYPI

### Creating wheels for PYPI

```bash
docker pull quay.io/pypa/manylinux2014_x86_64
docker run -it -v $(pwd):/io quay.io/pypa/manylinux2014_x86_64 /io/build-wheels.sh
python3 -m twine upload wheelhouse/*
```

For uploading to testpypi:
```bash
python3 -m twine upload --repository testpypi wheelhouse/*
```

For creating a local package:
```bash
TODO
```

# Why choose dyn_kdtree?

* It's faster than OMPL and simpler than NIGH.
* It's dynamic: you can add points individually, making it ideal for Motion Planning.
* Supports Euclidean, SO(2), SO(3), and various combinations.
* Offers performant and flexible C++ code based on Eigen.
* Available as a Single Header File.
* Comes with Python Bindings for seamless integration.
* Can be extended to accommodate custom spaces.

# Python Bindings

## Option 1

## Option 2

## Option 3

# Dependencies

For Python: No dependencies.

For C++ Code: Eigen

For Development: Eigen and Boost Testing

# Interface

Refer to the provided examples for more information.

# Benchmark

For benchmarking against other platforms, refer to this repository:
* [bucket-pr-kdtree](https://github.com/jkflying/bucket-pr-kdtree)
* OMPL
* Nigh
* ...

# Roadmap

The code is stable and is currently in use ...
