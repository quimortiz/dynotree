# dyn_kdtree

Dynamic Kd-Tree in C++ and Python: Euclidean, SO(2), SO(3) and more!

The underlying KD-tree implementation is based on [bucket-pr-kdtree](https://github.com/jkflying/bucket-pr-kdtree), but dynotree supports non-Euclidean spaces and custom compound spaces. State spaces can be defined at both compile and runtime for both efficiency and flexibility.

The C++ library is header-only, with Eigen as single dependency, and we provide Python bindings for fast development and prototyping.


# Try it out!

## PYTHON

```bash
pip3 install dynotree
```

A first example:
[rrt_free.py](test/python/rrt_free.py)

A second example:
[rrt.py](https://github.com/quimortiz/dyn_kdtree/blob/main/test/python/rrt.py)

## C++

Test and examples are in  [main.cpp](https://github.com/quimortiz/dyn_kdtree/blob/main/src/main.cpp)

```bash
git clone --recurse-submodules https://github.com/quimortiz/dyn_kdtree
cd dyn_kdtree
mkdir build
cd build
cmake ..
make
./main
```

use
```
./main --run_test=NAME_OF_TEST
```
to run only one example or test


# DEV


## Testing



## Create wheels

for publishing in pypi
```bash
docker pull quay.io/pypa/manylinux2014_x86_64
docker run -it -v $(pwd):/io quay.io/pypa/manylinux2014_x86_64 /io/build-wheels.sh
python3 -m twine upload wheelhouse/*
```

using testpypi
```bash
python3 -m twine upload --repository testpypi wheelhouse/*
```


for creating local package:
```bash
TODO
```


# Why dyn_kdtree?

* Faster than OMPL and simpler than NIGH
* Dynamic: Add points one by one -- Ideal for Motion Planning
* Support Euclidean, SO(2), SO(3) and any combination
* Performant and flexible C++ code based on Eigen.
* Single Header File
* Python Bindings for easy integration
* Extendable with custom spaces!


# Python Bindings

## Option 1

## Option 2

## Option 3


# Dependencies

Python: No dependencies


C++ Code: Eigen

Develop: Eigen and Boost Testing

# Interface

See this for examples.

# Benchmark

Check this repo for benchmark against
* [bucket-pr-kdtree](https://github.com/jkflying/bucket-pr-kdtree)
* OMPL
* Nigh
* ...


# Roadmap

Code is Stable, currently used in ...
