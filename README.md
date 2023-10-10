# dyn_kdtree

Dynamic Kd-Tree:  Euclidean, SO(2), SO(3) and more!
C++ with Python Bindings

The KD-tree implementation is based on [bucket-pr-kdtree](https://github.com/jkflying/bucket-pr-kdtree)


# Try out

```bash
pip3 install dynotree
```

A first example:
[rrt_free.py](test/python/rrt_free.py)

A second example:
[rrt.py](https://github.com/quimortiz/dyn_kdtree/blob/main/test/python/rrt.py)



# C++

instructions here


# DEV


## Testing



## Create wheels

for publishing in pypy
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
