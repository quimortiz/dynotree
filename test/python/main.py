

import pykdtree
import numpy as np
import  time 


num_points = 10000

a = pykdtree.TreeRX(2)

a.addPoint(np.array([0, 0]), 0, True)
a.addPoint(np.array([0, 1]), 1, True)
a.addPoint(np.array([0, 2]), 2, True)
nn = a.search(np.array([0, 2.1]))
print(nn)
print(nn.distance)
print(nn.payload)

b = pykdtree.TreeR2(-1)
b.addPoint(np.array([0, 0]), 0, True)
b.addPoint(np.array([0, 1]), 1, True)
b.addPoint(np.array([0, 2]), 2, True)
nn = b.search(np.array([0, 2.1]))
print(nn)
print(nn.distance)
print(nn.payload)

a = pykdtree.TreeRX(2)
b = pykdtree.TreeR2(-1)

for i in range(num_points):
    x = np.random.rand(2)
    a.addPoint(x, i, True)
    b.addPoint(x, i, True)

k = np.array([0.81, 0.15])

tic = time.time()
o = a.searchKnn(k,5)
toc = time.time()
print("elapsed time: ", toc - tic)
# print(o[0].payload)

tic = time.time()
oo = b.searchKnn(k,5)
toc = time.time()
print("elapsed time: ", toc - tic)

# print(o[0].payload)

a = pykdtree.TreeRX(4)
b = pykdtree.TreeR4(-1)

for i in range(num_points):
    x = np.random.rand(4)
    a.addPoint(x, i, True)
    b.addPoint(x, i, True)

k = np.array([0.81, 0.15, .1, .2])

tic = time.time()
o = a.searchKnn(k,2)
toc = time.time()
print("elapsed time: ", toc - tic)
print(o[0].payload)

tic = time.time()
o = b.searchKnn(k,2)
toc = time.time()
print("elapsed time: ", toc - tic)
print(o[0].payload)



a = pykdtree.TreeRX(7)
b = pykdtree.TreeR7(-1)

for i in range(num_points):
    x = np.random.rand(7)
    a.addPoint(x, i, True)
    b.addPoint(x, i, True)

k = np.array([0.81, 0.15, .1, .2, 0.1 , 0.1, .5])

tic = time.time()
o = a.searchKnn(k,2)
toc = time.time()
print("elapsed time: ", toc - tic)
print(o[0].payload)

tic = time.time()
o = b.searchKnn(k,2)
toc = time.time()
print("elapsed time: ", toc - tic)
print(o[0].payload)

# TODO: check against nearest neighbor 




