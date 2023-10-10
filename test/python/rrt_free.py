import sys
sys.path.append('.')

import numpy as np
from typing import Tuple, List
import matplotlib.pyplot as plt
import dynotree

np.random.seed(0)

start = np.array([1.5, 1.5])
xlim = [0, 3]
ylim = [0, 3]

num_steps = 100
max_expansion = 1.


def plot_robot(ax, x, color="black", alpha=1.0):
    ax.plot([x[0]], [x[1]], marker="o", color=color, alpha=alpha)

def get_x_rand(): return np.random.rand(2) * (ub - lb) + lb

ub = np.array([xlim[1], ylim[1]])
lb = np.array([xlim[0], ylim[0]])

tree = dynotree.TreeR2(-1)
interpolate_fun = tree.interpolate
valid_configs = []
tree.addPoint(start, 0, True)
valid_configs.append(start)
parents = []
parents.append(-1)
__xnew = np.zeros(2)
for i in range(num_steps):
    xrand = get_x_rand()
    nn = tree.searchKnn(xrand, 1)[0]
    xnear = valid_configs[nn.payload]
    advance_rate = min(max_expansion / nn.distance, 1.)
    interpolate_fun(xnear, xrand, advance_rate, __xnew)
    tree.addPoint(__xnew, len(valid_configs), True)
    valid_configs.append(np.copy(__xnew))
    parents.append(nn.payload)

fig, ax = plt.subplots()
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect('equal')

for i, p in enumerate(parents):
    if i != -1:
        ax.plot([valid_configs[i][0], valid_configs[p][0]], [
                valid_configs[i][1], valid_configs[p][1]], color="yellow", linestyle="dashed")

for i in range(len(valid_configs)):
    plot_robot(ax, valid_configs[i], color="gray")

plot_robot(ax, start, "green")
ax.set_title("build rrt in free space")
plt.show()
