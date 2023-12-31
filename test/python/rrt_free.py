import sys

sys.path.append(".")

import numpy as np
from typing import Tuple, List
import matplotlib.pyplot as plt
import pydynotree as dynotree
import os

np.random.seed(0)

start = np.array([1.5, 1.5])
xlim = [0, 3]
ylim = [0, 3]

num_steps = 100
max_expansion = 1.0


def plot_robot(ax, x, color="black", alpha=1.0):
    ax.plot([x[0]], [x[1]], marker="o", color=color, alpha=alpha)


ub = np.array([xlim[1], ylim[1]])
lb = np.array([xlim[0], ylim[0]])

dynotree.srand(1)
tree = dynotree.TreeR2()
tree.init_tree()
space = tree.getStateSpace()
print(lb, ub)
space.set_bounds(lb, ub)
interpolate_fun = space.interpolate
# interpolate_fun = tree.interpolate
valid_configs = []
tree.addPoint(start, 0, True)
valid_configs.append(start)
parents = []
parents.append(-1)
__xnew = np.zeros(2)
xrand = np.zeros(2)
for i in range(num_steps):
    space.sample_uniform(xrand)
    nn = tree.search(xrand)
    xnear = valid_configs[nn.id]
    advance_rate = min(max_expansion / nn.distance, 1.0)
    interpolate_fun(xnear, xrand, advance_rate, __xnew)
    tree.addPoint(__xnew, len(valid_configs), True)
    valid_configs.append(np.copy(__xnew))
    parents.append(nn.id)

fig, ax = plt.subplots()
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect("equal")

for i, p in enumerate(parents):
    if i != -1:
        ax.plot(
            [valid_configs[i][0], valid_configs[p][0]],
            [valid_configs[i][1], valid_configs[p][1]],
            color="yellow",
            linestyle="dashed",
        )

for i in range(len(valid_configs)):
    plot_robot(ax, valid_configs[i], color="gray")

plot_robot(ax, start, "green")
ax.set_title("build rrt in free space")


if os.environ.get("CI") != "1":
    plt.show()
