import sys

sys.path.append(".")

import numpy as np
from typing import Tuple, List
import matplotlib.pyplot as plt
import dynotree

# np.random.seed(0)
dynotree.srand(1)


# ref: https://paulbourke.net/geometry/pointlineplane/
def distance_point_to_segment(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """
    p1, p2: two points defining the segment
    p3: the point
    """
    u = np.dot(p3 - p1, p2 - p1) / np.dot(p2 - p1, p2 - p1)
    u = np.clip(u, 0, 1)
    return np.linalg.norm(p1 + u * (p2 - p1) - p3)


# R2xSO2
# a robot is a segment.


length = 0.5
radius = 0.01


def compute_two_points(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    x: 3D vector (x, y, theta)

    """
    p1 = x[0:2]
    p2 = p1 + length * np.array([np.cos(x[2]), np.sin(x[2])])
    return p1, p2


# env : list( Tuple(np.array, float ) = []

xlim = [0, 3]
ylim = [0, 3]

obstacles = [(np.array([1, 0.4]), 0.5), (np.array([1, 2]), 0.5)]


def plot_env(ax, env):
    for obs in obstacles:
        ax.add_patch(plt.Circle((obs[0]), obs[1], color="blue", alpha=0.5))


def plot_robot(ax, x, color="black", alpha=1.0):
    p1, p2 = compute_two_points(x)
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, alpha=alpha)
    ax.plot([p1[0]], [p1[1]], marker="o", color=color, alpha=alpha)


def is_collision(x: np.ndarray) -> bool:
    """
    x: 3D vector (x, y, theta)

    """
    p1, p2 = compute_two_points(x)
    for obs in obstacles:
        if distance_point_to_segment(p1, p2, obs[0]) < radius + obs[1]:
            return True
    return False


static = False

if static:
    tree = dynotree.TreeR2SO2(-1)
    state_space = tree.getStateSpace()
    state_space.set_bounds(np.array([0, 0]), np.array([3, 3]))
else:
    tree = dynotree.TreeX(3, ["Rn:2", "SO2"])
    state_space = tree.getStateSpace()
    state_space.set_bounds(
        [np.array([0, 0]), np.array([])], [np.array([3, 3]), np.array([])]
    )


fig, ax = plt.subplots()
start = np.array([0.1, 0.1, np.pi / 2])
goal = np.array([2.0, 0.2, 0])


plot_env(ax, obstacles)


plot_robot(ax, start, "green")
plot_robot(ax, goal, "red")


ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect("equal")
ax.set_title("env, start and goal configurations")


plt.show()

# new plot


fig, ax = plt.subplots()

N = 100

x = np.zeros(3)
for i in range(N):
    state_space.sample_uniform(x)
    # x = np.random.rand(3) * np.array([3, 3, 2 * np.pi]) - np.array([0, 0, np.pi])
    is_col = is_collision(x)
    if is_col:
        plot_robot(ax, x, color="gray")
    else:
        plot_robot(ax, x, color="black")


plot_robot(ax, start, "green")
plot_robot(ax, goal, "red")

plot_env(ax, obstacles)

ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect("equal")
ax.set_title("random configurations - collisions")

plt.show()

# now lets do rrt


fig, ax = plt.subplots()
plot_robot(ax, start, "green")
plot_robot(ax, goal, "red")
plot_env(ax, obstacles)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect("equal")
ax.set_title("interpolation")


interpolate_fun = state_space.interpolate
N = 10
for i in range(N):
    out = np.zeros(3)
    interpolate_fun(start, goal, i / N, out)
    plot_robot(ax, out, color="gray")

plt.show()

# now lets solve rrt!
goal_bias = 0.1
max_steps = 1000

check_cols = 10


valid_configs = []
parents = []

tree.addPoint(start, len(valid_configs), True)
valid_configs.append(start)
parents.append(-1)
# num_points_in_tree += 1

max_expansion = 1.0

solved = False
xrand = np.zeros(3)
for i in range(max_steps):
    state_space.sample_uniform(xrand)
    # xrand = np.random.rand(3) * np.array([3, 3, 2 * np.pi]) - np.array([0, 0, np.pi])
    if np.random.rand() < goal_bias:
        xrand = np.copy(goal)
    # find nearest neighbor
    nn = tree.searchKnn(xrand, 1)[0]
    xnear = valid_configs[nn.payload]
    advance_rate = min(max_expansion / nn.distance, 1.0)
    print("advance_rate", advance_rate)
    xnew = np.zeros(3)
    interpolate_fun(xnear, xrand, advance_rate, xnew)
    collision = False
    __tmp = np.zeros(3)
    for j in range(check_cols + 1):
        interpolate_fun(xnear, xnew, j / check_cols, __tmp)
        if is_collision(__tmp):
            collision = True
            break
    if not collision:
        print("adding point")
        tree.addPoint(xnew, len(valid_configs), True)
        valid_configs.append(xnew)
        parents.append(nn.payload)
        if np.linalg.norm(xnew - goal) < 0.01:
            solved = True
            break


fig, ax = plt.subplots()
plot_env(ax, obstacles)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect("equal")


print("done")
print("solved: ", solved)
for i in range(len(valid_configs)):
    print(valid_configs[i])
    plot_robot(ax, valid_configs[i], color="gray")


plot_robot(ax, start, "green")
plot_robot(ax, goal, "red")


for p, i in enumerate(parents):
    if i != -1:
        # plot_robot(ax, valid_configs[i], color="gray")
        # plot_robot(ax, valid_configs[p], color="gray")
        ax.plot(
            [valid_configs[i][0], valid_configs[p][0]],
            [valid_configs[i][1], valid_configs[p][1]],
            color="yellow",
            linestyle="dashed",
        )


# trace back the solution

i = len(valid_configs) - 1
path = []
path.append(np.copy(valid_configs[i]))
while i != -1:
    path.append(np.copy(valid_configs[i]))
    i = parents[i]

# interpolate the path

# reverse
path.reverse()

for i in range(len(path) - 1):
    _start = path[i]
    _goal = path[i + 1]
    for i in range(N):
        out = np.zeros(3)
        interpolate_fun(_start, _goal, i / N, out)
        plot_robot(ax, out, color="gray", alpha=0.5)
# add the last configuration
plot_robot(ax, path[-1], color="gray", alpha=0.5)


for p in path:
    plot_robot(ax, p, color="blue", alpha=1)


print("path", path)

ax.set_title("rrt solution")
plt.show()


# just build rrt
