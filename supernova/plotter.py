from __future__ import annotations

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np


def set_axes_equal(ax: Axes3D):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array(
        [
            ax.get_xlim3d(),
            ax.get_ylim3d(),
            ax.get_zlim3d(),
        ]
    )
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)


def _set_axes_radius(ax: Axes3D, origin: np.ndarray, radius: np.ndarray):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])


def plot_from_array(t: np.ndarray, arr: np.ndarray) -> None:
    fig = plt.figure(figsize=(12, 5), tight_layout=True)

    earth = [plt.Circle((0, 0), 6.378e6, color="lightblue") for i in range(3)]
    axes: list[plt.Axes] = []

    for i in range(3):
        axes.append(fig.add_subplot(131 + i))
        axes[i].add_patch(earth[i])
        axes[i].set_aspect("equal")
        axes[i].grid()

    x = arr[:, 0]
    y = arr[:, 1]
    z = arr[:, 2]
    axes[0].plot(x, y, color="green")
    axes[0].set_ylabel("y position (m)")
    axes[0].set_xlabel("x position (m)")
    axes[1].plot(y, z, color="green")
    axes[1].set_ylabel("z position (m)")
    axes[1].set_xlabel("y position (m)")
    axes[2].plot(x, z, color="green")
    axes[2].set_ylabel("z position (m)")
    axes[2].set_xlabel("x position (m)")
    axes[0].set_title("xy")
    axes[1].set_title("yz")
    axes[2].set_title("xz")

    plt.suptitle("Spacecraft Position, Inertial Coordinates")

    plt.show()


def plot_3d_from_array(
    t: np.ndarray, arr: np.ndarray, array_desampling_factor: int = 1, earth_n: int = 20
):
    fig = plt.figure(figsize=(6, 6), tight_layout=True)
    ax = plt.axes(projection="3d")

    x = arr[::array_desampling_factor, 0]
    y = arr[::array_desampling_factor, 1]
    z = arr[::array_desampling_factor, 2]

    ax.plot(x, y, z, color="green")

    # Set the labels
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")

    u = np.linspace(0, 2 * np.pi, earth_n)
    v = np.linspace(0, np.pi, earth_n)
    x = np.outer(np.cos(u), np.sin(v))  # np.outer() -> outer vector product
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x * 6378e3, y * 6378e3, z * 6378e3, alpha=0.5, color="lightblue")

    # Set the title
    ax.set_title("Spacecraft Position, Inertial Coordinates")

    ax.set_box_aspect([1, 1, 1])  # IMPORTANT - this is the new, key line
    # ax.set_proj_type('ortho') # OPTIONAL - default is perspective (shown in image above)
    set_axes_equal(ax)  # IMPORTANT - this is also required

    plt.show()


def plot_error(diff: np.ndarray):
    fig = plt.figure(figsize=(18, 7), tight_layout=True)

    axes = []

    for i in range(3):
        axes.append(fig.add_subplot(131 + i))
        axes[i].set_ylabel("Solver Error (m)")
        axes[i].set_xlabel("Integrator Step Number")
        axes[i].grid()

    for i in range(3):
        axes[i].plot(diff[:, i])

    axes[0].set_title("x error")
    axes[1].set_title("y error")
    axes[2].set_title("z error")

    plt.show()
