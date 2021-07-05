import matplotlib.pyplot as plt
import numpy as np


def plot_trajectory(filename: str):
    fig = plt.figure(figsize=(18, 7), tight_layout=True)

    earth = [plt.Circle((0, 0), 6.378e6, color='lightblue') for i in range(3)]
    axes = []

    for i in range(3):
        axes.append(fig.add_subplot(131 + i))
        axes[i].add_patch(earth[i])
        axes[i].set_ylabel("y position (m)")
        axes[i].set_xlabel("x position (m)")
        axes[i].set_aspect("equal")
        axes[i].grid()

    arr = np.loadtxt(filename, delimiter=",")
    x = arr[:, 1]
    y = arr[:, 2]
    z = arr[:, 3]
    axes[0].plot(x, y, label=filename)
    axes[1].plot(y, z, label=filename)
    axes[2].plot(x, z, label=filename)
    axes[0].set_title("xy")
    axes[1].set_title("yz")
    axes[2].set_title("xz")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    plot_trajectory("./results/RK8.csv")
