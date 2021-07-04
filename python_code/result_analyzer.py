import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    plt.figure(figsize=(7, 7))

    plt.title("Trajectory")
    plt.ylabel("y position (m)")
    plt.xlabel("x position (m)")

    earth = plt.Circle((0, 0), 6.378e6, color='lightblue')
    fig = plt.gcf()
    ax = fig.gca()
    ax.add_patch(earth)

    plt.axis("equal")
    plt.grid()

    # EULER
    arr = np.loadtxt("c_code/euler.csv", delimiter=",")
    x = arr[:, 1]
    y = arr[:, 2]
    plt.plot(x, y, label="Euler")

    # RK4
    arr = np.loadtxt("c_code/RK4.csv", delimiter=",")
    x = arr[:, 1]
    y = arr[:, 2]
    plt.plot(x, y, label="RK4")

    # RK8
    arr = np.loadtxt("c_code/RK8.csv", delimiter=",")
    x = arr[:, 1]
    y = arr[:, 2]
    plt.plot(x, y, label="RK8")

    plt.legend()


    plt.show()
