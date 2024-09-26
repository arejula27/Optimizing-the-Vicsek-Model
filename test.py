import numpy as np

"""
Create Your Own Active Matter Simulation (With Python)
Philip Mocz (2021) Princeton Univeristy, @PMocz

Simulate Viscek model for flocking birds

"""


def main():
    """Finite Volume simulation"""

    # Simulation parameters
    R = 1  # interaction radius
    Nt = 1  # number of time steps
    # number of birds

    x = np.array([0.0, 10.0, 100.0, 200.0])
    #  y = {1, 1.0, 0.0, 0.3};
    y = np.array([0.0, 10.0, 100.0, 2000.0])

    N = len(x)  # number of birds
    # bird velocities
    # theta = {1.0, 1.0, 1.0, 1};
    theta = np.ones(N)

    # find mean angle of neighbors within R
    mean_theta = theta
    for b in range(N):
        #        print((x-x[b])**2+(y-y[b])**2)
        neighbors = (x - x[b]) ** 2 + (y - y[b]) ** 2 < R**2
        sx = np.sum(np.cos(theta[neighbors]))
        # print(np.cos(theta[neighbors]))
        sy = np.sum(np.sin(theta[neighbors]))
        # print("=====")
        # print(sx)
        # print(sy)
        # print(np.arctan2(sy, sx))
        mean_theta[b] = np.arctan2(sy, sx)
    #    print(mean_theta)

    # add random perturbations
    theta = mean_theta
    print(theta)

    return 0


if __name__ == "__main__":
    main()
