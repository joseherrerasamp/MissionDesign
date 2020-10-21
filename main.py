import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from core import *

t_0 = 0.0

T_begin = 2459400.0     # 4th July 2021
T_end = 2459650.0       # 11th March 2022                 # 2459740.0       # 9th June 2022

T_tof_1 = 80            # 80 days of flight
T_tof_1_end = 325       # 325 days of flight

T_tof_2 = 650           # 650 days of flight
T_tof_2_end = 950      # 1100 days of flight


def optimization_task(planet_1, planet_2, r_planet2, K_planet2, planet_3, T_beg, T_fin, T_tof_1_beg, T_tof_1_fin,
                      T_tof_2_beg,
                      T_tof_2_fin):

    x = np.arange(T_beg, T_fin, 5)
    y = np.arange(T_tof_1_beg, T_tof_1_fin, 5)
    z = np.arange(T_tof_2_beg, T_tof_2_fin, 5)
    X, Y, Z = np.meshgrid(x, y, z)

    DV0_departure = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0])))
    DV_arrival_planet_1 = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0])))

    DV0_departure_1 = np.zeros(shape=(len(Y[:, 0, 0]), len(Z[0, 0, :])))
    DV_arrival_planet_2 = np.zeros(shape=(len(Y[:, 0, 0]), len(Z[0, 0, :])))

    DV_1_flyby = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))

    DV_total = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))

    ratio_radius = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))
    fi_flyby = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))
    fi_max = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))

    for i in range(len(X[0, :, 0])):
        for j in range(len(Y[:, 0, 0])):
            print("FECHA DE LANZAMIENTO -----------------------------------")
            print(X[0, i, 0])
            print("TIEMPO DE VUELO 1 ----------------------------------")
            print(Y[j, 0, 0])
            R_1dep, V_1dep, R_2arr, V_2arr, V_hdep, DV_dep, V_harr_1, DV_arr_1 = segment_trajectory(planet_1, planet_2,
                                                                                                    X[0, i, 0], t_0,
                                                                                                    Y[j, 0, 0],
                                                                                                    Y[j, 0, 0])
            DV0_departure[i, j] = np.linalg.norm(DV_dep)
            DV_arrival_planet_1[i, j] = np.linalg.norm(DV_arr_1)

            for k in range(len(Z[0, 0, :])):
                print("TIEMPO DE VUELO 2")
                print(Z[0, 0, k])
                R_2dep, V_2dep, R_3arr, V_3arr, V_hdep_1, DV_dep_1, V_harr_2, DV_arr_2 = segment_trajectory(planet_2,
                                                                                                            planet_3,
                                                                                                            X[0, i, 0],
                                                                                                            t_0 + Y[
                                                                                                                j, 0, 0],
                                                                                                            t_0 + Y[
                                                                                                                j, 0, 0]
                                                                                                            + Z[
                                                                                                                0, 0, k],
                                                                                                            Z[0, 0, k])

                DV0_departure_1[j, k] = np.linalg.norm(DV_dep_1)
                DV_arrival_planet_2[j, k] = np.linalg.norm(DV_arr_2)

                fi, r_p = r_periapse(DV_arr_1, DV_dep_1, K_planet2)
                ratio_radius[i, j, k] = np.absolute(r_p / r_planet2)
                fi_flyby[i, j, k] = fi

                if r_p > r_planet2:
                    DV_flyby = V_hdep_1 - V_harr_1
                    DV_flyby_mag = np.linalg.norm(DV_flyby)
                    DV_1_flyby[i, j, k] = DV_flyby_mag
                    fi_max[i, j, k] = 0
                    # print("radio de periapsis", r_p)
                    # print("Pericentral radius was BIGGER than planet radius")
                elif r_p <= r_planet2:
                    r_angle, V_flyby = fly_by(R_2arr, DV_arr_1, V_2dep, K_planet2, r_planet2)
                    DV_flyby = V_hdep_1 - V_flyby
                    DV_flyby_mag = np.linalg.norm(DV_flyby)
                    DV_1_flyby[i, j, k] = DV_flyby_mag
                    fi_max[i, j, k] = r_angle
                    # print("Maximo angulo de rotacion,", r_angle)
                    # print("radio periapsis", r_p)
                    # print("Pericentral radus was SMALLER than planet radius")
                DV_total[i, j, k] = DV0_departure[i, j] + DV_1_flyby[i, j, k] + DV_arrival_planet_2[j, k]
    """
    plt.figure()
    plt.subplot(1, 1, 1)
    vmin = np.min(DV0_departure)
    v = np.linspace(vmin, 7, 6)
    cs = plt.contour(X[0, :, 0], Y[:, 0, 0], np.transpose(DV0_departure), v, cmap='jet')
    plt.colorbar(cs, label="DV_1 [km/s]")
    plt.title("DV_1 at Departure from Earth")
    plt.xlabel("Departure Date (in Julian Time)")
    plt.ylabel("Time of Flight [days]")
    plt.show()

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.plot_trisurf(X.ravel(), Y.ravel(), Z.ravel(), DV_1_flyby.ravel(), cmap='jet')
    plt.title("DV_flyby")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    fig.colorbar(scat, shrink=0.5, aspect=5, label='DV_flyby [km/s]')
    plt.show()

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.scatter3D(X, Y, Z, c=DV_total.ravel(), cmap='jet', depthshade=0.1)
    plt.title("DV_total")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    fig.colorbar(scat, shrink=0.5, aspect=5, label='DV_total [km/s]')
    plt.show()

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    label = np.where(ratio_radius < 1, 'yellow', 'purple')
    label1 = np.where(ratio_radius < 1, 'purple', 'yellow')
    ax.scatter3D(X, Y, Z, c=label.ravel(), label='r_p < r_planet')
    ax.scatter3D(X, Y, Z, c=label1.ravel(), label='r_p >= r_planet')
    plt.title("Ratio R_flyby to R_planet")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    # fig.colorbar(scat, shrink=0.5, aspect=5, label='R_flyby / R_planet')
    ax.legend(title='R_p/R_planet ratio')
    plt.show()

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.scatter3D(X, Y, Z, c=fi_flyby.ravel(), cmap='viridis', depthshade=0.1)
    plt.title("Fi Angle")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    fig.colorbar(scat, shrink=0.5, aspect=5, label='Fi Angle')

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.scatter3D(X, Y, Z, c=fi_max.ravel(), cmap='coolwarm', depthshade=0)
    plt.title("Maximum Fi Angle")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    fig.colorbar(scat, shrink=0.5, aspect=5, label='Max Fi Angle')
    plt.show()
    """
    min_1 = np.min(DV_total)
    position_min_1 = np.unravel_index(np.argmin(DV_total), DV_total.shape)
    print("Minimum before optimization:", min_1)
    print("Minimum position before optimization:", position_min_1)
    print("Dates of minimum position before optimization")
    print("Departure Date:", X[0, position_min_1[0], 0], "Time of Flight 1:", Y[position_min_1[1], 0, 0],
          "Time of flight 2:", Z[0, 0, position_min_1[2]])

    """

    Optimization of the solution by grid search: Taking the previous minimum, next we reduce the grid search redefining 
    the ranges to look at

    """

    eps = 30

    T_begin_opt = X[0, position_min_1[0], 0] - eps
    T_end_opt = X[0, position_min_1[0], 0] + eps

    T_tof_1_opt = Y[position_min_1[1], 0, 0] - eps
    T_tof_1_end_opt = Y[position_min_1[1], 0, 0] + eps

    T_tof_2_opt = Z[0, 0, position_min_1[2]] - eps
    T_tof_2_end_opt = Z[0, 0, position_min_1[2]] + eps

    x = np.arange(T_begin_opt, T_end_opt, 1)
    y = np.arange(T_tof_1_opt, T_tof_1_end_opt, 1)
    z = np.arange(T_tof_2_opt, T_tof_2_end_opt, 1)
    X, Y, Z = np.meshgrid(x, y, z)

    DV0_departure = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0])))
    DV_arrival_planet_1 = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0])))

    DV0_departure_1 = np.zeros(shape=(len(Y[:, 0, 0]), len(Z[0, 0, :])))
    DV_arrival_planet_2 = np.zeros(shape=(len(Y[:, 0, 0]), len(Z[0, 0, :])))

    DV_1_flyby = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))

    DV_total = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))

    ratio_radius = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))
    fi_flyby = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))
    fi_max = np.zeros(shape=(len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])))

    for i in range(len(X[0, :, 0])):
        for j in range(len(Y[:, 0, 0])):
            # print("FECHA DE LANZAMIENTO -----------------------------------")
            # print(X[0, i, 0])
            # print("TIEMPO DE VUELO 1 ----------------------------------")
            # print(Y[j, 0, 0])
            R_1dep, V_1dep, R_2arr, V_2arr, V_hdep, DV_dep, V_harr_1, DV_arr_1 = segment_trajectory(planet_1, planet_2,
                                                                                                    X[0, i, 0],
                                                                                                    t_0, Y[j, 0, 0],
                                                                                                    Y[j, 0, 0])
            DV0_departure[i, j] = np.linalg.norm(DV_dep)
            DV_arrival_planet_1[i, j] = np.linalg.norm(DV_arr_1)

            for k in range(len(Z[0, 0, :])):
                # print("TIEMPO DE VUELO 2")
                # print(Z[0, 0, k])
                R_2dep, V_2dep, R_3arr, V_3arr, V_hdep_1, DV_dep_1, V_harr_2, DV_arr_2 = segment_trajectory(planet_2,
                                                                                                            planet_3,
                                                                                                            X[0, i, 0],
                                                                                                            t_0 + Y[
                                                                                                                j, 0, 0],
                                                                                                            t_0 + Y[
                                                                                                                j, 0, 0]
                                                                                                            + Z[
                                                                                                                0, 0, k],
                                                                                                            Z[0, 0, k])

                DV0_departure_1[j, k] = np.linalg.norm(DV_dep_1)
                DV_arrival_planet_2[j, k] = np.linalg.norm(DV_arr_2)

                fi, r_p = r_periapse(DV_arr_1, DV_dep_1, K_planet2)
                ratio_radius[i, j, k] = r_p / r_planet2
                fi_flyby[i, j, k] = fi

                if r_p > r_planet2:
                    DV_flyby = V_hdep_1 - V_harr_1
                    DV_flyby_mag = np.linalg.norm(DV_flyby)
                    DV_1_flyby[i, j, k] = DV_flyby_mag
                    fi_max[i, j, k] = 0
                    # print("radio de periapsis", r_p)
                    # print("Pericentral radius was BIGGER than planet radius")
                elif r_p <= r_planet2:
                    r_angle, V_flyby = fly_by(R_2arr, DV_arr_1, V_2dep, K_planet2, r_planet2)
                    DV_flyby = V_hdep_1 - V_flyby
                    DV_flyby_mag = np.linalg.norm(DV_flyby)
                    DV_1_flyby[i, j, k] = DV_flyby_mag
                    fi_max[i, j, k] = r_angle
                    # print("Maximo angulo de rotacion,", r_angle)
                    # print("radio periapsis", r_p)
                    # print("Pericentral radus was SMALLER than planet radius")
                DV_total[i, j, k] = DV0_departure[i, j] + DV_1_flyby[i, j, k] + DV_arrival_planet_2[j, k]
    """
    plt.figure()
    plt.subplot(1, 1, 1)
    vmin = np.min(DV0_departure)
    v = np.linspace(vmin, 10, 11)
    cs = plt.contourf(X[0, :, 0], Y[:, 0, 0], np.transpose(DV0_departure), v, cmap='jet')
    plt.colorbar(cs, label="DV_1 [km/s]")
    plt.title("DV_1 at Departure from Earth")
    plt.xlabel("Departure Date (in Julian Time)")
    plt.ylabel("Time of Flight")
    plt.show()
    

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.scatter3D(X, Y, Z, c=DV_1_flyby.ravel(), cmap='jet', depthshade=1)
    plt.title("DV_flyby Optimized")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    fig.colorbar(scat, shrink=0.5, aspect=5, label='DV_flyby [km/s]')
    plt.show()

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.scatter3D(X, Y, Z, c=DV_total.ravel(), cmap='jet', depthshade=1)
    plt.title("DV_total Optimized")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    fig.colorbar(scat, shrink=0.5, aspect=5, label='DV_total [km/s]')
    plt.show()

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    label = np.where(ratio_radius < 1, 'yellow', 'purple')
    label1 = np.where(ratio_radius < 1, 'purple', 'yellow')
    ax.scatter3D(X, Y, Z, c=label.ravel(), label='r_p < r_planet')
    ax.scatter3D(X, Y, Z, c=label1.ravel(), label='r_p >= r_planet')
    plt.title("Ratio R_flyby to R_planet")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    # fig.colorbar(scat, shrink=0.5, aspect=5, label='R_flyby / R_planet')
    ax.legend(title='R_p/R_planet ratio')

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.scatter3D(X, Y, Z, c=fi_flyby.ravel(), cmap='viridis', depthshade=1)
    plt.title("Fi Angle")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    fig.colorbar(scat, shrink=0.5, aspect=5, label='Fi Angle')

    fig = plt.figure(figsize=(35.0, 35.0))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.scatter3D(X, Y, Z, c=fi_max.ravel(), cmap='coolwarm', depthshade=1)
    plt.title("Maximum Fi Angle")
    ax.set_xlabel("Departure Date (in Julian Time)")
    ax.set_ylabel("Time of Flight [days]")
    ax.set_zlabel("Time of Flight 2 [days]")
    fig.colorbar(scat, shrink=0.5, aspect=5, label='Max Fi Angle')
    plt.show()
    """
    min_opt = np.min(DV_total)
    position_min_opt = np.unravel_index(np.argmin(DV_total), DV_total.shape)
    print("Minimum after optimization", min_opt)
    print("Minimum position after optimization", position_min_opt)
    print("Dates of minimum value after optimization")
    print("Departure Date:", X[0, position_min_opt[0], 0], "Time of Flight 1:", Y[position_min_opt[1], 0, 0],
          "Time of Flight 2:", Z[0, 0, position_min_opt[2]])
    print("-------------------------")
    print("Minimum DV_departure from Earth", DV0_departure[position_min_opt[0], position_min_opt[1]])
    print("Minimum DV_departure from Earth try",  DV0_departure[position_min_opt[1], position_min_opt[0]])
    print("Minimum DV_flyby", DV_1_flyby[position_min_opt[0], position_min_opt[1], position_min_opt[2]])
    print("Minimum DV_flyby try", DV_1_flyby[position_min_opt[1], position_min_opt[0], position_min_opt[2]])
    print("Minimum DV_arrival", DV_arrival_planet_2[position_min_opt[1], position_min_opt[2]])
    print("Minimum DV_arrival try", DV_arrival_planet_2[position_min_opt[0], position_min_opt[2]])



    return X[0, position_min_opt[0], 0], Y[position_min_opt[1], 0, 0], Z[0, 0, position_min_opt[2]]


depart_time, t_tof1, t_tof2 = optimization_task(Earth, Venus, r_Venus, K_Venus, Jupiter, T_begin, T_end, T_tof_1,
                                                T_tof_1_end, T_tof_2, T_tof_2_end)
