import numpy as np
from planet_consts import *

# Consts
# mu = 398600
mu = K_Sun

tolerance = 10 ** (-8)
au = 1.496 * 10 ** 8                            # km in 1 Astronomical Units
h_atm = 600                                     # [km]

# Decision parameter for Lambert's Problem
kind_orbit = 1                                  # 1: Prograde Orbit / 2: Retrograde Orbit

"""
Some useful generic functions
"""


def romin(r_planet):
    return r_planet + h_atm


"""
State Vector to Orbital Elements
    h       [km2 / s]
    i       [rad]
    RA      [rad]
    w       [rad]
    theta   [rad]
    e       []
    rp      [km]
    ra      [km]
"""


def vstate_to_coe(vect1, vect2):
    r = np.linalg.norm(vect1)
    # v = np.linalg.norm(vect2)
    v_rad = np.dot(vect1, vect2) / r

    H = np.cross(vect1, vect2)
    h = np.linalg.norm(H)

    inc = np.arccos(H[2] / h)

    N = np.cross([0, 0, 1], H)
    n = np.linalg.norm(N)

    if n != 0:
        RA = np.arccos(N[0] / n)
        if N[1] < 0:
            RA = 2 * np.pi - RA
    else:
        RA = 0

    e_vect = (1 / mu) * (np.cross(vect2, H) - mu * (vect1 / r))
    e = np.linalg.norm(e_vect)

    if n != 0:
        w = np.arccos(np.dot(N, e_vect) / (n * e))
        if e_vect[2] < 0:
            w = 2 * np.pi - w
    else:
        w = 0

    theta = np.arccos(np.dot(e_vect, vect1) / (e * r))
    if v_rad < 0:
        theta = 2 * np.pi - theta

    T = (2.0 * np.pi / mu ** 2.0) * (h / np.sqrt(1 - e ** 2.0)) ** 3.0

    return h, e, RA, inc, w, theta, T


"""

Orbital Elements to State Vector

"""


def coe_to_vstate(h, e, RA, inc, w, theta):
    R_pf = (h ** 2.0 / mu) * (1 / (1 + e * np.cos(theta))) * np.array([np.cos(theta), np.sin(theta), 0.0])
    V_pf = (mu / h) * np.array([-np.sin(theta), e + np.cos(theta), 0.0])

    R3_RA = np.array([[np.cos(RA), np.sin(RA), 0],
                      [-np.sin(RA), np.cos(RA), 0],
                      [0, 0, 1]])

    R1_inc = np.array([[1, 0, 0],
                       [0, np.cos(inc), np.sin(inc)],
                       [0, -np.sin(inc), np.cos(inc)]])

    R3_w = np.array([[np.cos(w), np.sin(w), 0],
                     [-np.sin(w), np.cos(w), 0],
                     [0, 0, 1]])

    conv_matrix = R3_w @ R1_inc @ R3_RA
    conv_matrix = np.transpose(conv_matrix)

    R = conv_matrix @ R_pf
    V = conv_matrix @ V_pf

    R = np.transpose(R)
    V = np.transpose(V)

    return R, V


"""

Kepler Equation


"""


def kepler_eliptic(e, M):
    if M < np.pi:
        E = M + e / 2.0
    else:
        E = M - e / 2.0

    ratio = 1
    while abs(ratio) > tolerance:
        fE = E - e * np.sin(E) - M
        df_E = 1 - e * np.cos(E)
        ratio = fE / df_E
        E = E - ratio

    theta = 2.0 * np.arctan(np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(E / 2.0))
    if theta < 0:
        theta = 2 * np.pi + theta

    return E, theta


"""

SOLAR SYSTEM 

"""


def julian_date(t_ephemeris, t):                            # t in days
    t_julian = (t_ephemeris - 2451545.0) / 36525.0
    t_julian = t_julian + (t / 36525.0)

    return t_julian


"""
planets_elements function returns Keplerian Elements for each planet
    h_p     [km2/s]
    a_p     [km]
    e_p     [rad]
    I_p     [rad]
    L_p     [rad]
    w1_p    [rad]
    w_p     [rad]
    o_p     [rad]
    M_p     [rad]
"""


def planets_elements(planet, julian_time):
    a_p = planet[0] + (planet[1] * julian_time)         # [au]
    e_p = planet[2] + (planet[3] * julian_time)         # [rad]
    I_p = planet[4] + (planet[5] * julian_time)         # [º]
    L_p = planet[6] + (planet[7] * julian_time)         # [º]
    w1_p = planet[8] + (planet[9] * julian_time)        # [º]
    O_p = planet[10] + (planet[11] * julian_time)       # [º]

    # Every parameter should be given in radians
    # Also a_p should be transformed to [km] in order to have the same units

    a_p = a_p * au                                      # [km]
    h_p = np.sqrt(mu * a_p * (1 - e_p ** 2))            # [km2/s]
    I_p = np.radians(I_p)                               # [rad]
    L_p = np.radians(L_p)                               # [rad]
    w1_p = np.radians(w1_p)                             # [rad]
    O_p = np.radians(O_p)                               # [rad]

    w_p = w1_p - O_p                                    # [rad]

    M_p = L_p - w1_p + np.radians(planet[12] * (julian_time ** 2.0)) + np.radians(planet[13]) * (
        np.cos(np.radians(planet[15] * julian_time))) + np.radians(planet[14]) * (
              np.sin(np.radians(planet[15] * julian_time)))

    return h_p, a_p, e_p, I_p, L_p, w1_p, w_p, O_p, M_p


def planet_state(planet, T_ephemeris, t):                   # t in days
    t_julian = julian_date(T_ephemeris, t)

    h_p, a_p, e_p, I_p, L_p, w1_p, w_p, O_p, M_p = planets_elements(planet, t_julian)

    E_p, theta_p = kepler_eliptic(e_p, M_p)

    R_planet, V_planet = coe_to_vstate(h_p, e_p, O_p, I_p, w_p, theta_p)

    return R_planet, V_planet


"""

LAMBERT'S PROBLEM: Calculation of the orbit between two given points

        R1 [km] Initial position of spacecraft in its trajectory. Coincides with planet position at time t1
        R2 [km] Desired position of spacecraft at the end. Coincides with planet position at time t2

"""


# Stumpff Functions

def stumpfS(z):
    if z > 0:
        S = (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z)) ** 3
    elif z < 0:
        S = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z)) ** 3
    else:
        S = 1 / 6

    return S


def stumpfC(z):
    if z > 0:
        C = (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        C = (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        C = 1 / 2

    return C


def lambert_problem(R1, R2, fly_time):                                  # fly_time is given in days
    R1_norm = np.linalg.norm(R1)                                        # [km]
    R2_norm = np.linalg.norm(R2)                                        # [km]

    cross_product = np.cross(R1, R2)
    theta_lambert = np.arccos(np.dot(R1, R2) / (R1_norm * R2_norm))     # [rad]

    if kind_orbit == 1:                                                 # 1: prograde orbit
        if cross_product[2] < 0:
            theta_lambert = 2 * np.pi - theta_lambert                   # [rad]
    elif kind_orbit == 2:                                               # 2: retrograde orbit
        if cross_product[2] >= 0:
            theta_lambert = 2 * np.pi - theta_lambert                   # [rad]

    A = np.sin(theta_lambert) * np.sqrt((R1_norm * R2_norm) / (1 - np.cos(theta_lambert)))

    # Additional Functions needed to calculate the variables involved in Lambert Problems
    def y_equation(z):
        S = stumpfS(z)
        C = stumpfC(z)

        y_z = R1_norm + R2_norm + (A * (z * S - 1)) / np.sqrt(C)

        return y_z

    def Fz(z):
        y = y_equation(z)
        S = stumpfS(z)
        C = stumpfC(z)

        if y < 0:
            F = -1
        else:
            F = (((y / C) ** 1.5) * S) + (A * np.sqrt(y)) - (np.sqrt(mu) * (fly_time * 24 * 3600))

        return F

    def dFz(z):
        y = y_equation(z)
        S = stumpfS(z)
        C = stumpfC(z)

        if z == 0:
            dF = (np.sqrt(2) / 40 * (y ** 1.5)) + A / 8 * (np.sqrt(y) + A * np.sqrt(1 / (2 * y)))
        elif z != 0:
            dF = ((y / C) ** 1.5) * ((1 / (2 * z)) * (C - ((3 * S) / (2 * C))) + ((3 * S ** 2) / (4 * C))) + (A / 8) * (
                    (3 * S / C * np.sqrt(y)) + (A * np.sqrt(C / y)))
        return dF

    z = -100
    F = Fz(z)

    while F < 0:
        z = z + 0.1
        F = Fz(z)

    tol = 10 ** (-8)
    nmax = 5000

    ratio = 1
    n = 0

    while (np.absolute(ratio) > tol) and (n <= nmax):
        n = n + 1
        F = Fz(z)
        dF = dFz(z)
        ratio = F / dF
        z = z - ratio

    if n >= nmax:
        print("Number of iterations exceeds", nmax)

    f = 1 - y_equation(z) / R1_norm
    g = A * np.sqrt(y_equation(z) / mu)
    dg = 1 - y_equation(z) / R2_norm

    V1 = 1 / g * (R2 - f * R1)
    V2 = 1 / g * (dg * R2 - R1)

    return V1, V2


'''

Function which returns one segment of the mission

    R_planet1, V_planet1       Defines the planet state vector at the beginning of the transference 
    R_planet2, V_planet2       Defines the planet state vector at the end of the transference
    V_departure                Defines the spacecraft heliocentric velocity when it leaves the planet at the beginning 
                                of the transference
    V_arrival                  Defines the spacecraft heliocentric velocity when it arrives the planet at the end of the
                                transference
    DV_departure               Defines the spacecraft planetocentric velocity when it leaves the planet at the beginning 
                                of the transference
    DV_arrival                 Defines the spacecraft planetocentric velocity when it arrives the planet at the end of
                                the transference
'''


def segment_trajectory(planet_1, planet_2, T_ephem_begin, t_0, t_1, fly_time):
    R_planet1, V_planet1 = planet_state(planet_1, T_ephem_begin, t_0)               # [km]
    R_planet2, V_planet2 = planet_state(planet_2, T_ephem_begin, t_1)               # [km]

    V_departure, V_arrival = lambert_problem(R_planet1, R_planet2, fly_time)        # [km/s]

    DV_departure = V_departure - V_planet1                                          # [km/s]
    DV_arrival = V_arrival - V_planet2                                              # [km/s]

    return R_planet1, V_planet1, R_planet2, V_planet2, V_departure, DV_departure, V_arrival, DV_arrival


"""

GRAVITY ASSIST MANEUVER: rotation matrix for velocity and position vectors are calculated

"""


# Function which returns the periapse radius


def r_periapse(v_arrival, v_departure, K_planet):
    # v_arrival defines the planetocentric velocity of the spacecraft at the inbound crossing point of the flyby planet
    # v_departure defines the planetocentric velocity of the spacecraft at the outbound crossing point of this planet
    rotated_angle = np.dot(v_arrival, v_departure) / (np.linalg.norm(v_arrival) * np.linalg.norm(v_departure))
    rotated_angle = np.arccos(rotated_angle)

    r_p = (K_planet * ((1 / np.sin(rotated_angle / 2)) - 1)) / (np.linalg.norm(v_arrival)) ** 2

    return rotated_angle, r_p


def rot_vel_angle(v_infinity, K_planet, r_planet):  # r, v are given in planetocentric coordinates

    r_min = romin(r_planet)

    rotation_v_angle = 2 * np.arcsin(K_planet / (K_planet + r_min * (v_infinity ** 2)))  # [rad]

    return rotation_v_angle


def rotation_vel_matrix(r, v, K_planet, r_planet):
    c1 = r[1] * v[2] - r[2] * v[1]                              # [km2/s]
    c2 = r[2] * v[0] - r[0] * v[2]                              # [km2/s]
    c3 = r[0] * v[1] - r[1] * v[0]                              # [km2/s]
    c = np.sqrt((c1 ** 2.0) + (c2 ** 2.0) + (c3 ** 2.0))        # [km2/s]

    v_infinity = np.linalg.norm(v)

    fi = rot_vel_angle(v_infinity, K_planet, r_planet)          # [rad]

    return np.array([[np.cos(fi), (c3 / c) * np.sin(fi), -(c2 / c) * np.sin(fi)],
                     [-(c3 / c) * np.sin(fi), np.cos(fi), (c1 / c) * np.sin(fi)],
                     [(c2 / c) * np.sin(fi), -(c1 / c) * np.sin(fi), np.cos(fi)]])  # no units


"""

fly_by function will describe the fly by 

        r_1 defines spacecraft position at the beginning of the fly by in planetocentric coors
        v_1 defines spacecraft velocity at the beginning of the fly by in planetocentric coors

        r_2 defines spacecraft position at the end of the fly by in planetocentric coors
        v_2 defines spacecraft position at the end of the fly by in planetocentric coors

        R_2 defines spacecraft position at the end of the fly by in heliocentric coors
        V_2 defines spacecraft velocity at the end of the fly by in heliocentric coors

"""


def fly_by(r_1, v_1, V_planet, K_planet, r_planet):
    v_infinity = np.linalg.norm(v_1)
    v_2 = rotation_vel_matrix(r_1, v_1, K_planet, r_planet) @ v_1.T
    matrix = rotation_vel_matrix(r_1, v_1, K_planet, r_planet)
    # R_2 = R_planet + r_2
    V_2 = V_planet + v_2

    rot_angle = rot_vel_angle(v_infinity, K_planet, r_planet)

    return rot_angle, V_2


"""
R_Earth, V_Earth = planet_state(Earth, 2450394.500)
R_Earth1, V_Earth1 = planet_state(Earth, 2452879.000)
R_Mars, V_Mars = planet_state(Mars, 2450703.500)

print("Earth, first case")
print(R_Earth, V_Earth)
print("Earth, second case")
print(R_Earth1, V_Earth1)
print("Mars, case")
print(R_Mars, V_Mars)
"""
