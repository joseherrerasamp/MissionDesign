import numpy as np
from core import vstate_to_coe, coe_to_vstate, kepler_eliptic, mu
from core import julian_date, planets_elements, lambert_problem, planet_state, segment_trajectory, r_periapse, fly_by
from planet_consts import *


"""
Mutual Validation of the function vstate_to_coe/coe_to_vstate
"""

R1 = np.array([-6045.0, -3490.0, 2500.0])
V1 = np.array([-3.457, 6.618, 2.533])

R2 = np.array([-4039.9, 4814.56, 3628.62])
V2 = np.array([-10.386, -4.77192, 1.74388])

R3 = np.array([5000, 10000, 2100])
V3 = np.array([-5.99249, 1.92536, 3.24564])

R4 = np.array([-14600, 2500, 7000])
V4 = np.array([-3.31246, -4.19662, -0.385288])

R5 = np.array([3830.68, -2216.47, 6605.09])
V5 = np.array([1.50357, -4.56099, -0.291536])


def valid_vstate_to_coe_coe_tovstate():
    val1 = vstate_to_coe(R1, V1)
    rval1 = coe_to_vstate(val1[0], val1[1], val1[2], val1[3], val1[4], val1[5])
    val2 = vstate_to_coe(R2, V2)
    rval2 = coe_to_vstate(val2[0], val2[1], val2[2], val2[3], val2[4], val2[5])
    val3 = vstate_to_coe(R3, V3)
    rval3 = coe_to_vstate(val3[0], val3[1], val3[2], val3[3], val3[4], val3[5])
    val4 = vstate_to_coe(R4, V4)
    rval4 = coe_to_vstate(val4[0], val4[1], val4[2], val4[3], val4[4], val4[5])
    val5 = vstate_to_coe(R5, V5)
    rval5 = coe_to_vstate(val5[0], val5[1], val5[2], val5[3], val5[4], val5[5])

    validation = [val1, val2, val3, val4, val5]
    rvalidation = [rval1, rval2, rval3, rval4, rval5]

    for i in range(len(validation)):
        print("Cartesian to Orbital Elements Case", i + 1, ":")
        print(validation[i])
        print("")
        print("Orbital Elements to Cartesian Case", i + 1, ":")
        print(rvalidation[i])
        print("")


"""
Validation of Kepler Function through kepler_equation function
"""


e1 = 0.37255
M1 = 3.6029

e2 = 0.0167088
M2 = 4.054539479


e3 = 0.98
M3 = 2.53072

e4 = 0.1
M4 = 0.3458

e5 = 0.5
M5 = 0.471238898


def validation_keplereq():

    val1 = kepler_eliptic(e1, M1)
    val2 = kepler_eliptic(e2, M2)
    val3 = kepler_eliptic(e3, M3)
    val4 = kepler_eliptic(e4, M4)
    val5 = kepler_eliptic(e5, M5)

    validation = [val1, val2, val3, val4, val5]

    for i in range(len(validation)):
        print("Validation of Kepler Equation", i + 1, ":")
        print("Eccentric Anomaly, True Anomaly =", validation[i])
        print("")


'''
Validation of Julian Date code which returns the time pass since the beginning in the Julian Calendar
t is given in seconds
'''

T_eph = 2458485.0

tJ1 = 100

tJ2 = 1504.5

tJ3 = 43200

tJ4 = 86400

T_eph_2 = 2455050.0
tJ5 = 900576.25


def validation_juliantime():
    val1 = julian_date(T_eph, tJ1)
    val2 = julian_date(T_eph, tJ2)
    val3 = julian_date(T_eph, tJ3)
    val4 = julian_date(T_eph, tJ4)
    val5 = julian_date(T_eph_2, tJ5)

    validation = [val1, val2, val3, val4, val5]

    for i in range(len(validation)):
        print("Julian time validation Case", i, ":")
        print(validation[i])


T_eph1 = 2452879.0
T_eph2 = 2504068.5
T_eph3 = 2453719.5
T_eph4 = 2457158.0
T_eph5 = 2451545.0

tjulian1 = julian_date(T_eph1, 0)
tjulian2 = julian_date(T_eph2, 0)
tjulian3 = julian_date(T_eph3, 0)
tjulian4 = julian_date(T_eph5, 0)
tjulian5 = julian_date(T_eph5, 0)


def validation_kepelements_equation():
    val1 = planets_elements(Earth, tjulian1)
    equ1 = kepler_eliptic(val1[2], val1[8])
    val2 = planets_elements(Mars, tjulian2)
    equ2 = kepler_eliptic(val2[2], val2[8])
    val3 = planets_elements(Mars, tjulian3)
    equ3 = kepler_eliptic(val3[2], val3[8])
    val4 = planets_elements(Jupiter, tjulian4)
    equ4 = kepler_eliptic(val4[2], val4[8])
    val5 = planets_elements(Pluto, tjulian5)
    equ5 = kepler_eliptic(val5[2], val5[8])

    validation = [val1, val2, val3, val4, val5]
    validequation = [equ1, equ2, equ3, equ4, equ5]

    for i in range(len(validation)):
        print("Validation", i + 1, ":")
        print(validation[i])
        print("Validation of the Kepler Equation", i + 1, ":")
        print(validequation[i])


planet1 = Earth
T_try_1 = 2452879.0

planet2 = Earth
T_try_2 = 2450394.5

planet3 = Mars
T_try_3 = 2450703.5


def validation_planet_state():
    val1 = planet_state(planet1, T_try_1, 0)
    val2 = planet_state(planet2, T_try_2, 0)
    val3 = planet_state(planet3, T_try_2, 309)          # 309 days

    validation = [val1, val2, val3]

    for i in range(len(validation)):
        print("Validation", i + 1, ":")
        print(validation[i])


R_1 = np.array([5000, 10000, 2100])
R_2 = np.array([-14600, 2500, 7000])
fly_time1 = 3600 / 86400                                                # days

R2_1 = np.array([1.04994*10**8, 1.04655*10**8, 988.331])
R2_2 = np.array([-2.08329*10**7, -2.18404*10**8, -4.06287*10**6])
fly_time2 = 309                                                         # days

R3_1 = np.array([7231.58074563487, 218.02523761425, 11.79251215952])
R3_2 = np.array([7357.06485698842, 253.55724281562, 38.81222241557])
fly_time3 = 12300 / 86400                                               # days

R4_1 = np.array([22592.145603, -1599.915239, -19783.950506])
R4_2 = np.array([1922.067697, 4054.157051, -8925.727465])
fly_time4 = 36000 / 86400                                               # days


def validation_lambert():
    v1 = lambert_problem(R_1, R_2, fly_time1)
    v2 = lambert_problem(R2_1, R2_2, fly_time2)
    v3 = lambert_problem(R3_1, R3_2, fly_time3)
    v4 = lambert_problem(R4_1, R4_2, fly_time4)

    validation = [v1, v2, v3, v4]

    for i in range(len(validation)):
        print("validation", i + 1, ":")
        print(validation[i])


# Validation of function which returns a segment trajectory. This functions uses other defined and validated functions
# This is just a try to verify it works


def validation_segment_trajectory():
    t_0 = 0
    R_1p, V_1p, R_2p, V_2p, V_d, DV_d, V_a, DV_a = segment_trajectory(Earth, Mars, 2450394.500, t_0, 309, 309)
    print("Earth position at beginning,", R_1p)
    print("Earth velocity at beginning,", V_1p)
    print("Mars position at the end,", R_2p)
    print("Mars velocity at the end,", V_2p)
    print("Heliocentric Velocity at departure;", V_d)
    print("Planetocentric Velocity at departure;", DV_d)
    print("Heliocentric Velocity at arrival;", V_a)
    print("Planetocentric Velocity at arrival;", DV_a)


# print("Validation of function which transforms Cartesian Elements into Orbital Elements and viceversa")
# print(valid_vstate_to_coe_coe_tovstate())
# print("")
# print("")
# print("Validation of function which solve Kepler Equation")
# validation_keplereq()
# print("")
# print("")
# print("Validation of function which return Julian time passed")
# validation_juliantime()
# print("")
# print("")
# print("Validation of function which returns Planets Elements")
# validation_kepelements_equation()
# print("")
# print("")
# print("Validation of function which returns planet position and velocity")
# print(validation_planet_state())
# print("")
# print("")
# print("Validation of funtcion which solve Lamberts problem")
# print(validation_lambert())
# print("")
# print("")
# print("Validation of function which returns a trajectory segment of the mission")
# print(validation_segment_trajectory())
'''
r_Earth, v_Earth = planet_state(Earth, 2457300.0, 0)
R_Mars, v_Mars = planet_state(Mars, 2457300.0, 235)
r_Jupiter, v_Jupiter = planet_state(Jupiter, 2457300.0, 695 + 235)

print("Earth", r_Earth, v_Earth)
print("Mars", R_Mars, v_Mars)
print("Jupiter", r_Jupiter, v_Jupiter)

V_depEarth, V_arrMars = lambert_problem(r_Earth, R_Mars, 235)
V_depMars, V_arrJup = lambert_problem(R_Mars, r_Jupiter, 695)

print(V_depEarth)
print(V_arrMars)
print(V_depMars)
print(V_arrJup)

DV_depEarth = V_depEarth - v_Earth
DV_arrMars = V_arrMars - v_Mars
DV_depMars = V_depMars - v_Mars
DV_arrJupiter = V_arrJup - v_Jupiter

print("Planetocentrica")
print(DV_depEarth)
print(DV_arrMars)
print(DV_depMars)
print(DV_arrJupiter)

print(np.linalg.norm(DV_depEarth))
print(np.linalg.norm(DV_arrMars))
print(np.linalg.norm(DV_depMars))
print(np.linalg.norm(DV_arrJupiter))

fi, r = r_periapse(DV_arrMars, DV_depMars, K_Mars)

print(fi, r)

R_inf = [1.01767010*10**8, 2.04495096*10**8, 2.930510*10**6]
rot, V = fly_by(R_inf, DV_arrMars, v_Mars, K_Mars, r_Mars)
print("rotated angle", "v salida")
print(rot, V)
'''

