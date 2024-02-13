from fencis import *
import numpy as np


# densities taken from fluoride salt coolant properties paper
def density_flibe1(T): # Janz 1974/1988 (3,14,17,21,22)
    return 2413 - 0.488 * T  # kg/m3

def density_flibe2(T): # Cantor 1968 (18,19)
    return 2214 - 0.42 * T

def density_flibe3(T): # Zaghloul 2003 (7,14,22)
    return 2415.6 - 0.49072 * T

def density_flibe4(T): # Ignat'ev et al. 2006
    if T < 973:
        return 2163 - 0.406*(T-601.4)
    else:
        return 2163 - 0.687*(T-601.4)

def density_flibe5(T): # Williams et al. 2006
    return 2280 - 0.488 * T

def density_flibe6(T): # Vidrio et al., 2022
    return 2245 - 0.424 * T

def density_flibe7(T): # Chapdelaine, 2017
    return 2241.6 - 0.42938 * T

density_prop_array = [(density_flibe1, "Janz 1974/1988"),
                      (density_flibe2, "Cantor 1968"),
                      (density_flibe3, "Zaghloul 2003"),
                      (density_flibe4, "Ignat'ev et al. 2006"),
                      (density_flibe5, "Williams et al. 2006"),
                      (density_flibe6, "Vidrio et al., 2022"),
                      (density_flibe7, "Chapdelaine, 2017")]


# viscosities taken from fluoride salt coolant properties paper
e = 2.718281828459045
def viscosity_flibe1(T): # Cantor (3–7,11,14,18)
    return 1.16e-04 * e ** (3755 / (T + DOLFIN_EPS))  # Pa.s

def viscosity_flibe2(T): # Romatoski (12,17,21,32)
    return 5.94e-05 * e ** (4605 / (T + DOLFIN_EPS))  # Pa.s

def viscosity_flibe3(T): # Gierszewski (20)
    return 1.16e-04 * e ** (3760 / (T + DOLFIN_EPS))  # Pa.s

def viscosity_flibe4(T): # Cohen and Jones (30,31,33)
    return 1.18e-04 * e ** (3624 / (T + DOLFIN_EPS))  # Pa.s

viscosity_prop_array = [(viscosity_flibe1, "Cantor"),
                        (viscosity_flibe2, "Romatoski"),
                        (viscosity_flibe3, "Gierszewski"),
                        (viscosity_flibe4, "Cohen and Jones")]


def beta_flibe1(T):
    return 1.8319e-04 + 5.55e-08 * T  # K-1

def beta_flibe2(T):
    return 2.3e-04 + 5.55e-08 * T  # K-1

beta_prop_array = [(beta_flibe1, '1'),
                   (beta_flibe2, '2')]