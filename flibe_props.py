import numpy as np
import matplotlib.pyplot as plt

DOLFIN_EPS = 1e-16


# densities taken from fluoride salt coolant properties paper
def density_flibe1(T):  # Janz 1974/1988 (3,14,17,21,22)
    return 2413 - 0.488 * T  # kg/m3


def density_flibe2(T):  # Cantor 1968 (18,19)
    return 2214 - 0.42 * T


def density_flibe3(T):  # Zaghloul 2003 (7,14,22)
    return 2415.6 - 0.49072 * T


def density_flibe4(T):  # Ignat'ev et al. 2006
    T = np.asarray(T)
    return np.where(T < 973, 2163 - 0.406 * (T - 601.4), 2163 - 0.687 * (T - 601.4))


def density_flibe5(T):  # Williams et al. 2006
    return 2280 - 0.488 * T


def density_flibe6(T):  # Vidrio et al., 2022
    return 2245 - 0.424 * T


def density_flibe7(T):  # Chapdelaine, 2017
    return 2241.6 - 0.42938 * T


density_prop_array = [
    (density_flibe1, "Janz 1974/1988"),
    (density_flibe2, "Cantor 1968"),
    (density_flibe3, "Zaghloul 2003"),
    (density_flibe4, "Ignat'ev et al. 2006"),
    (density_flibe5, "Williams et al. 2006"),
    (density_flibe6, "Vidrio et al., 2022"),
    (density_flibe7, "Chapdelaine, 2017"),
]


# viscosities taken from fluoride salt coolant properties paper

mPa_to_Pa = 1e-03


def viscosity_flibe1(T):  # Cantor (3–7,11,14,18)
    return mPa_to_Pa * 0.116 * np.exp(3755 / (T + DOLFIN_EPS))  # Pa.s


def viscosity_flibe2(T):  # Romatoski (12,17,21,32)
    return mPa_to_Pa * 0.0594 * np.exp(4605 / (T + DOLFIN_EPS))  # Pa.s


def viscosity_flibe3(T):  # Gierszewski (20)
    return mPa_to_Pa * 0.116 * np.exp(3760 / (T + DOLFIN_EPS))  # Pa.s


def viscosity_flibe4(T):  # Cohen and Jones (30,31,33)
    return mPa_to_Pa * 0.118 * np.exp(3624 / (T + DOLFIN_EPS))  # Pa.s


viscosity_prop_array = [
    (viscosity_flibe1, "Cantor"),
    (viscosity_flibe2, "Romatoski"),
    (viscosity_flibe3, "Gierszewski"),
    (viscosity_flibe4, "Cohen and Jones"),
]


def beta_flibe1(T):
    return 1.8319e-04 + 5.55e-08 * T  # K-1


def beta_flibe2(T):
    return 2.3e-04 + 5.55e-08 * T  # K-1


beta_prop_array = [(beta_flibe1, "1"), (beta_flibe2, "2")]


def heat_capacity_flibe(T):
    return 2386  # J/kg.K


def thermal_conductivity_flibe(T):
    return 1.1  # W/m.K


if __name__ == "__main__":
    T = np.linspace(700, 1200, 1000)

    for viscosity in viscosity_prop_array:
        plt.plot(T, viscosity[0](T), label=viscosity[1])
    plt.legend()
    plt.xlabel("Temperature (K)")
    plt.ylabel("Viscosity (Pa.s)")
    plt.ylim(bottom=0)
    plt.show()

    for density in density_prop_array:
        plt.plot(T, density[0](T), label=density[1])
    plt.legend()
    plt.xlabel("Temperature (K)")
    plt.ylabel("Density (kg/m3)")
    plt.show()

    for beta in beta_prop_array:
        plt.plot(T, beta[0](T), label=beta[1])
    plt.legend()
    plt.xlabel("Temperature (K)")
    plt.ylabel("Thermal expansion coefficient (K-1)")
    plt.show()
