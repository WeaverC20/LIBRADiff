from fenics import *
import matplotlib.pyplot as plt
import numpy as np


def density_flibe(T):
    return 2413 - 0.488 * T  # kg/m3


def viscosity_flibe(T):
    e = 2.718281828459045
    return 1.16e-04 * e ** (3755 / (T + DOLFIN_EPS))  # Pa.s


def beta_flibe(T):
    return 1.8319e-04 + 5.55e-08 * T  # K-1


def plot_properties():
    T = np.linspace(400, 1200, 1000)

    mu = viscosity_flibe(T)
    rho = density_flibe(T)
    beta = beta_flibe(T)

    plt.figure()
    mu *= 1e03
    plt.plot(T, mu, label="viscosity")
    plt.ylim(0, 30)
    plt.legend()

    plt.figure()
    plt.plot(T, rho, label="density")
    plt.ylim(1800, 2150)
    plt.legend()

    plt.figure()
    beta *= 1e04
    plt.plot(T, beta, label="Thermal expansion")
    plt.ylim(1.6, 2.8)
    plt.legend()

    plt.show()


def velocity_field(T_cold, T_hot, results_folder="Results/"):
    T_bulk = ((T_hot - T_cold) / 2) + T_cold

    N = 30
    my_mesh = RectangleMesh(Point(0, 0), Point(0.1, 0.1), N, N)
    V_ele = VectorElement("CG", my_mesh.ufl_cell(), 2)
    Q_ele = FiniteElement("CG", my_mesh.ufl_cell(), 1)
    T_ele = FiniteElement("CG", my_mesh.ufl_cell(), 1)
    W = FunctionSpace(my_mesh, MixedElement([V_ele, Q_ele, T_ele]))

    upT = Function(W)
    upT_old = Function(W)
    u, p, T = split(upT)
    v, q, S = TestFunctions(W)

    for factor in [1e-03, 1e-02, 1e-01, 1]:
        print("Running for factor={:.1e}".format(factor))

        g = Constant((0, -9.81))  # gravity acceleration in m/s2
        mu = viscosity_flibe(T_bulk)  # dynamic viscosity in kg/m/s
        rho = density_flibe(T_bulk)  # density in kg/m3
        rho_0 = density_flibe(T_bulk)  # density at T_cold
        cp = 2386  # heat capacity in J/(kg.K)
        thermal_cond = 1.1  # thermal conductivity in W/(m.K)
        beta = beta_flibe(T_bulk) * factor

        # CFD momentum
        F = (
            rho_0 * inner(dot(grad(u), u), v) * dx
            - inner(p, div(v)) * dx
            + mu * inner(grad(u), grad(v)) * dx
            + inner(rho_0 * beta * (T - T_bulk) * g, v) * dx
        )

        # CFD continuity
        F -= inner(q, div(u)) * dx

        # Heat transfer
        F += rho * cp * inner(dot(grad(T), u), S) * dx
        F += inner(thermal_cond * grad(T), grad(S)) * dx

        bcs = [
            DirichletBC(W.sub(0), Constant((0, 0)), "on_boundary"),
            DirichletBC(W.sub(2), T_hot, "on_boundary && x[0] == 0"),
            DirichletBC(W.sub(2), T_cold, "on_boundary && x[0] == 0.1"),
        ]

        solve(
            F == 0,
            upT,
            bcs=bcs,
            solver_parameters={
                "newton_solver": {
                    "linear_solver": "mumps",
                    "absolute_tolerance": 1e-09,
                    "relative_tolerance": 1e-09,
                    "maximum_iterations": 15,
                }
            },
        )

        upT_old.assign(upT)

    u, p, T = upT.split()

    XDMFFile(results_folder + "velocity_field.xdmf").write(u)
    XDMFFile(results_folder + "temperature_field.xdmf").write(T)


if __name__ == "__main__":
    # plot_properties()

    T_cold_temps = [1, 3, 5, 10, 20]
    T_hot = 700 + 273.15
    for dt in T_cold_temps:
        T_cold = T_hot - dt
        foldername = f"Results/dt={dt}/"
        velocity_field(T_cold=T_cold, T_hot=T_hot, results_folder=foldername)
