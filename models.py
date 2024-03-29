import fenics as f
import festim as F
from mshr import Rectangle, generate_mesh
from cyl_classes import *
from flibe_props import *


def mesh_2d():
    x1 = 0.00476
    x2 = 0.0162
    y1 = 0.01
    y2 = 0.0761

    p1 = f.Point(0, 0)
    p2 = f.Point(x1, y1)
    r1 = Rectangle(p1, p2)
    p1 = f.Point(x1, 0)
    p2 = f.Point(x1 + x2, y2)
    r2 = Rectangle(p1, p2)
    domain = r1 + r2
    mesh_fenics = generate_mesh(domain, 70)

    f.plot(mesh_fenics)

    # marking physical groups (volumes and surfaces)
    volume_markers = f.MeshFunction("size_t", mesh_fenics, mesh_fenics.topology().dim())
    volume_markers.set_all(1)

    left_surface_str = f"on_boundary && near(x[0], 0, tol)"
    left_surface = f.CompiledSubDomain(left_surface_str, tol=1e-14)

    right_surface_str = f"on_boundary && near(x[0], {x1 + x2}, tol)"
    right_surface = f.CompiledSubDomain(right_surface_str, tol=1e-14)

    bottom_surface_str = f"on_boundary && near(x[1], 0, tol)"
    bottom_surface = f.CompiledSubDomain(bottom_surface_str, tol=1e-14)

    top_surface_str = f"on_boundary && near(x[1], {y2}, tol)"
    top_surface = f.CompiledSubDomain(top_surface_str, tol=1e-14)

    upper_left_surface_str = f"on_boundary && near(x[0], {x1}, tol)"
    upper_left_surface = f.CompiledSubDomain(upper_left_surface_str, tol=1e-14)

    left_top_surface_str = f"on_boundary && near(x[1], {y1}, tol)"
    left_top_surface = f.CompiledSubDomain(left_top_surface_str, tol=1e-14)

    surface_markers = f.MeshFunction(
        "size_t", mesh_fenics, mesh_fenics.topology().dim() - 1
    )
    surface_markers.set_all(0)

    # Surface ids
    left_id = 1
    top_id = 2
    right_id = 3
    bottom_id = 4
    upper_left_id = 5
    left_top_id = 6
    left_surface.mark(surface_markers, left_id)
    right_surface.mark(surface_markers, right_id)
    top_surface.mark(surface_markers, top_id)
    bottom_surface.mark(surface_markers, bottom_id)
    upper_left_surface.mark(surface_markers, upper_left_id)
    left_top_surface.mark(surface_markers, left_top_id)

    f.plot(surface_markers, title="Surface Markers")
    correspondance_dict = {
        "left": left_id,
        "top": top_id,
        "right": right_id,
        "bottom": bottom_id,
        "upper_left": upper_left_id,
        "left_top": left_top_id,
    }
    return mesh_fenics, volume_markers, surface_markers, correspondance_dict


def velocity_field(T_cold, T_hot, my_mesh, surface_markers, correspondance_dict):
    """Computes the velocity field for a given mesh and temperature difference

    Args:
        T_cold (float): the cold temperature (right) in K
        T_hot (float): the hot temperature (left) in K
        my_mesh (fenics.Mesh): the mesh

    Returns:
        fenics.Function, fenics.Function, fenics.Function: velocity field (m/s), pressure (Pa), temperature (K)
    """
    T_bulk = ((T_hot - T_cold) / 2) + T_cold

    V_ele = f.VectorElement("CG", my_mesh.ufl_cell(), 2)
    Q_ele = f.FiniteElement("CG", my_mesh.ufl_cell(), 1)
    T_ele = f.FiniteElement("CG", my_mesh.ufl_cell(), 1)
    W = f.FunctionSpace(my_mesh, f.MixedElement([V_ele, Q_ele, T_ele]))

    x = f.SpatialCoordinate(my_mesh)

    def div_cyl(v, x):
        return (1 / x[0]) * (x[0] * v[0]).dx(0) + v[1].dx(1)

    upT = f.Function(W)
    upT_old = f.Function(W)
    u, p, T = f.split(upT)
    v, q, S = f.TestFunctions(W)

    # Cylindrical Coordinates
    for factor in [1e-03, 1e-02, 1e-01, 1]:
        print("Running for factor={:.1e}".format(factor))

        g = f.Constant((0, -9.81))  # gravity acceleration in m/s2
        mu = viscosity_flibe1(T_bulk)  # dynamic viscosity in kg/m/s
        rho = density_flibe1(T_bulk)  # density in kg/m3
        rho_0 = density_flibe1(T_cold)  # density at T_cold
        cp = 2386  # heat capacity in J/(kg.K)
        thermal_cond = 1.1  # thermal conductivity in W/(m.K)
        beta = beta_flibe1(T_bulk) * factor

        # CFD momentum
        F = (
            rho_0 * f.inner(f.dot(f.grad(u), u), v) * x[0] * f.dx
            - f.inner(p, div_cyl(v, x)) * x[0] * f.dx
            + mu * f.inner(f.grad(u), f.grad(v)) * x[0] * f.dx
            + f.inner(rho_0 * beta * (T - T_bulk) * g, v) * x[0] * f.dx
        )

        # CFD continuity
        F -= f.inner(q, div_cyl(u, x)) * x[0] * f.dx

        # Heat transfer
        F += rho * cp * f.inner(f.dot(f.grad(T), u), S) * x[0] * f.dx
        F += f.inner(thermal_cond * f.grad(T), f.grad(S)) * x[0] * f.dx

        # Temperature boundary conditions
        upper_left_surface_TBC = f.DirichletBC(
            W.sub(2), T_hot, surface_markers, correspondance_dict["upper_left"]
        )
        left_top_surface_TBC = f.DirichletBC(
            W.sub(2), T_hot, surface_markers, correspondance_dict["left_top"]
        )
        right_surface_TBC = f.DirichletBC(
            W.sub(2), T_cold, surface_markers, correspondance_dict["right"]
        )

        # Non-slip NS boundary conditions
        upper_left_surface_VBC = f.DirichletBC(
            W.sub(0),
            f.Constant((0, 0)),
            surface_markers,
            correspondance_dict["upper_left"],
        )
        left_top_surface_VBC = f.DirichletBC(
            W.sub(0),
            f.Constant((0, 0)),
            surface_markers,
            correspondance_dict["left_top"],
        )
        right_surface_VBC = f.DirichletBC(
            W.sub(0), f.Constant((0, 0)), surface_markers, correspondance_dict["right"]
        )
        bottom_surface_VBC = f.DirichletBC(
            W.sub(0), f.Constant((0, 0)), surface_markers, correspondance_dict["bottom"]
        )

        # Non-slip through boundary conditions
        left_surface_VBC = f.DirichletBC(
            W.sub(0).sub(0), f.Constant(0), surface_markers, correspondance_dict["left"]
        )
        top_surface_VBC = f.DirichletBC(
            W.sub(0).sub(1), f.Constant(0), surface_markers, correspondance_dict["top"]
        )

        bcs = [
            upper_left_surface_TBC,
            left_top_surface_TBC,
            right_surface_TBC,
            upper_left_surface_VBC,
            left_top_surface_VBC,
            right_surface_VBC,
            bottom_surface_VBC,
            left_surface_VBC,
            top_surface_VBC,
        ]

        # # Old boundary conditions
        # bcs = [
        #     f.DirichletBC(W.sub(0), Constant((0, 0)), "on_boundary"),
        #     f.DirichletBC(W.sub(2), T_hot, "on_boundary && x[0] == 0"),
        #     f.DirichletBC(W.sub(2), T_cold, "on_boundary && x[0] == 0.0162"),
        # ]

        f.solve(
            F == 0,
            upT,
            bcs=bcs,
            solver_parameters={
                "newton_solver": {
                    "linear_solver": "mumps",
                    "absolute_tolerance": 1e-09,
                    "relative_tolerance": 1e-09,
                    "maximum_iterations": 25,
                }
            },
        )

        upT_old.assign(upT)

    u, p, T = upT.split()

    return u, p, T


def t_transport_sim(
    temperature_field,
    mesh_fenics,
    velocity,
    volume_markers,
    surface_markers,
    correspondance_dict,
):
    """
    Takes in a list of temperatures and a set mesh and returns a list of diffusion coefficients that correspond to each temperature

    Args:
        temperature_field (fenics.Function): the temperature field in K
        mesh_fenics (fenics.Mesh): the mesh (should be the same as the one used to compute the temperature field)
        velocity (fenics.Function): the velocity field in m/s

    Returns:
        float: the mass transport coefficient in m/s
    """
    # setting up current simulation
    model_2d = F.Simulation()

    # D, E_d source: "nakamura_hydrogen_2015"
    # Thermal cond source: https://dspace.mit.edu/bitstream/handle/1721.1/123988/Romatoski_SaltPropertyReview02.pdf?sequence=1&isAllowed=y#:~:text=From%20the%20data%20collected%2C%20the,an%20uncertainty%20of%20%C2%B110%25.
    flibe_mat = F.Material(
        id=1,
        D_0=1.508521565198744e-08,
        E_D=0.23690444592353738,
    )
    model_2d.materials = F.Materials([flibe_mat])

    # creating mesh with festim
    model_2d.mesh = F.Mesh(
        mesh=mesh_fenics,  # TODO we should be able to get the mesh from the temperature field
        volume_markers=volume_markers,
        surface_markers=surface_markers,
    )

    # setting up steady state heat transfer problem

    # model_2d.T = F.TemperatureFromXDMF(temperature_file, label="temperature")
    model_2d.T = F.Temperature(
        value=973
    )  # dummy temperature, will be overwritten later

    # setting up T source
    rthetaz = f.SpatialCoordinate(mesh_fenics)
    salt_volume = 2 * np.pi * f.assemble(rthetaz[0] * f.dx())
    measured_tritium_source = 1.84e5  # T/s

    model_2d.sources = [
        F.Source(value=measured_tritium_source / salt_volume, volume=1, field=0)
    ]

    top_id = correspondance_dict["top"]
    bottom_id = correspondance_dict["bottom"]
    right_id = correspondance_dict["right"]
    left_id = correspondance_dict["left"]
    upper_left_id = correspondance_dict["upper_left"]
    left_top_id = correspondance_dict["left_top"]

    # setting up transport boundary conditions
    tritium_transport_bcs = [
        F.DirichletBC(
            surfaces=[top_id, bottom_id, right_id, left_top_id, left_id, upper_left_id],
            value=0,
            field=0,
        )
    ]

    model_2d.boundary_conditions = tritium_transport_bcs

    # simulation parameters and running model
    model_2d.settings = F.Settings(
        transient=False,
        absolute_tolerance=1e-09,
        relative_tolerance=1e-09,
    )

    # setting up exports
    export_folder = "BABY_2D_results"

    derived_quantities = F.DerivedQuantities(filename=export_folder + "/simulation.csv")

    derived_quantities.derived_quantities = [
        SurfaceFluxCylindrical(field="solute", surface=right_id),
        SurfaceFluxCylindrical(field="solute", surface=top_id),
        SurfaceFluxCylindrical(field="solute", surface=bottom_id),
        SurfaceFluxCylindrical(field="solute", surface=upper_left_id),
        SurfaceFluxCylindrical(field="solute", surface=left_id),
        SurfaceFluxCylindrical(field="solute", surface=left_top_id),
        AverageVolumeCylindrical(field="solute", volume=1),
    ]

    model_2d.exports = F.Exports(
        [
            F.XDMFExport("solute", folder=export_folder),
            F.XDMFExport("retention", folder=export_folder),
            F.XDMFExport("T", folder=export_folder),
            derived_quantities,
        ]
    )
    # adding advection
    model_2d.initialise()  # reinitialisation is needed

    model_2d.T.T = temperature_field
    # model_2d.T.T_n = temperature_field  # don't know if this is needed

    hydrogen_concentration = model_2d.h_transport_problem.mobile.solution
    test_function_solute = model_2d.h_transport_problem.mobile.test_function

    advection_term = (
        f.inner(f.dot(f.grad(hydrogen_concentration), velocity), test_function_solute)
        * model_2d.mesh.dx
    )

    model_2d.h_transport_problem.F += advection_term

    model_2d.run()

    plt.figure()
    plt.title("Hydrogen concentration")
    CS = f.plot(hydrogen_concentration)
    # f.plot(velocity, scale=1e-3, color="black", alpha=0.5)
    plt.colorbar(CS, label="H/m3")
    plt.show()

    # reading results
    my_data = np.genfromtxt(
        export_folder + "/simulation.csv", names=True, delimiter=","
    )

    flux_1 = my_data["Flux_surface_1_solute"]
    flux_2 = my_data["Flux_surface_2_solute"]
    flux_3 = my_data["Flux_surface_3_solute"]
    flux_4 = my_data["Flux_surface_4_solute"]
    flux_5 = my_data["Flux_surface_5_solute"]
    flux_6 = my_data["Flux_surface_6_solute"]

    # calculating diffusion coefficient
    total_flux = abs(flux_1 + flux_2 + flux_3 + flux_4 + flux_5 + flux_6)

    average_conc = my_data["Average_solute_volume_1"]

    total_surface = 2 * np.pi * f.assemble(rthetaz[0] * model_2d.mesh.ds)
    print(f"Total surface: {total_surface:.2e} m2")
    k = total_flux / (total_surface * average_conc)

    print(f"Total flux: {total_flux:.2e} H/s/m")
    print(f"Average concentration: {average_conc:.2e} H/m3")
    print(f"k: {k:.2e} m/s")

    return k
