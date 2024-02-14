from mshr import Rectangle, generate_mesh
from fenics import *
import fencis as f
import mshr
import gmsh

def rectMeshFenics():
    # creating a mesh with FEniCS
    nx = round(50*(0.0162/0.0761))
    ny = 50
    # nx2 = round(50*(0.00476/0.0761))
    # ny2 = round(50*(0.01/0.0761))
    # mesh_fenics = UnitSquareMesh(nx, ny)

    mesh_fenics = RectangleMesh(Point(0.00476, 0.0), Point(0.00476 + 0.0162, 0.0761), nx, ny)

    plot(mesh_fenics)

    # marking physical groups (volumes and surfaces)
    volume_markers = MeshFunction(
        "size_t", mesh_fenics, mesh_fenics.topology().dim())
    volume_markers.set_all(1)

    left_surface = CompiledSubDomain(
        'on_boundary && near(x[0], 0, tol)', tol=1e-14)
    right_surface = CompiledSubDomain(
        'on_boundary && near(x[0], 0.0162, tol)', tol=1e-14)
    bottom_surface = CompiledSubDomain(
        'on_boundary && near(x[1], 0, tol)', tol=1e-14)
    top_surface = CompiledSubDomain(
        'on_boundary && near(x[1], 0.0761, tol)', tol=1e-14)

    surface_markers = MeshFunction(
        "size_t", mesh_fenics, mesh_fenics.topology().dim() - 1)
    surface_markers.set_all(0)

    # Surface ids
    left_id = 1
    top_id = 2
    right_id = 3
    bottom_id = 4
    left_surface.mark(surface_markers, left_id)
    right_surface.mark(surface_markers, right_id)
    top_surface.mark(surface_markers, top_id)
    bottom_surface.mark(surface_markers, bottom_id)

    plot(mesh_fenics)



def LMeshmshr(x1, x2, y1, y2, resolution):
    """
    Creates backwards L shaped mesh
    Parameters:
        x1: float, location of staff from origin
        x2: float, location of end of staff
        y1: float, height of bottom part of L
        y2: float, height of staff
        resolution: int, density of mesh
    """
    p1 = f.Point(0, 0)
    p2 = f.Point(x1, y1)
    r1 = Rectangle(p1, p2)
    p1 = f.Point(x1, 0)
    p2 = f.Point(x1 + x2, y2)
    r2 = Rectangle(p1, p2)
    domain = r1 + r2
    mesh_fenics = generate_mesh(domain, resolution)

    plot(mesh_fenics)

    # marking physical groups (volumes and surfaces)
    volume_markers = MeshFunction(
        "size_t", mesh_fenics, mesh_fenics.topology().dim())
    volume_markers.set_all(1)

    left_surface_str = f'on_boundary && near(x[0], {x1}, tol)'
    left_surface = CompiledSubDomain(
        left_surface_str, tol=1e-14)

    right_surface_str = f'on_boundary && near(x[0], {x1 + x2}, tol)'
    right_surface = CompiledSubDomain(
        right_surface_str, tol=1e-14)

    bottom_surface_str = f'on_boundary && near(x[1], 0, tol)'
    bottom_surface = CompiledSubDomain(
        bottom_surface_str, tol=1e-14)

    top_surface_str = f'on_boundary && near(x[1], {y2}, tol)'
    top_surface = CompiledSubDomain(
        top_surface_str, tol=1e-14)

    upper_left_surface_str = f'on_boundary && near(x[0], {x2}, tol)'
    upper_left_surface = CompiledSubDomain(
        upper_left_surface_str, tol=1e-14)

    left_top_surface_str = f'on_boundary && near(x[1], {y2}, tol)'
    left_top_surface = CompiledSubDomain(
        left_top_surface_str, tol=1e-14)


    surface_markers = MeshFunction(
        "size_t", mesh_fenics, mesh_fenics.topology().dim() - 1)
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

    return mesh_fenics, surface_markers, left_id, top_id, right_id, bottom_id, upper_left_id, left_top_id


def LMeshgmsh():
    # Initialize Gmsh
    gmsh.initialize()

    # Create a new model
    gmsh.model.add("l_shape")

    # Define parameters
    lc = 0.1  # characteristic length

    # Add points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(2, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(2, 1, 0, lc)
    p4 = gmsh.model.geo.addPoint(1, 1, 0, lc)
    p5 = gmsh.model.geo.addPoint(1, 2, 0, lc)
    p6 = gmsh.model.geo.addPoint(0, 2, 0, lc)

    # Add lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p5)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p1)

    # Add line loop
    ll = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6])

    # Add plane surface
    ps = gmsh.model.geo.addPlaneSurface([ll])

    # Physical groups
    gmsh.model.addPhysicalGroup(1, [l1], tag=1)  # boundary 1
    gmsh.model.setPhysicalName(1, 1, "Dirichlet_boundary")
    gmsh.model.addPhysicalGroup(1, [l3], tag=2)  # boundary 2
    gmsh.model.setPhysicalName(1, 2, "Neumann_boundary")
    gmsh.model.addPhysicalGroup(1, [l5], tag=3)  # boundary 3
    gmsh.model.setPhysicalName(1, 3, "Other_boundary")
    gmsh.model.addPhysicalGroup(2, [ps], tag=4)  # interior
    gmsh.model.setPhysicalName(2, 4, "Interior")

    # Generate mesh
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # Get mesh and cleanup Gmsh
    mesh = dolfin.Mesh()
    gmsh_option = "-v 0"
    gmsh.write("l_shape.msh")
    gmsh.finalize()

    # Convert mesh to XML format
    with dolfin.XDMFFile("l_shape.xml") as xdmf:
        xdmf.write(mesh)


if __name__ == "__main__":
    mesh_fenics, surface_markers, left_id, top_id, right_id, bottom_id, upper_left_id, left_top_id = LMeshmshr