from mshr import Rectangle, generate_mesh
from fenics import *
import fencis as f
import mshr

def rectangleMeshFenics():
    # creating a mesh with FEniCS
    nx = round(50*(0.0162/0.0761))
    ny = 50
    nx2 = round(50*(0.00476/0.0761))
    ny2 = round(50*(0.01/0.0761))
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
    Generates backwards shaped L mesh using mshr
    Parameters:
        x1: int, lower left corner
        x2: int, location where L stem begins from origin
        y1: int, height of lower L piece
        y2: int, height of L stem
        resolution: int, specific to resolution of mesh
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


    plot(surface_markers, title="Surface Markers")

    return 



# def LMeshgmsh():
