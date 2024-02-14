import festim as F
import numpy as np
import fenics as f

class AverageVolumeCylindrical(F.VolumeQuantity):
    def __init__(self, field, volume: int) -> None:
        super().__init__(field, volume)
        self.title = "Average {} volume {}".format(self.field, self.volume)

    def compute(self):

        mesh = self.function.function_space().mesh()  # get the mesh from the function
        rthetaz = f.SpatialCoordinate(mesh)  # get the coordinates from the mesh
        r = rthetaz[0]  # only care about r here

        return f.assemble(r * self.function * self.dx(self.volume)) / f.assemble(
            r * self.dx(self.volume)
        )

class SurfaceFluxCylindrical(F.SurfaceFlux): # Inherets from class SurfaceFlux
    def __init__(self, field, surface) -> None:
        super().__init__(field, surface)
        self.r = None

    def compute(self, soret=False):
        field_to_prop = {
            "0": self.D,
            "solute": self.D,
            0: self.D,
            "T": self.thermal_cond,
        }
        self.prop = field_to_prop[self.field]
        if soret:
            raise NotImplementedError(
                "Soret effect not implemented for cylindrical coordinates"
            )

        if self.r is None:
            mesh = (
                self.function.function_space().mesh()
            )  # get the mesh from the function
            rthetaz = f.SpatialCoordinate(mesh)  # get the coordinates from the mesh
            self.r = rthetaz[0]  # only care about r here

        # dS_z = r dr dtheta , assuming axisymmetry dS_z = theta r dr
        # dS_r = r dz dtheta , assuming axisymmetry dS_r = theta r dz
        # in both cases the expression with self.ds is the same
        # we assume full cylinder theta = 2 pi
        flux = f.assemble(
            self.prop
            * self.r
            * f.dot(f.grad(self.function), self.n)
            * self.ds(self.surface)
        )
        theta = 2 * np.pi
        flux *= theta
        return flux