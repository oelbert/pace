import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.ice_cloud import IceCloud
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class TranslateIceCloud(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "ic_qv", "microph": True},
            "qliquid": {"serialname": "ic_ql", "microph": True},
            "qrain": {"serialname": "ic_qr", "microph": True},
            "qice": {"serialname": "ic_qi", "microph": True},
            "qsnow": {"serialname": "ic_qs", "microph": True},
            "qgraupel": {"serialname": "ic_qg", "microph": True},
            "temperature": {"serialname": "ic_pt", "microph": True},
            "density": {"serialname": "ic_den", "microph": True},
            "density_factor": {"serialname": "ic_denfac", "microph": True},
            "vterminal_water": {"serialname": "ic_vtw", "microph": True},
            "vterminal_rain": {"serialname": "ic_vtr", "microph": True},
            "vterminal_ice": {"serialname": "ic_vti", "microph": True},
            "vterminal_snow": {"serialname": "ic_vts", "microph": True},
            "vterminal_graupel": {"serialname": "ic_vtg", "microph": True},
            "h_var": {"serialname": "ic_h_var", "microph": True},
        }

        self.in_vars["parameters"] = [
            "dt",
        ]

        self.out_vars = {
            "qvapor": {"serialname": "ic_qv", "kend": namelist.npz - 1},
            "qliquid": {"serialname": "ic_ql", "kend": namelist.npz - 1},
            "qrain": {"serialname": "ic_qr", "kend": namelist.npz - 1},
            "qice": {"serialname": "ic_qi", "kend": namelist.npz - 1},
            "qsnow": {"serialname": "ic_qs", "kend": namelist.npz - 1},
            "qgraupel": {"serialname": "ic_qg", "kend": namelist.npz - 1},
            "temperature": {"serialname": "ic_pt", "kend": namelist.npz - 1},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = IceCloud(
            self.stencil_factory,
            self.config,
            timestep=inputs.pop("dt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
