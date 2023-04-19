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
            "qvapor": {"serialname": "ic_qv"},
            "qliquid": {"serialname": "ic_ql"},
            "qrain": {"serialname": "ic_qr"},
            "qice": {"serialname": "ic_qi"},
            "qsnow": {"serialname": "ic_qs"},
            "qgraupel": {"serialname": "ic_qg"},
            "temperature": {"serialname": "ic_pt"},
            "density": {"serialname": "ic_den"},
            "density_factor": {"serialname": "ic_denfac"},
            "vterminal_water": {"serialname": "ic_vtw"},
            "vterminal_rain": {"serialname": "ic_vtr"},
            "vterminal_ice": {"serialname": "ic_vti"},
            "vterminal_snow": {"serialname": "ic_vts"},
            "vterminal_graupel": {"serialname": "ic_vtg"},
            "h_var": {"serialname": "ic_h_var"},
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
