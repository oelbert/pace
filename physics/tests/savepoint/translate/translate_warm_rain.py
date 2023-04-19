import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.warm_rain import WarmRain
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class TranslateWarmRain(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "wr_qv"},
            "qliquid": {"serialname": "wr_ql"},
            "qrain": {"serialname": "wr_qr"},
            "qice": {"serialname": "wr_qi"},
            "qsnow": {"serialname": "wr_qs"},
            "qgraupel": {"serialname": "wr_qg"},
            "temperature": {"serialname": "wr_pt"},
            "delp": {"serialname": "wr_delp"},
            "delz": {"serialname": "wr_delz"},
            "density": {"serialname": "wr_den"},
            "density_factor": {"serialname": "wr_denfac"},
            "vterminal_water": {"serialname": "wr_vtw"},
            "vterminal_rain": {"serialname": "wr_vtr"},
            "cloud_condensation_nuclei": {"serialname": "wr_ccn"},
            "reevap": {"serialname": "wr_reevap"},
            "h_var": {"serialname": "wr_h_var"},
        }

        self.in_vars["parameters"] = [
            "dt",
        ]

        self.out_vars = {
            "qvapor": {"serialname": "wr_qv", "kend": namelist.npz - 1},
            "qliquid": {"serialname": "wr_ql", "kend": namelist.npz - 1},
            "qrain": {"serialname": "wr_qr", "kend": namelist.npz - 1},
            "qice": {"serialname": "wr_qi", "kend": namelist.npz - 1},
            "qsnow": {"serialname": "wr_qs", "kend": namelist.npz - 1},
            "qgraupel": {"serialname": "wr_qg", "kend": namelist.npz - 1},
            "temperature": {"serialname": "wr_pt", "kend": namelist.npz - 1},
            "cloud_condensation_nuclei": {
                "serialname": "wr_ccn",
                "kend": namelist.npz - 1,
            },
            "reevap": {"serialname": "wr_reevap", "kend": namelist.npz - 1},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = WarmRain(
            self.stencil_factory,
            self.config,
            timestep=inputs.pop("dt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
