import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.subgrid_z_proc import (
    VerticalSubgridProcesses,
)
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class TranslateSubgridZProc(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "sz_qv"},
            "qliquid": {"serialname": "sz_ql"},
            "qrain": {"serialname": "sz_qr"},
            "qice": {"serialname": "sz_qi"},
            "qsnow": {"serialname": "sz_qs"},
            "qgraupel": {"serialname": "sz_qg"},
            "temperature": {"serialname": "sz_pt"},
            "density": {"serialname": "sz_den"},
            "density_factor": {"serialname": "sz_denfac"},
            "delp": {"serialname": "wr_delp"},
            "rh_adj": {"serialname": "sz_rh_adj"},
            "cloud_condensation_nuclei": {"serialname": "sz_ccn"},
            "cloud_ice_nuclei": {"serialname": "sz_cin"},
            "condensation": {"serialname": "sz_cond"},
            "deposition": {"serialname": "sz_dep"},
            "evaporation": {"serialname": "sz_evap"},
            "sublimation": {"serialname": "sz_sub"},
        }

        self.in_vars["parameters"] = [
            "dt",
        ]

        self.out_vars = {
            "qvapor": {"serialname": "sz_qv", "kend": namelist.npz - 1},
            "qliquid": {"serialname": "sz_ql", "kend": namelist.npz - 1},
            "qrain": {"serialname": "sz_qr", "kend": namelist.npz - 1},
            "qice": {"serialname": "sz_qi", "kend": namelist.npz - 1},
            "qsnow": {"serialname": "sz_qs", "kend": namelist.npz - 1},
            "qgraupel": {"serialname": "sz_qg", "kend": namelist.npz - 1},
            "temperature": {"serialname": "sz_pt", "kend": namelist.npz - 1},
            "cloud_condensation_nuclei": {
                "serialname": "sz_ccn",
                "kend": namelist.npz - 1,
            },
            "cloud_ice_nuclei": {"serialname": "sz_cin", "kend": namelist.npz - 1},
            "condensation": {"serialname": "sz_cond", "kend": namelist.npz - 1},
            "deposition": {"serialname": "sz_dep", "kend": namelist.npz - 1},
            "evaporation": {"serialname": "sz_evap", "kend": namelist.npz - 1},
            "sublimation": {"serialname": "sz_sub", "kend": namelist.npz - 1},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = VerticalSubgridProcesses(
            self.stencil_factory,
            self.config,
            timestep=inputs.pop("dt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
