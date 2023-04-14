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
            "qvapor": {"serialname": "sz_qv", "microph": True},
            "qliquid": {"serialname": "sz_ql", "microph": True},
            "qrain": {"serialname": "sz_qr", "microph": True},
            "qice": {"serialname": "sz_qi", "microph": True},
            "qsnow": {"serialname": "sz_qs", "microph": True},
            "qgraupel": {"serialname": "sz_qg", "microph": True},
            "temperature": {"serialname": "sz_pt", "microph": True},
            "density": {"serialname": "sz_den", "microph": True},
            "density_factor": {"serialname": "sz_denfac", "microph": True},
            "delp": {"serialname": "wr_delp", "microph": True},
            "rh_adj": {"serialname": "sz_rh_adj", "microph": True},
            "cloud_condensation_nuclei": {"serialname": "sz_ccn", "microph": True},
            "cloud_ice_nuclei": {"serialname": "sz_cin", "microph": True},
            "condensation": {"serialname": "sz_cond", "microph": True},
            "deposition": {"serialname": "sz_dep", "microph": True},
            "evaporation": {"serialname": "sz_evap", "microph": True},
            "sublimation": {"serialname": "sz_sub", "microph": True},
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
