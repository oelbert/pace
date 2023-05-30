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
            "qvapor": {"serialname": "szs_qv", "mp3": True},
            "qliquid": {"serialname": "szs_ql", "mp3": True},
            "qrain": {"serialname": "szs_qr", "mp3": True},
            "qice": {"serialname": "szs_qi", "mp3": True},
            "qsnow": {"serialname": "szs_qs", "mp3": True},
            "qgraupel": {"serialname": "szs_qg", "mp3": True},
            "temperature": {"serialname": "szs_pt", "mp3": True},
            "density": {"serialname": "szs_den", "mp3": True},
            "density_factor": {"serialname": "szs_denfac", "mp3": True},
            "delp": {"serialname": "szs_delp", "mp3": True},
            "rh_adj": {"serialname": "szs_rh_adj", "mp3": True},
            "cloud_condensation_nuclei": {"serialname": "szs_ccn", "mp3": True},
            "cloud_ice_nuclei": {"serialname": "szs_cin", "mp3": True},
            "cond": {"serialname": "szs_cond", "mp3": True},
            "dep": {"serialname": "szs_dep", "mp3": True},
            "reevap": {"serialname": "szs_reevap", "mp3": True},
            "sub": {"serialname": "szs_sub", "mp3": True},
        }

        self.in_vars["parameters"] = [
            "dt",
        ]

        self.out_vars = {
            "qvapor": {"serialname": "szs_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "szs_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "szs_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "szs_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "szs_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "szs_qg", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "szs_pt", "kend": namelist.npz, "mp3": True},
            "cloud_condensation_nuclei": {
                "serialname": "szs_ccn",
                "kend": namelist.npz,
                "mp3": True,
            },
            "cloud_ice_nuclei": {
                "serialname": "szs_cin",
                "kend": namelist.npz,
                "mp3": True,
            },
            "cond": {"serialname": "szs_cond", "kend": namelist.npz, "mp3": True},
            "dep": {"serialname": "szs_dep", "kend": namelist.npz, "mp3": True},
            "reevap": {"serialname": "szs_reevap", "kend": namelist.npz, "mp3": True},
            "sub": {"serialname": "szs_sub", "kend": namelist.npz, "mp3": True},
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
