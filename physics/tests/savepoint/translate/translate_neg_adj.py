import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.neg_adj import AdjustNegativeTracers
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class TranslateNegAdjP(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "ne_qv", "microph": True},
            "qliquid": {"serialname": "ne_ql", "microph": True},
            "qrain": {"serialname": "ne_qr", "microph": True},
            "qice": {"serialname": "ne_qi", "microph": True},
            "qsnow": {"serialname": "ne_qs", "microph": True},
            "qgraupel": {"serialname": "ne_qg", "microph": True},
            "temperature": {"serialname": "ne_pt", "microph": True},
            "delp": {"serialname": "ne_delp", "microph": True},
            "condensation": {"serialname": "ne_cond", "microph": True},
        }

        self.in_vars["parameters"] = [
            "convt"
        ]

        self.out_vars = {
            "qvapor": {"serialname": "ne_qv", "kend": namelist.npz - 1},
            "qliquid": {"serialname": "ne_ql", "kend": namelist.npz - 1},
            "qrain": {"serialname": "ne_qr", "kend": namelist.npz - 1},
            "qice": {"serialname": "ne_qi", "kend": namelist.npz - 1},
            "qsnow": {"serialname": "ne_qs", "kend": namelist.npz - 1},
            "qgraupel": {"serialname": "ne_qg", "kend": namelist.npz - 1},
            "pt": {"serialname": "ne_pt", "kend": namelist.npz - 1},
            "delp": {"serialname": "ne_delp", "kend": namelist.npz - 1},
            "condensation": {"serialname": "ne_cond"},

        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = AdjustNegativeTracers(
            self.stencil_factory,
            self.config,
            self.config.ntimes,
            convert_mm_day=inputs.pop("convt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
