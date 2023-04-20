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
            "qvapor": {"serialname": "ne_qv", "mp3": True},
            "qliquid": {"serialname": "ne_ql", "mp3": True},
            "qrain": {"serialname": "ne_qr", "mp3": True},
            "qice": {"serialname": "ne_qi", "mp3": True},
            "qsnow": {"serialname": "ne_qs", "mp3": True},
            "qgraupel": {"serialname": "ne_qg", "mp3": True},
            "temperature": {"serialname": "ne_pt", "mp3": True},
            "delp": {"serialname": "ne_delp", "mp3": True},
            "condensation": {"serialname": "ne_cond", "mp3": True},
        }

        self.in_vars["parameters"] = ["convt"]

        self.out_vars = {
            "qvapor": {"serialname": "ne_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "ne_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "ne_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "ne_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "ne_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "ne_qg", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "ne_pt", "kend": namelist.npz, "mp3": True},
            "delp": {"serialname": "ne_delp", "kend": namelist.npz, "mp3": True},
            "condensation": {"serialname": "ne_cond", "mp3": True},
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
