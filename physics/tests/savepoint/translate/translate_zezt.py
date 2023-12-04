import pace.dsl
import pace.util
from pace.dsl.typing import FloatField, FloatFieldIJ
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.SHiELD_microphysics.sedimentation import (
    calc_edge_and_terminal_height,
)
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


# from pace.util import X_DIM, Y_DIM, Z_INTERFACE_DIM


class ZeZt:
    def __init__(
        self,
        stencil_factory: pace.dsl.StencilFactory,
        config,
        timestep,
    ):
        self._idx = stencil_factory.grid_indexing
        self.config = config
        self._calc_edge_and_terminal_height = stencil_factory.from_origin_domain(
            func=calc_edge_and_terminal_height,
            externals={"timestep": timestep},
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(add=(0, 0, 1)),
        )

    def __call__(
        self,
        z_surface: FloatFieldIJ,
        z_edge: FloatField,
        z_terminal: FloatField,
        delz: FloatField,
        v_terminal: FloatField,
    ):
        self._calc_edge_and_terminal_height(
            z_surface,
            z_edge,
            z_terminal,
            delz,
            v_terminal,
        )


class TranslateZeZt(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "z_surface": {"serialname": "zz_zs", "mp3": True},
            "z_edge": {"serialname": "zz_ze", "mp3": True},
            "z_terminal": {"serialname": "zz_zt", "mp3": True},
            "delz": {"serialname": "zz_dz", "mp3": True},
            "v_terminal": {"serialname": "zz_vt", "mp3": True},
        }

        self.in_vars["parameters"] = ["dt"]

        self.out_vars = {
            "z_edge": {"serialname": "zz_ze", "kend": namelist.npz + 1, "mp3": True},
            "z_terminal": {
                "serialname": "zz_zt",
                "kend": namelist.npz + 1,
                "mp3": True,
            },
            "z_surface": {"serialname": "zz_zs", "mp3": True},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = ZeZt(
            self.stencil_factory,
            self.config,
            timestep=inputs.pop("dt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
