from gt4py.cartesian.gtscript import FORWARD, computation, interval

import pace.dsl
import pace.physics.stencils.microphysics_v3.physical_functions as physfun
import pace.util
from pace.dsl.typing import FloatField
from pace.physics._config import PhysicsConfig
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


def init_tables(
    temp: FloatField,
    table0: FloatField,
    table2: FloatField,
):
    with computation(FORWARD), interval(...):
        table0 = physfun.table0(temp)
        table2 = physfun.table2(temp)


class InitTables:
    def __init__(self, stencil_factory: pace.dsl.StencilFactory, config):
        self._idx = stencil_factory.grid_indexing
        self.config = config
        self._init_tables = stencil_factory.from_origin_domain(
            func=init_tables,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

    def __call__(
        self,
        temp: FloatField,
        table0: FloatField,
        table2: FloatField,
    ):
        self._init_tables(
            temp,
            table0,
            table2,
        )


class TranslateTableComputation(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)

        self.in_vars["data_vars"] = {
            "temp": {"serialname": "tc_temp", "mp3": True},
            "table0": {"serialname": "tc_t0", "mp3": True},
            "table2": {"serialname": "tc_t2", "mp3": True},
        }

        self.out_vars = {
            "table0": {"serialname": "tc_t0", "kend": namelist.npz, "mp3": True},
            "table2": {"serialname": "tc_t2", "kend": namelist.npz, "mp3": True},
        }

        self.stencil_factory = stencil_factory
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = InitTables(
            self.stencil_factory,
            self.config,
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
