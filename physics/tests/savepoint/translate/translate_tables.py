from gt4py.cartesian.gtscript import FORWARD, computation, interval

import pace.dsl
import pace.physics.stencils.microphysics_v3.physical_functions as physfun
import pace.util
from pace.dsl.typing import FloatField
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.humidity_tables import (
    HumiditySaturationTables,
)
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


def init_tables(
    temp: FloatField,
    table0: FloatField,
    table2: FloatField,
):
    with computation(FORWARD), interval(...):
        table0 = physfun.table0(temp)
        table2 = physfun.table2(temp)


def calc_table_values(
    temp: FloatField,
    den: FloatField,
    iqs: FloatField,
    wqs: FloatField,
    didt: FloatField,
    dwdt: FloatField,
):
    with computation(FORWARD), interval(...):
        wqs, dwdt = physfun.sat_spec_hum_water(temp, den)
        iqs, didt = physfun.sat_spec_hum_water_ice(temp, den)


class CalcTables:
    def __init__(self, stencil_factory: pace.dsl.StencilFactory, config):
        self._idx = stencil_factory.grid_indexing
        self.config = config
        self._init_tables = stencil_factory.from_origin_domain(
            func=init_tables,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )
        self._calc_table_values = stencil_factory.from_origin_domain(
            func=calc_table_values,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

    def __call__(
        self,
        temp: FloatField,
        table0: FloatField,
        table2: FloatField,
        iqs: FloatField,
        wqs: FloatField,
        didt: FloatField,
        dwdt: FloatField,
        temp2: FloatField,
        den: FloatField,
    ):
        self._init_tables(
            temp,
            table0,
            table2,
        )
        self._calc_table_values(
            temp2,
            den,
            iqs,
            wqs,
            didt,
            dwdt,
        )


class LookupPython:
    def __init__(self, length):
        self.sat_tables = HumiditySaturationTables(length)

    def __call__(
        self,
        index: FloatField,
        table0: FloatField,
        table2: FloatField,
        iqs: FloatField,
        wqs: FloatField,
        didt: FloatField,
        dwdt: FloatField,
        temp: FloatField,
        den: FloatField,
    ):
        table0.view[:] = self.sat_tables.table0[index.view[:] - 1]
        table2.view[:] = self.sat_tables.table2[index.view[:] - 1]
        wqs.view[:], dwdt.view[:] = self.sat_tables.sat_water(temp.view[:], den.view[:])
        iqs.view[:], didt.view[:] = self.sat_tables.sat_ice_water(
            temp.view[:], den.view[:]
        )


class TranslatePythonTables(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)

        self.in_vars["data_vars"] = {
            "index": {"serialname": "tc_index", "mp3": True},
            "table0": {"serialname": "tc_t0", "mp3": True},
            "table2": {"serialname": "tc_t2", "mp3": True},
            "wqs": {"serialname": "tab_wq", "mp3": True},
            "dwdt": {"serialname": "tab_dwq", "mp3": True},
            "iqs": {"serialname": "tab_iq", "mp3": True},
            "didt": {"serialname": "tab_diq", "mp3": True},
            "temp": {"serialname": "tab_pt", "mp3": True},
            "den": {"serialname": "tab_den", "mp3": True},
        }

        self.out_vars = {
            "table0": {"serialname": "tc_t0", "kend": namelist.npz, "mp3": True},
            "table2": {"serialname": "tc_t2", "kend": namelist.npz, "mp3": True},
            "wqs": {"serialname": "tab_wq", "kend": namelist.npz, "mp3": True},
            "dwdt": {"serialname": "tab_dwq", "kend": namelist.npz, "mp3": True},
            "iqs": {"serialname": "tab_iq", "kend": namelist.npz, "mp3": True},
            "didt": {"serialname": "tab_diq", "kend": namelist.npz, "mp3": True},
        }

        self.stencil_factory = stencil_factory
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = LookupPython(2621)

        compute_func(**inputs)

        return self.slice_output(inputs)


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
            "wqs": {"serialname": "tab_wq", "mp3": True},
            "dwdt": {"serialname": "tab_dwq", "mp3": True},
            "iqs": {"serialname": "tab_iq", "mp3": True},
            "didt": {"serialname": "tab_diq", "mp3": True},
            "temp2": {"serialname": "tab_pt", "mp3": True},
            "den": {"serialname": "tab_den", "mp3": True},
        }

        self.out_vars = {
            "table0": {"serialname": "tc_t0", "kend": namelist.npz, "mp3": True},
            "table2": {"serialname": "tc_t2", "kend": namelist.npz, "mp3": True},
            "wqs": {"serialname": "tab_wq", "kend": namelist.npz, "mp3": True},
            "dwdt": {"serialname": "tab_dwq", "kend": namelist.npz, "mp3": True},
            "iqs": {"serialname": "tab_iq", "kend": namelist.npz, "mp3": True},
            "didt": {"serialname": "tab_diq", "kend": namelist.npz, "mp3": True},
        }

        self.max_error = 1.5e-14  # 10^-25 absolute errors at the top of the tables

        self.stencil_factory = stencil_factory
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = CalcTables(
            self.stencil_factory,
            self.config,
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
