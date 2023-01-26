from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import __INLINED, FORWARD, PARALLEL, computation, interval

import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import StencilFactory
from pace.dsl.typing import Float, FloatField, FloatFieldIJ
from pace.util import X_DIM, Y_DIM, Z_DIM
from pace.util.grid import GridData

from .._config import PhysicsConfig


@gtscript.function
def moist_total_energy(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temp,
    delp,
    moist_q,
    rgrav,
    c_air,
    c1_vapor,
    c1_liquid,
    c1_ice,
):
    q_liq = qrain + qliquid
    q_solid = qice + qsnow + qgraupel
    q_cond = q_liq + q_solid
    con = 1.0 - (qvapor + q_cond)
    if moist_q is True:
        cvm = con + qvapor * c1_vapor + q_liq * c1_liquid + q_solid * c1_ice
    else:
        cvm = 1.0 + qvapor * c1_vapor + q_liq * c1_liquid + q_solid * c1_ice
    return rgrav * cvm * c_air * temp * delp


def calc_total_energy(
    total_energy: FloatField,
    temp: FloatField,
    delp: FloatField,
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    c_air: Float,
    grav: Float,
    c1_vapor: Float,
    c1_liquid: Float,
    c1_ice: Float,
):
    from __externals__ import hydrostatic

    with computation(PARALLEL), interval(...):
        if __INLINED(hydrostatic is True):
            total_energy = -c_air * temp * delp
        else:
            total_energy = (
                -moist_total_energy(
                    qvapor,
                    qliquid,
                    qrain,
                    qice,
                    qsnow,
                    qgraupel,
                    temp,
                    delp,
                    True,
                    1.0 / grav,
                    c_air,
                    c1_vapor,
                    c1_liquid,
                    c1_ice,
                )
                * grav
            )


def moist_total_energy_water(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    temp: FloatField,
    ua: FloatField,
    va: FloatField,
    wa: FloatField,
    delp: FloatField,
    gsize: FloatFieldIJ,
    column_energy_change: FloatFieldIJ,
    vapor: FloatField,
    water: FloatFieldIJ,
    rain: FloatFieldIJ,
    ice: FloatFieldIJ,
    snow: FloatFieldIJ,
    graupel: FloatFieldIJ,
    sen: Float,
    stress: Float,
    timestep: Float,
    tot_energy: FloatField,
    tot_water: FloatField,
    total_energy_bot: FloatFieldIJ,
    total_water_bot: FloatFieldIJ,
    c1_vapor: Float,
    c1_liquid: Float,
    c1_ice: Float,
    rgrav: Float,
    lv00: Float,
    li00: Float,
    c_air: Float,
):
    """
    mtetw in Fortran
    """
    from __externals__ import hydrostatic, moist_q

    with computation(PARALLEL), interval(...):
        q_liq = qliquid + qrain
        q_solid = qice + qsnow + qgraupel
        q_cond = q_liq + q_solid
        con = 1.0 - (qvapor + q_cond)
        if __INLINED(moist_q is True):
            cvm = con + qvapor * c1_vapor + q_liq * c1_liquid + q_solid * c1_ice
        else:
            cvm = 1.0 + qvapor * c1_vapor + q_liq * c1_liquid + q_solid * c1_ice
        tot_energy = (cvm * temp + lv00 * qvapor - li00 * q_solid) * c_air
        if __INLINED(hydrostatic is True):
            tot_energy = tot_energy + 0.5 * (ua ** 2 + va ** 2)
        else:
            tot_energy = tot_energy + 0.5 * (ua ** 2 + va ** 2 + wa ** 2)
        tot_energy = rgrav * tot_energy * delp * gsize ** 2.0
        tot_water = rgrav * (qvapor + q_cond) * delp * gsize ** 2.0

    with computation(FORWARD), interval(-1, None):
        total_energy_bot = (
            column_energy_change
            + (lv00 * c_air * vapor - li00 * c_air * (ice + snow + graupel))
            * timestep
            / 86400
            + sen * timestep
            + stress * timestep
        ) * gsize ** 2.0
        total_water_bot = (
            (vapor + water + rain + ice + snow + graupel)
            * timestep
            / 86400
            * gsize ** 2.0
        )


def calc_sedimentation_energy_loss(
    energy_loss: FloatFieldIJ, column_energy_change: FloatFieldIJ, gsize: FloatFieldIJ
):
    with computation(FORWARD), interval(-1, None):
        energy_loss = column_energy_change * gsize ** 2.0


class MicrophysicsState:
    def __init__(
        self,
        qvapor: pace.util.Quantity,
    ):
        self.qvapor = qvapor
        pass


class Microphysics:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: pace.util.QuantityFactory,
        grid_data: GridData,
        namelist: PhysicsConfig,
    ):
        self.namelist = namelist
        self._max_timestep = self.namelist.mp_time
        self._ntimes = self.namelist.ntimes

        # atmospheric physics constants
        if self.namelist.hydrostatic:
            self._c_air = constants.CP_AIR
            self._c_vap = constants.CP_VAP
        else:
            self._c_air = constants.CV_AIR
            self._c_vap = constants.CV_VAP

        self._d0_vap = self._c_vap - constants.C_LIQ

        # scaled constants to reduce 32 bit floating errors
        self._lv00 = (constants.HLV - self._d0_vap * constants.TICE) / self._c_air
        self._li00 = constants.LI00 / self._c_air
        self._li20 = self._lv00 + self._li00

        self._d1_vap = self._d0_vap / self._c_air
        self._d1_ice = constants.DC_ICE / self._c_air

        self._c1_vap = self._c_vap / self._c_air
        self._c1_liquid = constants.C_LIQ / self._c_air
        self._c1_ice = constants.C_ICE / self._c_air

        # allocate memory, compile stencils, etc.
        self._cond = 0.0

        def make_quantity(**kwargs):
            return quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown")

        self._tot_energy_change = quantity_factory.zeros(
            dims=[X_DIM, Y_DIM], units="unknown"
        )
        self._adj_vmr = quantity_factory.ones(
            dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"
        )
        self._rain = make_quantity()

        self._set_timestepping(self.namelist.dt_atmos)  # will change from dt_atmos
        # when we inline microphysics

        self._convert_mm_day = 86400.0 * constants.RGRAV / self._split_timestep

        pass

    def gfdl_cloud_microphys_init(self, dt_atmos: float):
        pass

    def _update_timestep_if_needed(self, timestep: float):
        if timestep != self._full_timestep:
            self._set_timestepping(timestep)

    def _set_timestepping(self, full_timestep: float):
        self._ntimes = int(
            max(self._ntimes, full_timestep / min(full_timestep, self._max_timestep))
        )
        self._split_timestep = full_timestep / self._ntimes
        self._full_timestep = full_timestep
        pass

    def __call__(self, state: MicrophysicsState, timestep: float):
        pass
