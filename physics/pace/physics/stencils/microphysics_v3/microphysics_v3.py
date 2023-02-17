from gt4py.cartesian.gtscript import (
    __INLINED,
    FORWARD,
    PARALLEL,
    computation,
    interval,
    sqrt,
)
from moist_total_energy import calc_total_energy

import pace.fv3core.stencils.basic_operations as basic
import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import Float, FloatField, FloatFieldIJ
from pace.physics.stencils.microphysics_v3.mp_fast import FastMicrophysics
from pace.physics.stencils.microphysics_v3.mp_full import FullMicrophysics
from pace.physics.stencils.microphysics_v3.neg_adj import AdjustNegativeTracers
from pace.util import X_DIM, Y_DIM, Z_DIM, Timer
from pace.util.grid import GridData

from ..._config import MicroPhysicsConfig


def calc_sedimentation_energy_loss(
    energy_loss: FloatFieldIJ, column_energy_change: FloatFieldIJ, gsize: FloatFieldIJ
):
    """
    Last calculation of mtetw, split into a separate stencil to cover the
    if (present(te_loss)) conditional
    """
    with computation(FORWARD), interval(-1, None):
        energy_loss = column_energy_change * gsize ** 2.0


def moist_total_energy_and_water(
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
    tot_energy: FloatField,
    tot_water: FloatField,
    total_energy_bot: FloatFieldIJ,
    total_water_bot: FloatFieldIJ,
):
    """
    Fortran name is mtetw
    """
    from __externals__ import (
        c1_ice,
        c1_liq,
        c1_vap,
        c_air,
        hydrostatic,
        li00,
        lv00,
        moist_q,
        timestep,
    )

    with computation(PARALLEL), interval(...):
        q_liq = qliquid + qrain
        q_solid = qice + qsnow + qgraupel
        q_cond = q_liq + q_solid
        if __INLINED(moist_q is True):
            con = 1.0 - (qvapor + q_cond)
            cvm = con + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
        else:
            cvm = 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
        tot_energy = (cvm * temp + lv00 * qvapor - li00 * q_solid) * c_air
        if __INLINED(hydrostatic is True):
            tot_energy = tot_energy + 0.5 * (ua ** 2 + va ** 2)
        else:
            tot_energy = tot_energy + 0.5 * (ua ** 2 + va ** 2 + wa ** 2)
        tot_energy = constants.RGRAV * tot_energy * delp * gsize ** 2.0
        tot_water = constants.RGRAV * (qvapor + q_cond) * delp * gsize ** 2.0

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


def convert_specific_to_mass_mixing_ratios_and_calculate_densities(
    qzvapor: FloatField,
    qzliquid: FloatField,
    qzrain: FloatField,
    qzice: FloatField,
    qzsnow: FloatField,
    qzgraupel: FloatField,
    qza: FloatField,
    delp_specific: FloatField,
    density: FloatField,
    pz: FloatField,
    bottom_density: FloatFieldIJ,
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    qa: FloatField,
    delp: FloatField,
    delz: FloatField,
    temp: FloatField,
):
    from __externals__ import inline_mp

    with computation(PARALLEL):
        with interval(...):
            if __INLINED(inline_mp is True):
                con_r8 = 1.0 - (qvapor + qliquid + qrain + qice + qsnow + qgraupel)
            else:
                con_r8 = 1.0 - qvapor
            delp_specific = delp * con_r8
            rcon_r8 = 1.0 / con_r8
            qzvapor = qvapor * rcon_r8
            qzliquid = qliquid * rcon_r8
            qzrain = qrain * rcon_r8
            qzice = qice * rcon_r8
            qzsnow = qsnow * rcon_r8
            qzgraupel = qgraupel * rcon_r8
            qza = qa

            # Dry air density and layer-mean pressure thickness
            density = -delp_specific * (constants.GRAV * delz)
            pz = density * constants.RDGAS * temp
        with interval(-1, None):
            bottom_density = density


def calc_density_factor(
    den_fac: FloatField,
    density: FloatField,
    bottom_density: FloatFieldIJ,
):
    with computation(PARALLEL), interval(...):
        den_fac = sqrt(bottom_density / density)


def cloud_nuclei(
    geopotential_height: FloatFieldIJ,
    qnl: FloatField,
    qni: FloatField,
    density: FloatField,
    cloud_condensation_nuclei: FloatField,
    cloud_ice_nuclei: FloatField,
):
    from __externals__ import ccn_l, ccn_o, prog_ccn

    with computation(PARALLEL), interval(...):
        if __INLINED(prog_ccn is True):
            # boucher and lohmann (1995)
            nl = min(1.0, abs(geopotential_height / (10.0 * constants.GRAV))) * (
                10.0 ** 2.24 * (qnl * density * 1.0e9) ** 0.257
            ) + (1.0 - min(1.0, abs(geopotential_height) / (10 * constants.GRAV))) * (
                10.0 ** 2.06 * (qnl * density * 1.0e9) ** 0.48
            )
            ni = qni
            cloud_condensation_nuclei = (max(10.0, nl) * 1.0e6) / density
            cloud_ice_nuclei = (max(10.0, ni) * 1.0e6) / density
        else:
            cloud_condensation_nuclei = (
                (
                    ccn_l * min(1.0, abs(geopotential_height) / (10.0 * constants.GRAV))
                    + ccn_o
                    * (
                        1.0
                        - min(1.0, abs(geopotential_height) / (10.0 * constants.GRAV))
                    )
                )
                * 1.0e6
                / density
            )
            cloud_ice_nuclei = 0.0 / density


def subgrid_deviation_and_relative_humidity(
    gsize: FloatFieldIJ,
    geopotential_height: FloatFieldIJ,
    h_var: FloatFieldIJ,
    rh_adj: FloatFieldIJ,
    rh_rain: FloatFieldIJ,
):
    from __externals__ import dw_land, dw_ocean, rh_inc, rh_inr

    with computation(PARALLEL), interval(-1, None):
        t_lnd = dw_land * sqrt(gsize / 1.0e5)
        t_ocn = dw_ocean * sqrt(gsize / 1.0e5)
        tmp = min(1.0, abs(geopotential_height) / (10.0 * constants.GRAV))
        hvar = t_lnd * tmp + t_ocn * (1.0 - tmp)
        h_var = min(0.20, max(0.01, hvar))

        rh_adj = 1.0 - h_var - rh_inc
        rh_rain = max(0.35, rh_adj - rh_inr)


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
        config: MicroPhysicsConfig,
        do_mp_fast: bool = False,
        do_mp_full: bool = True,
    ):
        self.config = config
        self._idx: GridIndexing = stencil_factory.grid_indexing
        self._max_timestep = self.config.mp_time
        self._ntimes = self.config.ntimes
        self.do_mp_fast = do_mp_fast
        self.do_mp_full = do_mp_fast

        # allocate memory, compile stencils, etc.
        def make_quantity(**kwargs):
            return quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown")

        self._tot_energy_change = quantity_factory.zeros(
            dims=[X_DIM, Y_DIM], units="unknown"
        )
        self._adj_vmr = quantity_factory.ones(
            dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"
        )
        self._rain = make_quantity()

        self._rh_adj = quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="unknown")
        self._rh_rain = quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="unknown")

        self._set_timestepping(self.config.dt_atmos)  # will change from dt_atmos
        # when we inline microphysics

        self._convert_mm_day = 86400.0 * constants.RGRAV / self.split_timestep

        self._copy_stencil = stencil_factory.from_origin_domain(
            basic.copy_defn,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._calc_total_energy = stencil_factory.from_origin_domain(
            func=calc_total_energy,
            externals={
                "hydrostatic": self.config.hydrostatic,
                "c_air": self.config.c_air,
                "c1_vap": self.config.c1_vap,
                "c1_liq": self.config.c1_liq,
                "c1_ice": self.config.c1_ice,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        if self.config.consv_checker is True:
            self._moist_total_energy_and_water_mq = stencil_factory.from_origin_domain(
                func=moist_total_energy_and_water,
                externals={
                    "hydrostatic": self.config.hydrostatic,
                    "moist_q": True,
                    "c1_vap": self.config.c1_vap,
                    "c1_liq": self.config.c1_liq,
                    "c1_ice": self.config.c1_ice,
                    "lv00": self.config.lv00,
                    "li00": self.config.li00,
                    "c_air": self.config.c_air,
                    "timestep": self.full_timestep,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

            self._moist_total_energy_and_water = stencil_factory.from_origin_domain(
                func=moist_total_energy_and_water,
                externals={
                    "hydrostatic": self.config.hydrostatic,
                    "moist_q": False,
                    "c1_vap": self.config.c1_vap,
                    "c1_liq": self.config.c1_liq,
                    "c1_ice": self.config.c1_ice,
                    "lv00": self.config.lv00,
                    "li00": self.config.li00,
                    "c_air": self.config.c_air,
                    "timestep": self.full_timestep,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

            self._calc_sedimentation_energy_loss = stencil_factory.from_origin_domain(
                func=calc_sedimentation_energy_loss,
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        self._convert_specific_to_mass_mixing_ratios_and_calculate_densities = (
            stencil_factory.from_origin_domain(
                func=convert_specific_to_mass_mixing_ratios_and_calculate_densities,
                externals={
                    "inline_mp": self.config.do_inline_mp,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )
        )

        self._calc_density_factor = stencil_factory.from_origin_domain(
            func=calc_density_factor,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._cloud_nuclei = stencil_factory.from_origin_domain(
            func=cloud_nuclei,
            externals={
                "prog_ccn": self.config.prog_ccn,
                "ccn_l": self.config.ccn_l,
                "ccn_o": self.config.ccn_o,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._subgrid_deviation_and_relative_humidity = (
            stencil_factory.from_origin_domain(
                func=subgrid_deviation_and_relative_humidity,
                externals={
                    "dw_land": self.config.dw_land,
                    "dw_ocean": self.config.dw_ocean,
                    "rh_inc": self.config.rh_inc,
                    "rh_inr": self.config.rh_inr,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )
        )

        self._neg_adj = AdjustNegativeTracers(
            stencil_factory,
            self.config.adjustnegative,
            self._ntimes,
            self._convert_mm_day,
        )

        if self.do_mp_fast:
            self._mp_fast = FastMicrophysics(
                stencil_factory,
                self.config.fastmp,
                self.full_timestep,
                self._convert_mm_day,
            )

        if self.do_mp_full:
            self._mp_full = FullMicrophysics(
                stencil_factory,
                quantity_factory,
                config,
                self.split_timestep,
                self._ntimes,
                self._convert_mm_day,
            )

        pass

    def _update_timestep_if_needed(self, timestep: float):
        if timestep != self.full_timestep:
            self._set_timestepping(timestep)

    def _set_timestepping(self, full_timestep: float):
        """
        sets full and split timesteps
        full_timestep is equivalent to dtm
        split_timestep is equivalent to dts
        """
        self._ntimes = int(
            max(self._ntimes, full_timestep / min(full_timestep, self._max_timestep))
        )
        self.split_timestep = full_timestep / self._ntimes
        self.full_timestep = full_timestep

    def __call__(
        self,
        state: MicrophysicsState,
        timer: Timer = pace.util.NullTimer(),
    ):
        pass
