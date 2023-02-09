import physical_functions as physfun
from gt4py.cartesian import gtscript  # noqa
from gt4py.cartesian.gtscript import (
    __INLINED,
    BACKWARD,
    FORWARD,
    PARALLEL,
    computation,
    exp,
    interval,
    log,
)

import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ
from pace.util import X_DIM, Y_DIM, Z_DIM

from ..._config import MicroPhysicsConfig


def init_heat_cap_latent_heat(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    temperature: FloatField,
    q_liq: FloatField,
    q_solid: FloatField,
    cvm: FloatField,
    te: FloatField,
    lcpk: FloatField,
    icpk: FloatField,
    tcpk: FloatField,
    tcp3: FloatField,
):

    with computation(PARALLEL), interval(...):
        (
            q_liq,
            q_solid,
            cvm,
            te,
            lcpk,
            icpk,
            tcpk,
            tcp3,
        ) = physfun.calc_heat_cap_and_latent_heat_coeff(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temperature,
        )


def calc_terminal_velocity_rsg(
    q: FloatField,
    density: FloatField,
    density_factor: FloatField,
    v_fac: float,
    tva: float,
    tvb: float,
    mu: float,
    blin: float,
    v_max: float,
):
    """
    Calculate terminal velocity for rain, snow, and graupel, Lin et al. (1983)
    Fortran name is term_rsg
    """

    from __externals__ import const_v

    with computation(PARALLEL), interval(...):
        if __INLINED(const_v is True):
            v_terminal = v_fac
        else:
            if q < constants.QFMIN:
                v_terminal = 0.0
            else:
                v_terminal = physfun.calc_terminal_velocity(
                    q, density, tva, tvb, mu, blin
                )
                v_terminal = v_fac * v_terminal * density_factor
                v_terminal = min(v_max, max(0.0, v_terminal))

        return v_terminal


def calc_terminal_velocity_ice(
    qice: FloatField,
    temperature: FloatField,
    density: FloatField,
    v_terminal: FloatField,
):
    """
    Fortran name is term_ice
    """

    from __externals__ import aa, bb, cc, constant_v, dd, ee, ifflag, tice, v_fac, v_max

    with computation(PARALLEL), interval(...):

        if __INLINED(constant_v is True):
            v_terminal = v_fac
        else:
            if qice < constants.QFMIN:
                v_terminal = 0.0
            else:
                tc = temperature - tice
                if ifflag == 1:
                    v_terminal = (
                        (3.0 + (log(qice * density) / log(10)))
                        * (tc * (aa * tc + bb) + cc)
                        + dd * tc
                        + ee
                    )
                    v_terminal + 0.01 * v_fac * exp(v_terminal * log(10.0))
                elif ifflag == 2:
                    v_terminal = v_fac * 3.29 * exp(0.16 * log(qice * density))
                else:
                    raise ValueError(f"Ice Formation Flag must be 1 or 2 not {ifflag}")
                v_terminal = min(v_max, max(0.0, v_terminal))


# TODO: Can this be made a stencil?
def sedi_melt(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    cvm,
    temperature,
    delp,
    ze,
    zt,
    zs,
    timestep,
    v_terminal,
    r1,
    tau_mlt,
    icpk,
    ks,
    ke,
    is_,
    ie,
    js,
    je,
    mode,
):
    li00 = -5
    if mode == "ice":
        q_melt = qice
    elif mode == "snow":
        q_melt = qsnow
    elif mode == "graupel":
        q_melt = qgraupel
    else:
        raise ValueError(f"sedi_melt mode {mode} not ice, snow, or graupel")
    for i in range(is_, ie + 1):
        for j in range(js, je + 1):
            for k in range(ke, ks - 1, -1):
                if v_terminal[i, j, k] >= 1.0e-10:
                    continue
                if q_melt[i, j, k] > constants.QCMIN:
                    for m in range(k + 1, ke + 1):
                        if zt[i, j, k + 1] >= ze[i, j, m]:
                            break
                        if (zt[i, j, k] < ze[i, j, m + 1]) and (
                            temperature[i, j, m] > constants.TICE
                        ):
                            cvm[i, j, k] = physfun.moist_heat_capacity(
                                qvapor[i, j, k],
                                qliquid[i, j, k],
                                qrain[i, j, k],
                                qice[i, j, k],
                                qsnow[i, j, k],
                                qgraupel[i, j, k],
                            )
                            cvm[i, j, m] = physfun.moist_heat_capacity(
                                qvapor[i, j, m],
                                qliquid[i, j, m],
                                qrain[i, j, m],
                                qice[i, j, m],
                                qsnow[i, j, m],
                                qgraupel[i, j, m],
                            )
                            dtime = min(
                                timestep,
                                (ze[i, j, m] - ze[i, j, m + 1]) / v_terminal[i, j, k],
                            )
                            dtime = min(1.0, dtime / tau_mlt)
                            sink = min(
                                q_melt[i, j, k] * delp[i, j, k] / delp[i, j, m],
                                dtime
                                * (temperature[i, j, m] - constants.TICE)
                                / icpk[i, j, m],
                            )
                            q_melt[i, j, k] -= sink * delp[i, j, m] / delp[i, j, k]
                            if zt[i, j, k] < zs[i, j]:
                                r1[i, j] += sink * delp[i, j, m]
                            else:
                                qrain[i, j, m] += sink
                            if mode == "ice":
                                qice[i, j, k] = q_melt[i, j, k]
                            elif mode == "snow":
                                qsnow[i, j, k] = q_melt[i, j, k]
                            else:
                                qgraupel[i, j, k] = q_melt[i, j, k]
                            temperature[i, j, k] = (
                                temperature[i, j, k] * cvm[i, j, k]
                                - li00 * sink * delp[i, j, m] / delp[i, j, k]
                            ) / physfun.moist_heat_capacity(
                                qvapor[i, j, k],
                                qliquid[i, j, k],
                                qrain[i, j, k],
                                qice[i, j, k],
                                qsnow[i, j, k],
                                qgraupel[i, j, k],
                            )
                            temperature[i, j, m] = (
                                temperature[i, j, m] * cvm[i, j, m]
                            ) / physfun.moist_heat_capacity(
                                qvapor[i, j, m],
                                qliquid[i, j, m],
                                qrain[i, j, m],
                                qice[i, j, m],
                                qsnow[i, j, m],
                                qgraupel[i, j, m],
                            )
                        if q_melt[i, j, k] < constants.QCMIN:
                            break

    return q_melt, qrain, r1, temperature, cvm


def ze_zt(
    zs: FloatFieldIJ,
    ze: FloatField,
    zt: FloatField,
    delz: FloatField,
    v_terminal: FloatField,
):
    """
    Calculate ze zt for sedimentation
    Forttan name is zezt
    """
    from __externals__ import timestep

    with computation(FORWARD), interval(-1, None):
        zs = 0.0
        ze = zs
    with computation(BACKWARD), interval(0, -1):
        ze = ze[0, 0, -1] - delz
    with computation(FORWARD):
        with interval(0, 1):
            zt = ze
        with interval(1, -1):
            zt = ze - (0.5 * timestep * (v_terminal[0, 0, -1] - v_terminal))
        with interval(-1, None):
            zt = zs - timestep * v_terminal[0, 0, -1]
        with interval(1, None):
            if zt > zt[0, 0, -1]:
                zt = zt[0, 0, -1] - constants.DZ_MIN_FLIP


class Sedimentation:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: pace.util.QuantityFactory,
        config: MicroPhysicsConfig,
        timestep: float,
        convert_mm_day: float,
    ):
        self._idx: GridIndexing = stencil_factory.grid_indexing

        self._timestep = timestep
        self._convert_mm_day = convert_mm_day
        self.config = config

        # allocate internal storages
        def make_quantity():
            return quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="unknown")

        self._zs = make_quantity()
        self._ze = make_quantity()
        self._zt = make_quantity()

        # compile stencils
        self._init_heat_cap_latent_heat = stencil_factory.from_origin_domain(
            func=init_heat_cap_latent_heat,
            externals={
                "c1_vapor": config.c1_vap,
                "c1_liq": config.c1_liq,
                "c1_ice": config.c1_ice,
                "lv00": config.lv00,
                "li00": config.li00,
                "li20": config.li20,
                "d1_vap": config.d1_vap,
                "d1_ice": config.d1_ice,
                "tice": config.tice,
                "t_wfr": config.t_wfr,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        if config.do_psd_ice_fall is False:
            self._calc_terminal_ice_velocity = stencil_factory.from_origin_domain(
                func=calc_terminal_velocity_ice,
                externals={
                    "ifflag": config.ifflag,
                    "tice": config.tice,
                    "constant_v": config.const_vi,
                    "v_fac": config.vi_fac,
                    "v_max": config.vi_max,
                    "aa": -4.14122e-5,
                    "bb": -0.00538922,
                    "cc": -0.0516344,
                    "dd": 0.00216078,
                    "ee": 1.9714,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        if False in [
            config.const_vr,
            config.const_vi,
            config.const_vs,
            config.const_vg,
        ]:
            self._calc_terminal_rsg_velocity = stencil_factory.from_origin_domain(
                func=calc_terminal_velocity_rsg,
                externals={"const_v": False},
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        if True in [config.const_vr, config.const_vi, config.const_vs, config.const_vg]:
            self._calc_terminal_rsg_velocity_const = stencil_factory.from_origin_domain(
                func=calc_terminal_velocity_rsg,
                externals={"const_v": True},
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        self._zezt = stencil_factory.from_origin_domain(
            func=ze_zt,
            externals={"timestep": self._timestep},
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

    def __call__(
        self,
        qvapor: FloatField,
        qliquid: FloatField,
        qrain: FloatField,
        qice: FloatField,
        qsnow: FloatField,
        qgraupel: FloatField,
        ua: FloatField,
        va: FloatField,
        wa: FloatField,
        temperature: FloatField,
        delp: FloatField,
        delz: FloatField,
        density: FloatField,
        density_factor: FloatField,
        preflux_water: FloatField,
        preflux_rain: FloatField,
        preflux_ice: FloatField,
        preflux_snow: FloatField,
        preflux_graupel: FloatField,
        vterminal_water: FloatField,
        vterminal_rain: FloatField,
        vterminal_ice: FloatField,
        vterminal_snow: FloatField,
        vterminal_graupel: FloatField,
        column_energy_change: FloatFieldIJ,
        column_water: FloatFieldIJ,
        column_rain: FloatFieldIJ,
        column_ice: FloatFieldIJ,
        column_snow: FloatFieldIJ,
        column_graupel: FloatFieldIJ,
    ):

        if self.config.do_psd_ice_fall:
            self._calc_terminal_ice_velocity(qice, temperature, density, vterminal_ice)
        else:
            if self.config.const_vi is False:
                self._calc_terminal_rsg_velocity(
                    qice,
                    density,
                    density_factor,
                    self.config.vi_fac,
                    self.config.tvai,
                    self.config.tvbi,
                    self.config.mui,
                    self.config.blini,
                    self.config.vi_max,
                )
            else:
                self._calc_terminal_rsg_velocity_const(
                    qice,
                    density,
                    density_factor,
                    self.config.vi_fac,
                    self.config.tvai,
                    self.config.tvbi,
                    self.config.mui,
                    self.config.blini,
                    self.config.vi_max,
                )

        self._zezt(
            self._zs,
            self._ze,
            self._zt,
            delz,
            vterminal_ice,
        )

        return (
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            ua,
            va,
            wa,
            temperature,
            preflux_water,
            preflux_rain,
            preflux_ice,
            preflux_snow,
            preflux_graupel,
            vterminal_water,
            vterminal_rain,
            vterminal_ice,
            vterminal_snow,
            vterminal_graupel,
            column_energy_change,
            column_water,
            column_rain,
            column_ice,
            column_snow,
            column_graupel,
        )
