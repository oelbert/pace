from gt4py.cartesian.gtscript import (
    __INLINED,
    BACKWARD,
    FORWARD,
    PARALLEL,
    computation,
    exp,
    interval,
    log,
    log10,
)

import pace.physics.stencils.SHiELD_microphysics.physical_functions as physfun
import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ
from pace.physics.stencils.SHiELD_microphysics.terminal_fall import TerminalFall
from pace.util import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM

from ..._config import MicroPhysicsConfig


def moist_heat_capacity(
    qvapor, qliquid, qrain, qice, qsnow, qgraupel, c1_ice, c1_liq, c1_vap
):
    q_liq = qliquid + qrain
    q_solid = qice + qsnow + qgraupel
    return 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice


def init_zeros_heat_cap_latent_heat_precip(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    temperature: FloatField,
    icpk: FloatField,
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
    column_water: FloatFieldIJ,
    column_rain: FloatFieldIJ,
    column_ice: FloatFieldIJ,
    column_snow: FloatFieldIJ,
    column_graupel: FloatFieldIJ,
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

        preflux_water = 0.0
        preflux_rain = 0.0
        preflux_ice = 0.0
        preflux_snow = 0.0
        preflux_graupel = 0.0
        vterminal_water = 0.0
        vterminal_rain = 0.0
        vterminal_ice = 0.0
        vterminal_snow = 0.0
        vterminal_graupel = 0.0

    with computation(FORWARD), interval(-1, None):
        column_water = 0.0
        column_rain = 0.0
        column_ice = 0.0
        column_snow = 0.0
        column_graupel = 0.0


def calc_terminal_velocity_rsg(
    q: FloatField,
    density: FloatField,
    density_factor: FloatField,
    v_terminal: FloatField,
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
        if __INLINED(const_v):
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


def calc_terminal_velocity_ice(
    qice: FloatField,
    temperature: FloatField,
    density: FloatField,
    v_terminal: FloatField,
):
    """
    Fortran name is term_ice
    """

    from __externals__ import aa, bb, cc, constant_v, dd, ee, ifflag, v_fac, v_max

    with computation(PARALLEL), interval(...):

        if __INLINED(constant_v):
            v_terminal = v_fac
        else:
            if qice < constants.QFMIN:
                v_terminal = 0.0
            else:
                tc = temperature - constants.TICE0
                if ifflag == 1:
                    v_terminal = (
                        (3.0 + log10(qice * density)) * (tc * (aa * tc + bb) + cc)
                        + dd * tc
                        + ee
                    )
                    v_terminal = 0.01 * v_fac * exp(v_terminal * log(10.0))
                else:  # ifflag == 2:
                    v_terminal = v_fac * 3.29 * exp(0.16 * log(qice * density))
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
    z_edge,
    z_terminal,
    z_surface,
    timestep,
    v_terminal,
    r1,
    tau_mlt,
    icpk,
    li00,
    c1_vap,
    c1_liq,
    c1_ice,
    ks,
    ke,
    is_,
    ie,
    js,
    je,
    mode,
):
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
            for k in range(ke - 1, ks - 1, -1):
                if v_terminal[i, j, k] < 1.0e-10:
                    continue
                if q_melt[i, j, k] > constants.QCMIN:
                    for m in range(k + 1, ke + 1):
                        if z_terminal[i, j, k + 1] >= z_edge[i, j, m]:
                            break
                        if (z_terminal[i, j, k] < z_edge[i, j, m + 1]) and (
                            temperature[i, j, m] > constants.TICE0
                        ):
                            cvm[i, j, k] = moist_heat_capacity(
                                qvapor[i, j, k],
                                qliquid[i, j, k],
                                qrain[i, j, k],
                                qice[i, j, k],
                                qsnow[i, j, k],
                                qgraupel[i, j, k],
                                c1_ice,
                                c1_liq,
                                c1_vap,
                            )
                            cvm[i, j, m] = moist_heat_capacity(
                                qvapor[i, j, m],
                                qliquid[i, j, m],
                                qrain[i, j, m],
                                qice[i, j, m],
                                qsnow[i, j, m],
                                qgraupel[i, j, m],
                                c1_ice,
                                c1_liq,
                                c1_vap,
                            )
                            # cvm[i, j, k] = physfun.moist_heat_capacity(
                            #     qvapor[i, j, k],
                            #     qliquid[i, j, k],
                            #     qrain[i, j, k],
                            #     qice[i, j, k],
                            #     qsnow[i, j, k],
                            #     qgraupel[i, j, k],
                            # )
                            # cvm[i, j, m] = physfun.moist_heat_capacity(
                            #     qvapor[i, j, m],
                            #     qliquid[i, j, m],
                            #     qrain[i, j, m],
                            #     qice[i, j, m],
                            #     qsnow[i, j, m],
                            #     qgraupel[i, j, m],
                            # )
                            dtime = min(
                                timestep,
                                (z_edge[i, j, m] - z_edge[i, j, m + 1])
                                / v_terminal[i, j, k],
                            )
                            dtime = min(1.0, dtime / tau_mlt)
                            sink = min(
                                q_melt[i, j, k] * delp[i, j, k] / delp[i, j, m],
                                dtime
                                * (temperature[i, j, m] - constants.TICE0)
                                / icpk[i, j, m],
                            )
                            q_melt[i, j, k] -= sink * delp[i, j, m] / delp[i, j, k]
                            if z_terminal[i, j, k] < z_surface[i, j]:
                                r1[i, j] += sink * delp[i, j, m]
                            else:
                                qrain[i, j, m] += sink

                            # these may be redundant depending on how dace copies
                            if mode == "ice":
                                qice[i, j, k] = q_melt[i, j, k]
                            elif mode == "snow":
                                qsnow[i, j, k] = q_melt[i, j, k]
                            else:
                                qgraupel[i, j, k] = q_melt[i, j, k]

                            cvm_tmp = moist_heat_capacity(
                                qvapor[i, j, k],
                                qliquid[i, j, k],
                                qrain[i, j, k],
                                qice[i, j, k],
                                qsnow[i, j, k],
                                qgraupel[i, j, k],
                                c1_ice,
                                c1_liq,
                                c1_vap,
                            )
                            # physfun.moist_heat_capacity(
                            #     qvapor[i, j, k],
                            #     qliquid[i, j, k],
                            #     qrain[i, j, k],
                            #     qice[i, j, k],
                            #     qsnow[i, j, k],
                            #     qgraupel[i, j, k],
                            # )
                            temperature[i, j, k] = (
                                temperature[i, j, k] * cvm[i, j, k]
                                - li00 * sink * delp[i, j, m] / delp[i, j, k]
                            ) / cvm_tmp
                            cvm_tmp = moist_heat_capacity(
                                qvapor[i, j, m],
                                qliquid[i, j, m],
                                qrain[i, j, m],
                                qice[i, j, m],
                                qsnow[i, j, m],
                                qgraupel[i, j, m],
                                c1_ice,
                                c1_liq,
                                c1_vap,
                            )
                            # physfun.moist_heat_capacity(
                            #     qvapor[i, j, m],
                            #     qliquid[i, j, m],
                            #     qrain[i, j, m],
                            #     qice[i, j, m],
                            #     qsnow[i, j, m],
                            #     qgraupel[i, j, m],
                            # )
                            temperature[i, j, m] = (
                                temperature[i, j, m] * cvm[i, j, m]
                            ) / cvm_tmp
                        if q_melt[i, j, k] < constants.QCMIN:
                            break


def calc_edge_and_terminal_height(
    z_surface: FloatFieldIJ,
    z_edge: FloatField,
    z_terminal: FloatField,
    delz: FloatField,
    v_terminal: FloatField,
):
    """
    Calculate grid cell edge heights and terminal fall heights for sedimentation
    Fortran name is zezt
    """
    from __externals__ import timestep

    with computation(FORWARD), interval(-1, None):
        z_surface = 0.0
        z_edge = z_surface
    with computation(BACKWARD), interval(0, -1):
        z_edge = z_edge[0, 0, 1] - delz
    with computation(FORWARD):
        with interval(0, 1):
            z_terminal = z_edge
        with interval(1, -1):
            z_terminal = z_edge - (
                (0.5 * timestep) * (v_terminal[0, 0, -1] + v_terminal)
            )
        with interval(-1, None):
            z_terminal = z_surface - timestep * v_terminal[0, 0, -1]
    with computation(FORWARD):
        with interval(1, None):
            if z_terminal >= z_terminal[0, 0, -1]:
                z_terminal = z_terminal[0, 0, -1] - constants.DZ_MIN_FLIP


def adjust_fluxes(flux: FloatField):
    with computation(FORWARD):
        with interval(0, 1):
            flux = max(0.0, flux)
    with computation(BACKWARD):
        with interval(1, None):
            flux = max(0.0, flux - flux[0, 0, -1])


class Sedimentation:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: pace.util.QuantityFactory,
        config: MicroPhysicsConfig,
        timestep: float,
    ):
        self._idx: GridIndexing = stencil_factory.grid_indexing
        self._is_ = self._idx.isc
        self._ie = self._idx.iec
        self._js = self._idx.jsc
        self._je = self._idx.jec
        self._ks = 0
        self._ke = config.npz - 1
        self.c1_vap = config.c1_vap
        self.c1_liq = config.c1_liq
        self.c1_ice = config.c1_ice

        if config.do_hail:
            self.tvag = config.tvah
            self.tvbg = config.tvbh
            self.mug = config.muh
            self.bling = config.blinh
        else:
            self.tvag = config.tvag
            self.tvbg = config.tvbg
            self.mug = config.mug
            self.bling = config.bling

        assert config.ifflag in [
            1,
            2,
        ], f"Ice Formation Flag must be 1 or 2 not {config.ifflag}"

        self._timestep = timestep
        self.config = config
        self.li00 = config.li00

        # allocate internal storages
        def make_quantity():
            return quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="unknown")

        self._z_surface = quantity_factory.zeros([X_DIM, Y_DIM], units="unknown")
        self._z_edge = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM], units="unknown"
        )
        self._z_terminal = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM], units="unknown"
        )
        self._icpk = make_quantity()
        self._cvm = make_quantity()

        # compile stencils
        self._init_zeros_heat_cap_latent_heat_precip = (
            stencil_factory.from_origin_domain(
                func=init_zeros_heat_cap_latent_heat_precip,
                externals={
                    "c1_vap": config.c1_vap,
                    "c1_liq": config.c1_liq,
                    "c1_ice": config.c1_ice,
                    "lv00": config.lv00,
                    "li00": config.li00,
                    "li20": config.li20,
                    "d1_vap": config.d1_vap,
                    "d1_ice": config.d1_ice,
                    "t_wfr": config.t_wfr,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )
        )

        if config.do_psd_ice_fall is False:
            self._calc_terminal_ice_velocity = stencil_factory.from_origin_domain(
                func=calc_terminal_velocity_ice,
                externals={
                    "ifflag": config.ifflag,
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

        self._calc_edge_and_terminal_height = stencil_factory.from_origin_domain(
            func=calc_edge_and_terminal_height,
            externals={"timestep": self._timestep},
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(add=(0, 0, 1)),
        )

        self._terminal_fall = TerminalFall(
            stencil_factory, quantity_factory, config, timestep
        )

        self._adjust_fluxes = stencil_factory.from_origin_domain(
            func=adjust_fluxes,
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
        """
        Sedimentation of cloud ice, snow, graupel or hail, and rain
        Args:
            qvapor (inout):
            qliquid (inout):
            qrain (inout):
            qice (inout):
            qsnow (inout):
            qgraupel (inout):
            ua (inout):
            va (inout):
            wa (inout):
            temperature (inout):
            delp (in):
            delz (in):
            density (in):
            density_factor (in):
            preflux_water (inout):
            preflux_rain (inout):
            preflux_ice (inout):
            preflux_snow (inout):
            preflux_graupel (inout):
            vterminal_water (inout):
            vterminal_rain (inout):
            vterminal_ice (inout):
            vterminal_snow (inout):
            vterminal_graupel (inout):
            column_energy_change (inout):
            column_water (inout):
            column_rain (inout):
            column_ice (inout):
            column_snow (inout):
            column_graupel (inout):
        """

        self._init_zeros_heat_cap_latent_heat_precip(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temperature,
            self._icpk,
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
            column_water,
            column_rain,
            column_ice,
            column_snow,
            column_graupel,
        )

        # Terminal fall and melting of falling cloud ice into rain:
        if self.config.do_psd_ice_fall:
            if self.config.const_vi is False:
                self._calc_terminal_rsg_velocity(
                    qice,
                    density,
                    density_factor,
                    vterminal_ice,
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
                    vterminal_ice,
                    self.config.vi_fac,
                    self.config.tvai,
                    self.config.tvbi,
                    self.config.mui,
                    self.config.blini,
                    self.config.vi_max,
                )
        else:
            self._calc_terminal_ice_velocity(qice, temperature, density, vterminal_ice)

        self._calc_edge_and_terminal_height(
            self._z_surface,
            self._z_edge,
            self._z_terminal,
            delz,
            vterminal_ice,
        )

        if self.config.do_sedi_melt:
            sedi_melt(
                qvapor.data[:],
                qliquid.data[:],
                qrain.data[:],
                qice.data[:],
                qsnow.data[:],
                qgraupel.data[:],
                self._cvm.data[:],
                temperature.data[:],
                delp.data[:],
                self._z_edge.data[:],
                self._z_terminal.data[:],
                self._z_surface.data[:],
                self._timestep,
                vterminal_ice.data[:],
                column_rain.data[:],
                self.config.tau_imlt,
                self._icpk.data[:],
                self.li00,
                self.c1_vap,
                self.c1_liq,
                self.c1_ice,
                self._ks,
                self._ke,
                self._is_,
                self._ie,
                self._js,
                self._je,
                "ice",
            )

        self._terminal_fall(
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
            delp,
            delz,
            vterminal_ice,
            self._z_edge,
            self._z_terminal,
            preflux_ice,
            column_ice,
            column_energy_change,
            "ice",
        )

        self._adjust_fluxes(preflux_ice)

        # Terminal fall and melting of falling snow into rain:
        if self.config.const_vs is False:
            self._calc_terminal_rsg_velocity(
                qsnow,
                density,
                density_factor,
                vterminal_snow,
                self.config.vs_fac,
                self.config.tvas,
                self.config.tvbs,
                self.config.mus,
                self.config.blins,
                self.config.vs_max,
            )
        else:
            self._calc_terminal_rsg_velocity_const(
                qsnow,
                density,
                density_factor,
                vterminal_snow,
                self.config.vs_fac,
                self.config.tvas,
                self.config.tvbs,
                self.config.mus,
                self.config.blins,
                self.config.vs_max,
            )

        self._calc_edge_and_terminal_height(
            self._z_surface,
            self._z_edge,
            self._z_terminal,
            delz,
            vterminal_snow,
        )

        if self.config.do_sedi_melt:
            sedi_melt(
                qvapor.data[:],
                qliquid.data[:],
                qrain.data[:],
                qice.data[:],
                qsnow.data[:],
                qgraupel.data[:],
                self._cvm.data[:],
                temperature.data[:],
                delp.data[:],
                self._z_edge.data[:],
                self._z_terminal.data[:],
                self._z_surface.data[:],
                self._timestep,
                vterminal_snow.data[:],
                column_rain.data[:],
                self.config.tau_smlt,
                self._icpk.data[:],
                self.li00,
                self.c1_vap,
                self.c1_liq,
                self.c1_ice,
                self._ks,
                self._ke,
                self._is_,
                self._ie,
                self._js,
                self._je,
                "snow",
            )

        self._terminal_fall(
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
            delp,
            delz,
            vterminal_snow,
            self._z_edge,
            self._z_terminal,
            preflux_snow,
            column_snow,
            column_energy_change,
            "snow",
        )

        self._adjust_fluxes(preflux_snow)

        # Terminal fall and melting of falling graupel into rain
        if self.config.const_vg is False:
            self._calc_terminal_rsg_velocity(
                qgraupel,
                density,
                density_factor,
                vterminal_graupel,
                self.config.vg_fac,
                self.tvag,
                self.tvbg,
                self.mug,
                self.bling,
                self.config.vg_max,
            )
        else:
            self._calc_terminal_rsg_velocity_const(
                qgraupel,
                density,
                density_factor,
                vterminal_graupel,
                self.config.vg_fac,
                self.tvag,
                self.tvbg,
                self.mug,
                self.bling,
                self.config.vg_max,
            )

        self._calc_edge_and_terminal_height(
            self._z_surface,
            self._z_edge,
            self._z_terminal,
            delz,
            vterminal_graupel,
        )

        if self.config.do_sedi_melt:
            sedi_melt(
                qvapor.data[:],
                qliquid.data[:],
                qrain.data[:],
                qice.data[:],
                qsnow.data[:],
                qgraupel.data[:],
                self._cvm.data[:],
                temperature.data[:],
                delp.data[:],
                self._z_edge.data[:],
                self._z_terminal.data[:],
                self._z_surface.data[:],
                self._timestep,
                vterminal_graupel.data[:],
                column_rain.data[:],
                self.config.tau_gmlt,
                self._icpk.data[:],
                self.li00,
                self.c1_vap,
                self.c1_liq,
                self.c1_ice,
                self._ks,
                self._ke,
                self._is_,
                self._ie,
                self._js,
                self._je,
                "graupel",
            )

        self._terminal_fall(
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
            delp,
            delz,
            vterminal_graupel,
            self._z_edge,
            self._z_terminal,
            preflux_graupel,
            column_graupel,
            column_energy_change,
            "graupel",
        )

        self._adjust_fluxes(preflux_graupel)

        # Terminal fall of cloud water
        if self.config.do_psd_water_fall:
            if self.config.const_vw is False:
                self._calc_terminal_rsg_velocity(
                    qliquid,
                    density,
                    density_factor,
                    vterminal_water,
                    self.config.vw_fac,
                    self.config.tvaw,
                    self.config.tvbw,
                    self.config.muw,
                    self.config.blinw,
                    self.config.vw_max,
                )
            else:
                self._calc_terminal_rsg_velocity_const(
                    qliquid,
                    density,
                    density_factor,
                    vterminal_water,
                    self.config.vw_fac,
                    self.config.tvaw,
                    self.config.tvbw,
                    self.config.muw,
                    self.config.blinw,
                    self.config.vw_max,
                )

            self._calc_edge_and_terminal_height(
                self._z_surface,
                self._z_edge,
                self._z_terminal,
                delz,
                vterminal_water,
            )

            self._terminal_fall(
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
                delp,
                delz,
                vterminal_water,
                self._z_edge,
                self._z_terminal,
                preflux_water,
                column_water,
                column_energy_change,
                "liquid",
            )

            self._adjust_fluxes(preflux_water)

        # terminal fall of rain
        if self.config.const_vr is False:
            self._calc_terminal_rsg_velocity(
                qrain,
                density,
                density_factor,
                vterminal_rain,
                self.config.vr_fac,
                self.config.tvar,
                self.config.tvbr,
                self.config.mur,
                self.config.blinr,
                self.config.vr_max,
            )
        else:
            self._calc_terminal_rsg_velocity_const(
                qrain,
                density,
                density_factor,
                vterminal_rain,
                self.config.vr_fac,
                self.config.tvar,
                self.config.tvbr,
                self.config.mur,
                self.config.blinr,
                self.config.vr_max,
            )

        self._calc_edge_and_terminal_height(
            self._z_surface,
            self._z_edge,
            self._z_terminal,
            delz,
            vterminal_rain,
        )

        self._terminal_fall(
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
            delp,
            delz,
            vterminal_rain,
            self._z_edge,
            self._z_terminal,
            preflux_rain,
            column_rain,
            column_energy_change,
            "rain",
        )

        self._adjust_fluxes(preflux_rain)
