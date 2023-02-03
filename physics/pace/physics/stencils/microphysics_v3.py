from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import (
    __INLINED,
    BACKWARD,
    FORWARD,
    PARALLEL,
    computation,
    exp,
    interval,
    sqrt,
)

import pace.fv3core.stencils.basic_operations as basic
import pace.physics.functions.table_functions as tables
import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import Float, FloatField, FloatFieldIJ
from pace.util import X_DIM, Y_DIM, Z_DIM
from pace.util.grid import GridData

from .._config import PhysicsConfig


QCMIN = 1.0e-15
DT_FR = 8.0


@gtscript.function
def positive_diff(x, y):
    """
    equivalent of Fortran dim
    """
    return max(x - y, 0.0)


@gtscript.function
def calc_moist_total_energy(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temp,
    delp,
    moist_q,
):
    from __externals__ import c1_ice, c1_liquid, c1_vapor, c_air

    q_liq = qrain + qliquid
    q_solid = qice + qsnow + qgraupel
    q_cond = q_liq + q_solid
    con = 1.0 - (qvapor + q_cond)
    if moist_q is True:
        cvm = con + qvapor * c1_vapor + q_liq * c1_liquid + q_solid * c1_ice
    else:
        cvm = 1.0 + qvapor * c1_vapor + q_liq * c1_liquid + q_solid * c1_ice
    return constants.RGRAV * cvm * c_air * temp * delp


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
):
    from __externals__ import c_air, hydrostatic

    with computation(PARALLEL), interval(...):
        if __INLINED(hydrostatic is True):
            total_energy = -c_air * temp * delp
        else:
            total_energy = (
                -calc_moist_total_energy(
                    qvapor,
                    qliquid,
                    qrain,
                    qice,
                    qsnow,
                    qgraupel,
                    temp,
                    delp,
                    True,
                )
                * constants.GRAV
            )


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
    mtetw in Fortran
    """
    from __externals__ import (
        c1_ice,
        c1_liquid,
        c1_vapor,
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
            cvm = con + qvapor * c1_vapor + q_liq * c1_liquid + q_solid * c1_ice
        else:
            cvm = 1.0 + qvapor * c1_vapor + q_liq * c1_liquid + q_solid * c1_ice
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


def calc_sedimentation_energy_loss(
    energy_loss: FloatFieldIJ, column_energy_change: FloatFieldIJ, gsize: FloatFieldIJ
):
    """
    Last calculation of mtetw, split into a separate stencil to cover the
    if (present(te_loss)) conditional
    """
    with computation(FORWARD), interval(-1, None):
        energy_loss = column_energy_change * gsize ** 2.0


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


@gtscript.function
def calc_heat_cap_and_latent_heat_coeff(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temp,
):
    """
    Fortran name is cal_mhc_lhc
    """

    from __externals__ import (
        c1_ice,
        c1_liq,
        c1_vapor,
        d1_ice,
        d1_vap,
        li00,
        li20,
        lv00,
        t_wfr,
        tice,
    )

    q_liq = qliquid + qrain
    q_solid = qice + qsnow + qgraupel
    cvm = 1.0 + qvapor * c1_vapor + q_liq * c1_liq + q_solid * c1_ice
    te = cvm + temp + lv00 * qvapor - li00 * q_solid
    lcpk = (lv00 + d1_vap * temp) / cvm
    icpk = (li00 + d1_ice * temp) / cvm
    tcpk = (li20 + (d1_vap + d1_ice) * temp) / cvm
    tcp3 = lcpk + icpk * min(1.0, positive_diff(tice, temp) / (tice - t_wfr))

    return q_liq, q_solid, cvm, te, lcpk, icpk, tcpk, tcp3


@gtscript.function
def update_hydrometeors_and_temperatures(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    delta_vapor,
    delta_liquid,
    delta_rain,
    delta_ice,
    delta_snow,
    delta_graupel,
    te,
):
    """
    Fortran name is update_qt
    """

    from __externals__ import (
        c1_ice,
        c1_liq,
        c1_vapor,
        d1_ice,
        d1_vap,
        li00,
        li20,
        lv00,
        t_wfr,
        tice,
    )

    qvapor += delta_vapor
    qliquid += delta_liquid
    qrain += delta_rain
    qice += delta_ice
    qsnow += delta_snow
    qgraupel += delta_graupel

    q_l = qrain + qliquid
    q_solid = qice + qsnow + qgraupel
    cvm = 1.0 + qvapor * c1_vapor + q_l * c1_liq + q_solid * c1_ice

    tk = (te - lv00 + qvapor + li00 * (qice + qsnow + qgraupel)) / cvm
    lcpk = (lv00 + d1_vap * tk) / cvm
    icpk = (li00 + d1_ice * tk) / cvm
    tcpk = (li20 + (d1_vap + d1_ice) * tk) / cvm
    tcp3 = lcpk + icpk * min(1.0, positive_diff(tice, tk) / (tice - t_wfr))

    return (
        qvapor,
        qliquid,
        qrain,
        qice,
        qsnow,
        qgraupel,
        cvm,
        tk,
        lcpk,
        icpk,
        tcpk,
        tcp3,
    )


def adjust_negative_tracers(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    temp: FloatField,
    delp: FloatField,
    condensation: FloatFieldIJ,
):
    """
    Fortran name is neg_adj
    """
    from __externals__ import convt, ntimes

    with computation(PARALLEL), interval(...):
        upper_fix = 0.0  # type: FloatField
        lower_fix = 0.0  # type: FloatField

    with computation(FORWARD), interval(...):
        (
            q_liq,
            q_solid,
            cvm,
            te,
            lcpk,
            icpk,
            tcpk,
            tcp3,
        ) = calc_heat_cap_and_latent_heat_coeff
        (qvapor, qliquid, qrain, qice, qsnow, qgraupel, temp)

        # if cloud ice < 0, borrow from snow
        if qice < 0.0:
            sink = min(-qice, max(0, qsnow))
            qice += sink
            qsnow -= sink

        # if snow < 0, borrow from graupel
        if qsnow < 0.0:
            sink = min(-qsnow, max(0, qgraupel))
            qsnow += sink
            qgraupel -= sink

        # if graupel < 0, borrow from rain
        if qgraupel < 0.0:
            sink = min(-qgraupel, max(0, qrain))
            (
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                cvm,
                temp,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            ) = update_hydrometeors_and_temperatures(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                0.0,
                0.0,
                -sink,
                0.0,
                0.0,
                sink,
                te,
            )

        # if rain < 0, borrow from cloud water
        if qrain < 0.0:
            sink = min(-qrain, max(0, qliquid))
            qrain += sink
            qliquid -= sink

        # if cloud water < 0, borrow from vapor
        if qliquid < 0.0:
            sink = min(-qliquid, max(0, qvapor))
            condensation += (sink * delp) * convt * ntimes
            (
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                cvm,
                temp,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            ) = update_hydrometeors_and_temperatures(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                -sink,
                sink,
                0.0,
                0.0,
                0.0,
                0.0,
                te,
            )

    with computation(BACKWARD):
        with interval(1, 2):
            if qvapor[0, 0, -1] < 0:  # top level is negative
                # reduce level 1 by that amount to compensate:
                qvapor = qvapor + qvapor[0, 0, -1] * delp[0, 0, -1] / delp
        with interval(0, 1):
            if qvapor < 0.0:
                qvapor = 0.0  # top level is now 0

    with computation(FORWARD):
        with interval(1, -1):
            # if we borrowed from this level to fix the upper level, account for that:
            if lower_fix[0, 0, -1] != 0:
                qvapor += lower_fix[0, 0, -1]
            if qvapor < 0:  # If negative borrow from below
                lower_fix = qvapor * delp / delp[0, 0, 1]
                qvapor = 0
        with interval(-1, None):
            # if we borrowed from this level to fix the upper level, account for that:
            if lower_fix[0, 0, -1] != 0:
                qvapor += lower_fix[0, 0, -1]

    with computation(BACKWARD):
        with interval(-1, None):
            if (qvapor < 0) and (qvapor[0, 0, -1] > 0):
                # If we can, fix the bottom layer by borrowing from above
                upper_fix = min(-qvapor * delp, qvapor[0, 0, -1] * delp[0, 0, -1])
                qvapor += upper_fix / delp
        with interval(-2, -1):
            if upper_fix[0, 0, 1] != 0.0:
                qvapor -= upper_fix / delp


@gtscript.function
def melt_cloud_ice(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temp,
    cvm,
    te,
    lcpk,
    icpk,
    tcpk,
    tcp3,
):
    """
    Cloud ice melting to form cloud water and rain
    Fortran name is pimlt
    """

    from __externals__ import ql_mlt, tau_imlt, tice_mlt, timestep

    fac_imlt = 1.0 - exp(-timestep / tau_imlt)
    tc = temp - tice_mlt
    if (tc > 0.0) and (qice > QCMIN):
        sink = fac_imlt * tc / icpk
        sink = min(qice, sink)
        tmp = min(sink, positive_diff(ql_mlt, qliquid))

        (
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            cvm,
            temp,
            lcpk,
            icpk,
            tcpk,
            tcp3,
        ) = update_hydrometeors_and_temperatures(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            0.0,
            tmp,
            sink - tmp,
            -sink,
            0,
            0,
            te,
        )
    return (
        qvapor,
        qliquid,
        qrain,
        qice,
        qsnow,
        qgraupel,
        temp,
        cvm,
        lcpk,
        icpk,
        tcpk,
        tcp3,
    )


@gtscript.function
def complete_freeze(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temp,
    cvm,
    te,
    lcpk,
    icpk,
    tcpk,
    tcp3,
):
    """
    Enforces complete freezing below t_wfr
    Fortran name is pcomp
    """
    from __externals__ import t_wfr

    tc = t_wfr - temp
    if (tc > 0.0) and (qliquid > QCMIN):
        sink = qliquid * tc / DT_FR
        sink = min(qliquid, sink, tc / icpk)
        (
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            cvm,
            temp,
            lcpk,
            icpk,
            tcpk,
            tcp3,
        ) = update_hydrometeors_and_temperatures(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            0.0,
            -sink,
            0.0,
            sink,
            0,
            0,
            te,
        )
    return (
        qvapor,
        qliquid,
        qrain,
        qice,
        qsnow,
        qgraupel,
        temp,
        cvm,
        lcpk,
        icpk,
        tcpk,
        tcp3,
    )


@gtscript.function
def cloud_condensation_evaporation(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temp,
    delp,
    density,
    te,
    tcp3,
):
    """
    cloud water condensation and evaporation, Hong and Lim (2006)
    pcond_pevap in Fortran
    """

    from __externals__ import (
        do_cond_timescale,
        rh_fac,
        rhc_cevap,
        tau_l2v,
        tau_v2l,
        timestep,
        use_rhc_cevap,
    )

    fac_l2v = 1.0 - exp(-timestep / tau_l2v)
    fac_v2l = 1.0 - exp(-timestep / tau_v2l)

    qsw, dqdt = tables.sat_spec_hum_water(temp, density)
    qpz = qvapor + qliquid + qice
    rh_tem = qpz / qsw
    dq = qsw - qvapor

    if dq > 0.0:
        fac = min(1.0, fac_l2v * (rh_fac * dq / qsw))
        sink = min(qliquid, fac * dq / (1 + tcp3 * dqdt))
        if (use_rhc_cevap is True) and (rh_tem > rhc_cevap):
            sink = 0.0
        reevaporation = sink * delp
    elif do_cond_timescale:
        fac = min(1.0, fac_v2l * (rh_fac * (-dq) / qsw))
        sink = -min(qvapor, fac * (-dq) / (1.0 + tcp3 * dqdt))
        condensation = -sink * delp
    else:
        sink = -min(qvapor, -dq / (1.0 + tcp3 * dqdt))
        condensation = -sink * delp

    (
        qvapor,
        qliquid,
        qrain,
        qice,
        qsnow,
        qgraupel,
        cvm,
        temp,
        lcpk,
        icpk,
        tcpk,
        tcp3,
    ) = update_hydrometeors_and_temperatures(
        qvapor,
        qliquid,
        qrain,
        qice,
        qsnow,
        qgraupel,
        sink,
        -sink,
        0.0,
        0.0,
        0.0,
        0.0,
        te,
    )

    return (
        qvapor,
        qliquid,
        qrain,
        qice,
        qsnow,
        qgraupel,
        temp,
        cvm,
        lcpk,
        icpk,
        tcpk,
        tcp3,
        reevaporation,
        condensation,
    )


@gtscript.function
def freeze_cloud_water():
    """
    Cloud water freezing to form cloud ice and snow, Lin et al. (1983)
    Fortran name is pifr
    """
    pass


def fast_microphysics(
    temp: FloatField,
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    delp: FloatField,
    density: FloatField,
    cloud_condensation_nuclei: FloatField,
    cloud_ice_nuclei: FloatField,
    condensation: FloatFieldIJ,
    deposition: FloatFieldIJ,
    evaporation: FloatFieldIJ,
    sublimation: FloatFieldIJ,
):
    from __externals__ import convt, do_warm_rain_mp

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
        ) = calc_heat_cap_and_latent_heat_coeff
        (qvapor, qliquid, qrain, qice, qsnow, qgraupel, temp)

        if __INLINED(do_warm_rain_mp is False):
            (
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                temp,
                cvm,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            ) = melt_cloud_ice(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                temp,
                cvm,
                te,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            )

            (
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                temp,
                cvm,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            ) = complete_freeze(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                temp,
                cvm,
                te,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            )

        (
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temp,
            cvm,
            lcpk,
            icpk,
            tcpk,
            tcp3,
            reevap,
            cond,
        ) = cloud_condensation_evaporation(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temp,
            delp,
            density,
            te,
            tcp3,
        )

    with computation(FORWARD), interval(...):
        condensation += cond * convt
        evaporation += reevap * convt

    with computation(PARALLEL), interval(...):
        if __INLINED(do_warm_rain_mp is False):
            pass
        pass


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
        self._idx: GridIndexing = stencil_factory.grid_indexing
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

        self._n_min = 1600
        self._delt = 0.1
        self._tice = 273.15
        self._tice_mlt = self.namelist.tice
        self._esbasw = 1013246.0
        self._tbasw = self._tice + 100.0
        self._esbasi = 6107.1
        self._tmin = self._tice - self._n_min * self._delt
        if self.namelist.do_warm_rain:  # unsupported
            self._t_wfr = self._tmin
        else:
            self._t_wfr = self._tice - 40.0

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

        self._set_timestepping(self.namelist.dt_atmos)  # will change from dt_atmos
        # when we inline microphysics

        self._convert_mm_day = 86400.0 * constants.RGRAV / self._split_timestep

        self._copy_stencil = stencil_factory.from_origin_domain(
            basic.copy_defn,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._calc_total_energy = stencil_factory.from_origin_domain(
            func=calc_total_energy,
            externals={
                "hydrostatic": self.namelist.hydrostatic,
                "c_air": self._c_air,
                "c1_vapor": self._c1_vap,
                "c1_liquid": self._c1_liquid,
                "c1_ice": self._c1_ice,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        if self.namelist.consv_checker is True:
            self._moist_total_energy_and_water_mq = stencil_factory.from_origin_domain(
                func=moist_total_energy_and_water,
                externals={
                    "hydrostatic": self.namelist.hydrostatic,
                    "moist_q": True,
                    "c1_vapor": self._c1_vap,
                    "c1_liquid": self._c1_liquid,
                    "c1_ice": self._c1_ice,
                    "lv00": self._lv00,
                    "li00": self._li00,
                    "c_air": self._c_air,
                    "timestep": self._full_timestep,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

            self._moist_total_energy_and_water = stencil_factory.from_origin_domain(
                func=moist_total_energy_and_water,
                externals={
                    "hydrostatic": self.namelist.hydrostatic,
                    "moist_q": False,
                    "c1_vapor": self._c1_vap,
                    "c1_liquid": self._c1_liquid,
                    "c1_ice": self._c1_ice,
                    "lv00": self._lv00,
                    "li00": self._li00,
                    "c_air": self._c_air,
                    "timestep": self._full_timestep,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

            self._calc_sedimentation_energy_loss = stencil_factory.from_origin_domain(
                func=calc_sedimentation_energy_loss,
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        self._calc_total_energy = stencil_factory.from_origin_domain(
            func=calc_total_energy,
            externals={
                "hydrostatic": self.namelist.hydrostatic,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._convert_specific_to_mass_mixing_ratios_and_calculate_densities = (
            stencil_factory.from_origin_domain(
                func=convert_specific_to_mass_mixing_ratios_and_calculate_densities,
                externals={
                    "inline_mp": self.namelist.do_inline_mp,
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
                "prog_ccn": self.namelist.prog_ccn,
                "ccn_l": self.namelist.ccn_l,
                "ccn_o": self.namelist.ccn_o,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._subgrid_deviation_and_relative_humidity = (
            stencil_factory.from_origin_domain(
                func=subgrid_deviation_and_relative_humidity,
                externals={
                    "dw_land": self.namelist.dw_land,
                    "dw_ocean": self.namelist.dw_ocean,
                    "rh_inc": self.namelist.rh_inc,
                    "rh_inr": self.namelist.rh_inr,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )
        )

        if self.namelist.fix_negative:
            self._adjust_negative_tracers = stencil_factory.from_origin_domain(
                func=adjust_negative_tracers,
                externals={
                    "c1_vapor": self._c1_vap,
                    "c1_liq": self._c1_liquid,
                    "c1_ice": self._c1_ice,
                    "lv00": self._lv00,
                    "li00": self._li00,
                    "li20": self._li20,
                    "d1_vap": self._d1_vap,
                    "d1_ice": self._d1_ice,
                    "tice": self._tice,
                    "t_wfr": self._t_wfr,
                    "convt": self._convert_mm_day,
                    "ntimes": self._ntimes,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

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
