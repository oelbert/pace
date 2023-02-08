import math

import table_functions as tables
from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import (
    __INLINED,
    BACKWARD,
    FORWARD,
    PARALLEL,
    computation,
    exp,
    interval,
    log,
    sqrt,
)
from moist_total_energy import calc_total_energy

import pace.fv3core.stencils.basic_operations as basic
import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import Float, FloatField, FloatFieldIJ
from pace.util import X_DIM, Y_DIM, Z_DIM
from pace.util.grid import GridData

from ..._config import PhysicsConfig


QCMIN = 1.0e-15  # min value for cloud condensates (kg/kg)
QFMIN = 1.0e-8  # min value for sedimentation (kg/kg)
DT_FR = 8.0  # t_wfr - dt_fr: minimum temperature water can exist
# (Moore and Molinero 2011)

DZ_MIN_FLIP = 1.0e-2  # used for correcting flipped height (m)


@gtscript.function
def positive_diff(x, y):
    """
    equivalent of Fortran dim
    """
    return max(x - y, 0.0)


@gtscript.function
def calc_particle_concentration(tracer, density, pca, pcb, mu):
    """
    pc Part of cal_pc_ed_oe_rr_tv in Fortran
    """

    return pca / pcb * exp(mu / (mu + 3) * log(6 * density * tracer))


@gtscript.function
def calc_effective_diameter(tracer, density, eda, edb, mu):
    """
    ed Part of cal_pc_ed_oe_rr_tv in Fortran
    """

    return eda / edb * exp(1.0 / (mu + 3) * log(6 * density * tracer))


@gtscript.function
def calc_optical_extinction(tracer, density, oea, oeb, mu):
    """
    oe Part of cal_pc_ed_oe_rr_tv in Fortran
    """

    return oea / oeb * exp((mu + 2) / (mu + 3) * log(6 * density * tracer))


@gtscript.function
def calc_radar_reflectivity(tracer, density, rra, rrb, mu):
    """
    rr Part of cal_pc_ed_oe_rr_tv in Fortran
    """

    return rra / rrb * exp((mu + 6) / (mu + 3) * log(6 * density * tracer))


@gtscript.function
def calc_terminal_velocity(tracer, density, tva, tvb, mu, blin):
    """
    mass-weighted terminal velocity
    tv Part of cal_pc_ed_oe_rr_tv in Fortran
    """

    return tva / tvb * exp(blin / (mu + 3) * log(6 * density * tracer))


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
def freeze_cloud_water(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temp,
    density,
    cvm,
    te,
    lcpk,
    icpk,
    tcpk,
    tcp3,
):
    """
    Cloud water freezing to form cloud ice and snow, Lin et al. (1983)
    Fortran name is pifr
    """

    from __externals__ import qi0_crt, t_wfr

    tc = t_wfr - temp
    if (tc > 0.0) and (qliquid > QCMIN):
        sink = qliquid * tc / DT_FR
        sink = min(qliquid, sink, tc / icpk)
        qim = qi0_crt / density
        tmp = min(sink, positive_diff(qim, qice))

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
            tmp,
            sink - tmp,
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
        cvm,
        temp,
        lcpk,
        icpk,
        tcpk,
        tcp3,
    )


@gtscript.function
def wegener_bergeron_findeisen(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temp,
    density,
    cvm,
    te,
    lcpk,
    icpk,
    tcpk,
    tcp3,
):
    """
    Wegener Bergeron Findeisen process, Storelvmo and Tan (2015)
    Fortran name is pwbf
    """

    from __externals__ import qi0_crt, tau_wbf, tice, timestep

    tc = tice - temp
    qsw, dqdt = tables.sat_spec_hum_water(temp, density)
    qsi, dqdt = tables.sat_spec_hum_water_ice(temp, density)

    if (
        (tc > 0.0)
        and (qliquid > QCMIN)
        and (qice > QCMIN)
        and (qvapor > qsi)
        and (qvapor < qsw)
    ):
        fac_wbf = 1.0 - exp(-timestep / tau_wbf)
        sink = min(fac_wbf * qliquid, tc / icpk)
        qim = qi0_crt / density
        tmp = min(sink, positive_diff(qim, qice))

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
            tmp,
            sink - tmp,
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
        cvm,
        temp,
        lcpk,
        icpk,
        tcpk,
        tcp3,
    )


@gtscript.function
def freeze_bigg(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    cloud_condensation_nuclei,
    temp,
    density,
    cvm,
    te,
    lcpk,
    icpk,
    tcpk,
    tcp3,
):
    """
    Bigg freezing mechanism, Bigg (1953)
    Fortran name is pbigg
    """
    from __externals__ import do_psd_water_num, muw, pcaw, pcbw, tice, timestep

    tc = tice - temp
    if (tc > 0.0) and (qliquid > QCMIN):
        if do_psd_water_num is True:
            cloud_condensation_nuclei = calc_particle_concentration(
                qliquid, density, pcaw, pcbw, muw
            )
            cloud_condensation_nuclei = cloud_condensation_nuclei / density

        sink = (
            100.0
            / (constants.RHO_W * cloud_condensation_nuclei)
            * timestep
            * (exp(0.66 * tc) - 1.0)
            * qliquid ** 2.0
        )
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
        cloud_condensation_nuclei,
        cvm,
        temp,
        lcpk,
        icpk,
        tcpk,
        tcp3,
    )


@gtscript.function
def freeze_rain_to_graupel_simple(
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
    Fortran name is pgfr_simp
    """

    from __externals__ import tau_r2g, tice, timestep

    fac_r2g = 1.0 - exp(-timestep / tau_r2g)
    tc = temp - tice
    if (tc < 0.0) and (qrain > QCMIN):
        sink = (-tc * 0.025) ** 2 * qrain
        sink = min(qrain, sink, -fac_r2g * tc / icpk)

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

    return (
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
    )


@gtscript.function
def melt_snow_simple(
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
    Fortran name is psmlt_simp
    """

    from __externals__ import qs_mlt, tau_smlt, tice, timestep

    fac_smlt = 1.0 - exp(-timestep / tau_smlt)

    tc = temp - tice
    if (tc > 0.0) and (qsnow > QCMIN):
        sink = (tc * 0.1) ** 2 * qsnow
        sink = min(qsnow, sink, fac_smlt * tc / icpk)
        tmp = min(sink, positive_diff(qs_mlt, qliquid))

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
            0.0,
            sink,
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
        cvm,
        temp,
        lcpk,
        icpk,
        tcpk,
        tcp3,
    )


@gtscript.function
def autoconvert_water_to_rain_simple(
    qliquid,
    qrain,
    temp,
):
    """
    Cloud water to rain autoconversion, simple version
    Fortran name is praut_simp
    """

    from __externals__ import ql0_max, t_wfr, tau_l2r, timestep

    fac_l2r = 1.0 - exp(-timestep / tau_l2r)

    tc = temp - t_wfr
    if (tc > 0) and (qliquid > ql0_max):
        sink = fac_l2r * (qliquid - ql0_max)
        qliquid -= sink
        qrain += sink

    return qliquid, qrain


@gtscript.function
def deposit_and_sublimate_ice(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    cloud_ice_nuclei,
    temp,
    delp,
    density,
    cvm,
    te,
    lcpk,
    icpk,
    tcpk,
    tcp3,
):
    """
    Cloud ice deposition and sublimation, Hong et al. (2004)
    Fortran name is pidep_pisub
    """

    from __externals__ import (
        do_psd_ice_num,
        igflag,
        inflag,
        is_fac,
        mui,
        pcai,
        pcbi,
        prog_ccn,
        qi_lim,
        t_sub,
        tice,
        timestep,
    )

    sub = 0.0
    dep = 0.0

    if temp < tice:
        pidep = 0
        qsi, dqdt = tables.sat_spec_hum_water_ice(temp, density)
        dq = qvapor - qsi
        tmp = dq / (1.0 + tcpk * dqdt)

        if qice > QCMIN:
            if prog_ccn is False:
                if inflag == 1:
                    cloud_ice_nuclei = 5.38e7 * exp(0.75 * log(qice * density))
                elif inflag == 2:
                    cloud_ice_nuclei = exp(-2.80 + 0.262 * (tice - temp)) * 1000.0
                elif inflag == 3:
                    cloud_ice_nuclei = (
                        exp(-0.639 + 12.96 * (qvapor / qsi - 1.0)) * 1000.0
                    )
                elif inflag == 4:
                    cloud_ice_nuclei = 5.0e-3 * exp(0.304 * (tice - temp)) * 1000.0
                elif inflag == 5:
                    cloud_ice_nuclei = 1.0e-5 * exp(0.5 * (tice - temp)) * 1000.0
                else:
                    raise ValueError(
                        f"Ice Nucleation Flag must be an integer from 1 to 5"
                        f"not {inflag}"
                    )
            if do_psd_ice_num:
                cloud_ice_nuclei = calc_particle_concentration(
                    qice, density, pcai, pcbi, mui
                )
                cloud_ice_nuclei = cloud_ice_nuclei / density

            pidep = (
                timestep
                * dq
                * 4.0
                * 11.9
                * exp(0.5 * log(qice * density * cloud_ice_nuclei))
                / (
                    qsi
                    * density
                    * (tcpk * cvm) ** 2
                    / (constants.TCOND * constants.RVGAS * temp ** 2)
                    + 1.0 / constants.VDIFU
                )
            )
        if dq > 0:
            tc = tice - temp
            qi_gen = 4.92e-11 * exp(1.33 * log(1.0e3 * exp(0.1 * tc)))
            if igflag == 1:
                qi_crt = qi_gen / density
            elif igflag == 2:
                qi_crt = qi_gen * min(qi_lim, 0.1 * tc) / density
            elif igflag == 3:
                qi_crt = 1.82e-6 * min(qi_lim, 0.1 * tc) / density
            elif igflag == 4:
                qi_crt = max(qi_gen, 1.82e-6) * min(qi_lim, 0.1 * tc) / density
            else:
                raise ValueError(
                    f"Ice Generation Flag must be an integer from 1 to 4 not {igflag}"
                )
            sink = min(tmp, max(qi_crt - qice, pidep), tc / tcpk)
            dep = sink * delp
        else:
            pidep = pidep * min(1, positive_diff(temp, t_sub) * is_fac)
            sink = max(pidep, tmp - qice)
            sub = -sink * delp

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
            0.0,
            0.0,
            sink,
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
        cvm,
        temp,
        lcpk,
        icpk,
        tcpk,
        tcp3,
        dep,
        sub,
    )


@gtscript.function
def autoconvert_ice_to_snow_simple(qice, qsnow, temp, density):
    """
    Cloud ice to snow autoconversion, simple version
    Fortran name is psaut_simp
    """
    from __externals__ import qi0_max, tau_i2s, tice, timestep

    fac_i2s = 1.0 - exp(-timestep / tau_i2s)

    tc = temp - tice
    qim = qi0_max / density

    if (tc < 0) and (qice > qim):
        sink = fac_i2s * (qice - qim)
        qice -= sink
        qsnow += sink
    return qice, qsnow


def fast_microphysics(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    temp: FloatField,
    delp: FloatField,
    density: FloatField,
    cloud_condensation_nuclei: FloatField,
    cloud_ice_nuclei: FloatField,
    condensation: FloatFieldIJ,
    deposition: FloatFieldIJ,
    evaporation: FloatFieldIJ,
    sublimation: FloatFieldIJ,
):
    from __externals__ import convt, do_warm_rain_mp, do_wbf

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

        if __INLINED(do_warm_rain_mp is False):
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
            ) = freeze_cloud_water(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                temp,
                density,
                cvm,
                te,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            )
            if __INLINED(do_wbf is True):
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
                ) = wegener_bergeron_findeisen(
                    qvapor,
                    qliquid,
                    qrain,
                    qice,
                    qsnow,
                    qgraupel,
                    temp,
                    density,
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
                cloud_condensation_nuclei,
                cvm,
                temp,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            ) = freeze_bigg(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                cloud_condensation_nuclei,
                temp,
                density,
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
                cvm,
                temp,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            ) = freeze_rain_to_graupel_simple(
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
                cvm,
                temp,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            ) = melt_snow_simple(
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

        qliquid, qrain = autoconvert_water_to_rain_simple(qliquid, qrain, temp)

        if __INLINED(do_warm_rain_mp is False):
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
                dep,
                sub,
            ) = deposit_and_sublimate_ice(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                cloud_ice_nuclei,
                temp,
                delp,
                density,
                cvm,
                te,
                lcpk,
                icpk,
                tcpk,
                tcp3,
            )

            qice, qsnow = autoconvert_ice_to_snow_simple(qice, qsnow, temp, density)

    with computation(FORWARD), interval(...):
        if __INLINED(do_warm_rain_mp is False):
            condensation += cond * convt
            evaporation += reevap * convt
            deposition += dep * convt
            sublimation += sub * convt


@gtscript.function
def calc_terminal_velocity_rsg(
    q, density, density_factor, const_v, v_fac, tva, tvb, mu, blin, v_max
):
    """
    Calculate terminal velocity for rain, snow, and graupel, Lin et al. (1983)
    Fortran name is term_rsg
    """
    if const_v is True:
        v_terminal = v_fac
    else:
        if q < QFMIN:
            v_terminal = 0.0
        else:
            v_terminal = calc_terminal_velocity(q, density, tva, tvb, mu, blin)
            v_terminal = v_fac * v_terminal * density_factor
            v_terminal = min(v_max, max(0.0, v_terminal))

    return v_terminal


@gtscript.function
def calc_terminal_velocity_ice(qice, temp, density, constant_v, v_fac, v_max):
    """
    Fortran name is term_ice
    """

    from __externals__ import ifflag, tice

    aa = -4.14122e-5
    bb = -0.00538922
    cc = -0.0516344
    dd = 0.00216078
    ee = 1.9714

    if constant_v is True:
        v_terminal = v_fac
    else:
        if qice < QFMIN:
            v_terminal = 0.0
        else:
            tc = temp - tice
            if ifflag == 1:
                v_terminal = (
                    (3.0 + (log(qice * density) / log(10))) * (tc * (aa * tc + bb) + cc)
                    + dd * tc
                    + ee
                )
                v_terminal + 0.01 * v_fac * exp(v_terminal * log(10.0))
            elif ifflag == 2:
                v_terminal = v_fac * 3.29 * exp(0.16 * log(qice * density))
            else:
                raise ValueError(f"Ice Formation Flag must be 1 or 2 not {ifflag}")
            v_terminal = min(v_max, max(0.0, v_terminal))

    return v_terminal


def sedimentation(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    ua: FloatField,
    va: FloatField,
    wa: FloatField,
    temp: FloatField,
    delp: FloatField,
    delz: FloatField,
    density: FloatField,
    density_factor: FloatField,
    column_energy_change: FloatFieldIJ,
):
    """
    Sedimentation of cloud ice, snow, graupel or hail, and rain
    """

    from __externals__ import do_psd_water_fall  # noqa
    from __externals__ import timestep  # noqa
    from __externals__ import (
        blini,
        const_vi,
        do_psd_ice_fall,
        do_sedi_melt,
        mui,
        tvai,
        tvbi,
        vi_fac,
        vi_max,
    )

    with computation(PARALLEL), interval(...):
        w1 = 0.0
        r1 = 0.0
        i1 = 0.0
        s1 = 0.0
        g1 = 0.0

        vtw = 0.0
        vtr = 0.0
        vti = 0.0
        vts = 0.0
        vtg = 0.0

        pfw = 0.0
        pfr = 0.0
        pfi = 0.0
        pfs = 0.0
        pfg = 0.0

        (
            q_liq,
            q_solid,
            cvm,
            te,
            lcpk,
            icpk,
            tcpk,
            tcp3,
        ) = calc_heat_cap_and_latent_heat_coeff(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temp,
        )

        # terminal fall and melting of falling cloud ice into rain
        if __INLINED(do_psd_ice_fall is True):
            v_terminal_ice = calc_terminal_velocity_rsg(
                qice,
                density,
                density_factor,
                const_vi,
                vi_fac,
                tvai,
                tvbi,
                mui,
                blini,
                vi_max,
            )
        else:
            v_terminal_ice = calc_terminal_velocity_ice(
                qice, temp, density, const_vi, vi_fac, vi_max
            )

        if __INLINED(do_sedi_melt is True):
            pass


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
                zt = zt[0, 0, -1] - DZ_MIN_FLIP


def sedi_melt(qice: FloatField):
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
        do_mp_fast: bool = False,
    ):
        self.namelist = namelist
        self._idx: GridIndexing = stencil_factory.grid_indexing
        self._max_timestep = self.namelist.mp_time
        self._ntimes = self.namelist.ntimes
        self._do_mp_fast = do_mp_fast

        # atmospheric physics constants
        # TODO these should live in a separate object that can be passed
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

        self._calculate_particle_parameters()

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

        self._convert_mm_day = 86400.0 * constants.RGRAV / self.split_timestep

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
                    "timestep": self.full_timestep,
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

        if self._do_mp_fast:
            self._fast_microphysics = stencil_factory.from_origin_domain(
                func=fast_microphysics,
                externals={
                    "do_warm_rain_mp": self.namelist.do_warm_rain,
                    "do_wbf": self.namelist.do_wbf,
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
                    "ql_mlt": self.namelist.ql_mlt,
                    "qs_mlt": self.namelist.qs_mlt,
                    "tau_imlt": self.namelist.tau_imlt,
                    "tice_mlt": self._tice_mlt,
                    "timestep": self.full_timestep,
                    "do_cond_timescale": self.namelist.do_cond_timescale,
                    "rh_fac": self.namelist.rh_fac,
                    "rhc_cevap": self.namelist.rhc_cevap,
                    "tau_l2v": self.namelist.tau_l2v,
                    "tau_v2l": self.namelist.tau_v2l,
                    "tau_r2g": self.namelist.tau_r2g,
                    "tau_smlt": self.namelist.tau_smlt,
                    "tau_l2r": self.namelist.tau_l2r,
                    "use_rhc_cevap": self.namelist.use_rhc_cevap,
                    "qi0_crt": self.namelist.qi0_crt,
                    "qi0_max": self.namelist.qi0_max,
                    "ql0_max": self.namelist.ql0_max,
                    "tau_wbf": self.namelist.tau_wbf,
                    "do_psd_water_num": self.namelist.do_psd_water_num,
                    "do_psd_ice_num": self.namelist.do_psd_ice_num,
                    "muw": self.namelist.muw,
                    "mui": self.namelist.mui,
                    "pcaw": self._pcaw,
                    "pcbw": self._pcbw,
                    "pcai": self._pcai,
                    "pcbi": self._pcbi,
                    "prog_ccn": self.namelist.prog_ccn,
                    "inflag": self.namelist.inflag,
                    "igflag": self.namelist.igflag,
                    "qi_lim": self.namelist.qi_lim,
                    "t_sub": self.namelist.t_sub,
                    "is_fac": self.namelist.is_fac,
                    "tau_i2s": self.namelist.tau_i2s,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        pass

    def gfdl_cloud_microphys_init(self, dt_atmos: float):
        pass

    def _update_timestep_if_needed(self, timestep: float):
        if timestep != self.full_timestep:
            self._set_timestepping(timestep)

    def _set_timestepping(self, full_timestep: float):
        self._ntimes = int(
            max(self._ntimes, full_timestep / min(full_timestep, self._max_timestep))
        )
        self.split_timestep = full_timestep / self._ntimes
        self.full_timestep = full_timestep
        pass

    def _calculate_particle_parameters(self):
        """
        Calculate parameters for particle concentration, effective diameter,
        optical extinction, radar reflectivity, and terminal velocity
        for each tracer species
        """
        muw = self.namelist.muw
        mui = self.namelist.mui
        n0w_exp = self.namelist.n0w_exp
        n0i_exp = self.namelist.n0i_exp
        n0w_sig = self.namelist.n0w_sig
        n0i_sig = self.namelist.n0i_sig
        alini = self.namelist.alini
        blini = self.namelist.blini
        # Particle Concentration:
        self._pcaw = (
            math.exp(3 / (muw + 3) * math.log(n0w_sig))
            * math.gamma(muw)
            * math.exp(3 * n0w_exp / (muw + 3) * math.log(10.0))
        )
        self._pcbw = math.exp(
            muw / (muw + 3) * math.log(math.pi * constants.RHO_W * math.gamma(muw + 3))
        )
        self._pcai = (
            math.exp(3 / (mui + 3) * math.log(n0i_sig))
            * math.gamma(mui)
            * math.exp(3 * n0i_exp / (mui + 3) * math.log(10.0))
        )
        self._pcbi = math.exp(
            mui / (mui + 3) * math.log(math.pi * constants.RHO_I * math.gamma(mui + 3))
        )

        # Effective Diameter

        # Optical Extinction

        # Radar Reflectivity

        # Terminal Velocity
        self._tvai = (
            math.exp(blini / (mui + 3) * math.log(n0i_sig))
            * alini
            * math.gamma(mui + blini + 3)
            * math.exp(-blini * n0i_exp / (mui + 3) * math.log(10.0))
        )
        self._tvbi = math.exp(
            blini
            / (mui + 3)
            * math.log(math.pi * constants.RHO_I * math.gamma(mui + 3))
        ) * math.gamma(mui + 3)

    def __call__(self, state: MicrophysicsState, timestep: float):
        pass
