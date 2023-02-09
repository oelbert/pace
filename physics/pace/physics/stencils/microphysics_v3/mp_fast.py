import physical_functions as physfun
from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import (
    __INLINED,
    FORWARD,
    PARALLEL,
    computation,
    exp,
    interval,
    log,
)

import pace.fv3core.stencils.basic_operations as basic
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ

from ..._config import FastMPConfig


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
    if (tc > 0.0) and (qice > constants.QCMIN):
        sink = fac_imlt * tc / icpk
        sink = min(qice, sink)
        tmp = min(sink, basic.dim(ql_mlt, qliquid))

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
        ) = physfun.update_hydrometeors_and_temperatures(
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
    if (tc > 0.0) and (qliquid > constants.QCMIN):
        sink = qliquid * tc / constants.DT_FR
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
        ) = physfun.update_hydrometeors_and_temperatures(
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

    qsw, dqdt = physfun.sat_spec_hum_water(temp, density)
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
    ) = physfun.update_hydrometeors_and_temperatures(
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
    if (tc > 0.0) and (qliquid > constants.QCMIN):
        sink = qliquid * tc / constants.DT_FR
        sink = min(qliquid, sink, tc / icpk)
        qim = qi0_crt / density
        tmp = min(sink, basic.dim(qim, qice))

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
        ) = physfun.update_hydrometeors_and_temperatures(
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
    qsw, dqdt = physfun.sat_spec_hum_water(temp, density)
    qsi, dqdt = physfun.sat_spec_hum_water_ice(temp, density)

    if (
        (tc > 0.0)
        and (qliquid > constants.QCMIN)
        and (qice > constants.QCMIN)
        and (qvapor > qsi)
        and (qvapor < qsw)
    ):
        fac_wbf = 1.0 - exp(-timestep / tau_wbf)
        sink = min(fac_wbf * qliquid, tc / icpk)
        qim = qi0_crt / density
        tmp = min(sink, basic.dim(qim, qice))

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
        ) = physfun.update_hydrometeors_and_temperatures(
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
    if (tc > 0.0) and (qliquid > constants.QCMIN):
        if do_psd_water_num is True:
            cloud_condensation_nuclei = physfun.calc_particle_concentration(
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
        ) = physfun.update_hydrometeors_and_temperatures(
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
    if (tc < 0.0) and (qrain > constants.QCMIN):
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
        ) = physfun.update_hydrometeors_and_temperatures(
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
    if (tc > 0.0) and (qsnow > constants.QCMIN):
        sink = (tc * 0.1) ** 2 * qsnow
        sink = min(qsnow, sink, fac_smlt * tc / icpk)
        tmp = min(sink, basic.dim(qs_mlt, qliquid))

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
        ) = physfun.update_hydrometeors_and_temperatures(
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
        qsi, dqdt = physfun.sat_spec_hum_water_ice(temp, density)
        dq = qvapor - qsi
        tmp = dq / (1.0 + tcpk * dqdt)

        if qice > constants.QCMIN:
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
                cloud_ice_nuclei = physfun.calc_particle_concentration(
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
            pidep = pidep * min(1, basic.dim(temp, t_sub) * is_fac)
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
        ) = physfun.update_hydrometeors_and_temperatures(
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
        ) = physfun.calc_heat_cap_and_latent_heat_coeff
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


class FastMicrophysics:
    """
    Fast microphysics loop, mp_fast in Fortran
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        config: FastMPConfig,
        timestep: float,
        ntimes: int,
        convert_mm_day: float,
    ):

        self._idx: GridIndexing = stencil_factory.grid_indexing

        self._fast_microphysics = stencil_factory.from_origin_domain(
            func=fast_microphysics,
            externals={
                "do_warm_rain_mp": config.do_warm_rain,
                "do_wbf": config.do_wbf,
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
                "convt": convert_mm_day,
                "ntimes": ntimes,
                "ql_mlt": config.ql_mlt,
                "qs_mlt": config.qs_mlt,
                "tau_imlt": config.tau_imlt,
                "tice_mlt": config.tice_mlt,
                "timestep": timestep,
                "do_cond_timescale": config.do_cond_timescale,
                "rh_fac": config.rh_fac,
                "rhc_cevap": config.rhc_cevap,
                "tau_l2v": config.tau_l2v,
                "tau_v2l": config.tau_v2l,
                "tau_r2g": config.tau_r2g,
                "tau_smlt": config.tau_smlt,
                "tau_l2r": config.tau_l2r,
                "use_rhc_cevap": config.use_rhc_cevap,
                "qi0_crt": config.qi0_crt,
                "qi0_max": config.qi0_max,
                "ql0_max": config.ql0_max,
                "tau_wbf": config.tau_wbf,
                "do_psd_water_num": config.do_psd_water_num,
                "do_psd_ice_num": config.do_psd_ice_num,
                "muw": config.muw,
                "mui": config.mui,
                "pcaw": config.pcaw,
                "pcbw": config.pcbw,
                "pcai": config.pcai,
                "pcbi": config.pcbi,
                "prog_ccn": config.prog_ccn,
                "inflag": config.inflag,
                "igflag": config.igflag,
                "qi_lim": config.qi_lim,
                "t_sub": config.t_sub,
                "is_fac": config.is_fac,
                "tau_i2s": config.tau_i2s,
            },
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

        self._fast_microphysics(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temp,
            delp,
            density,
            cloud_condensation_nuclei,
            cloud_ice_nuclei,
            condensation,
            deposition,
            evaporation,
            sublimation,
        )
