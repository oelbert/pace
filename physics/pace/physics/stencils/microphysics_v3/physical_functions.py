from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import exp, floor, log, sqrt

import pace.fv3core.stencils.basic_operations as basic
import pace.util.constants as constants


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


@gtscript.function
def accretion_2d(qden, denfac, c, blin, mu):
    return denfac * c * exp((2 + mu + blin) / (mu + 3) * log(6 * qden))


@gtscript.function
def accretion_3d(
    v1,
    v2,
    q1,
    q2,
    density,
    c,
    acco1,
    acco2,
    acco3,
    acc1,
    acc2,
):
    """
    Accretion function, Lin et al. (1983)
    Fortran name is acr3d
    """
    from __externals__ import vdiffflag

    t1 = exp(1.0 / (acc1 + 3) * log(6 * q1 * density))
    t2 = exp(1.0 / (acc2 + 3) * log(6 * q2 * density))

    if vdiffflag == 1:
        vdiff = abs(v1 - v2)
    elif vdiffflag == 2:
        vdiff = sqrt((1.20 * v1 - 0.95 * v2) ** 2.0 + 0.08 * v1 * v2)
    else:  # vdiffflag == 3:
        vdiff = sqrt((1.00 * v1 - 1.00 * v2) ** 2.0 + 0.04 * v1 * v2)

    accrete = c * vdiff / density
    tmp = acco1 * exp((6 + acc1 - 1) * log(t1)) * exp((acc2 + 1 - 1) * log(t2))
    tmp += acco2 * exp((6 + acc1 - 2) * log(t1)) * exp((acc2 + 2 - 1) * log(t2))
    tmp += acco3 * exp((6 + acc1 - 3) * log(t1)) * exp((acc2 + 3 - 1) * log(t2))

    return accrete * tmp


@gtscript.function
def calc_heat_cap_and_latent_heat_coeff(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temperature,
):
    """
    Fortran name is cal_mhc_lhc
    """

    from __externals__ import (
        c1_ice,
        c1_liq,
        c1_vap,
        d1_ice,
        d1_vap,
        li00,
        li20,
        lv00,
        t_wfr,
    )

    q_liq = qliquid + qrain
    q_solid = qice + qsnow + qgraupel
    cvm = 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
    te = cvm * temperature + lv00 * qvapor - li00 * q_solid
    lcpk = (lv00 + d1_vap * temperature) / cvm
    icpk = (li00 + d1_ice * temperature) / cvm
    tcpk = (li20 + (d1_vap + d1_ice) * temperature) / cvm
    tcp3 = lcpk + icpk * min(
        1.0, basic.dim(constants.TICE0, temperature) / (constants.TICE0 - t_wfr)
    )

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
        c1_vap,
        d1_ice,
        d1_vap,
        li00,
        li20,
        lv00,
        t_wfr,
    )

    qvapor += delta_vapor
    qliquid += delta_liquid
    qrain += delta_rain
    qice += delta_ice
    qsnow += delta_snow
    qgraupel += delta_graupel

    q_l = qrain + qliquid
    q_solid = qice + qsnow + qgraupel
    cvm = 1.0 + qvapor * c1_vap + q_l * c1_liq + q_solid * c1_ice

    tk = (te - lv00 * qvapor + li00 * (qice + qsnow + qgraupel)) / cvm
    lcpk = (lv00 + d1_vap * tk) / cvm
    icpk = (li00 + d1_ice * tk) / cvm
    tcpk = (li20 + (d1_vap + d1_ice) * tk) / cvm
    tcp3 = lcpk + icpk * min(
        1.0, basic.dim(constants.TICE0, tk) / (constants.TICE0 - t_wfr)
    )

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


@gtscript.function
def table0(temp):
    """
    Saturation water vapor pressure table 0, water only
    useful for idealized experiments
    it can also be used in warm rain microphyscis only
    """
    return constants.E00 * exp(
        (
            constants.DC_VAP0 * log(temp / constants.TICE0)
            + constants.LV0_0 * (temp - constants.TICE0) / (temp * constants.TICE0)
        )
        / constants.RVGAS
    )


@gtscript.function
def table2(temp):
    """
    Saturation water vapor pressure table 2, water and ice
    same as table 1, but the blending is replaced with smoothing around 0 deg C
    it is not designed for mixed-phase cloud microphysics
    used for ice microphysics (< 0 deg C) or warm rain microphysics (> 0 deg C)
    """
    if temp < constants.TICE0:
        # Over ice between -160 degrees Celsius and 0 degrees Celsius
        return_val = constants.E00 * exp(
            (
                constants.D2ICE0 * log(temp / constants.TICE0)
                + constants.LI2_0 * (temp - constants.TICE0) / (temp * constants.TICE0)
            )
            / constants.RVGAS
        )

    else:
        # Over water between 0 degrees Celsius and 102 degrees Celsius
        return_val = table0(temp)

    return return_val


@gtscript.function
def sat_spec_hum_water(temp, density):
    """
    qs_core with table 0 in microphysics
    compute the saturated specific humidity, core function
    """
    q = table0(temp) / (constants.RVGAS * temp * density)
    dqdt = q * (constants.DC_VAP0 + constants.LV0_0 / temp) / (constants.RVGAS * temp)
    return q, dqdt


@gtscript.function
def sat_spec_hum_water_ice(temperature, density):
    temp = max(constants.TICE0 - 160.0, min(temperature, constants.TICE0 + 102.0))
    q = table2(temp) / (constants.RVGAS * temperature * density)
    if temp < constants.TICE0:
        dqdt = (
            q
            * (constants.D2ICE0 + constants.LI2_0 / temp)
            / (constants.RVGAS * temperature)
        )
    else:
        dqdt = (
            q
            * (constants.DC_VAP0 + constants.LV0_0 / temp)
            / (constants.RVGAS * temperature)
        )
    return q, dqdt


@gtscript.function
def temperature_index(temperature):
    tmin = constants.TICE0 - 160.0
    return floor(10.0 * (temperature - tmin)) / 10.0 + tmin


@gtscript.function
def table0_delta(int_temperature):
    tmax = constants.TICE0 - 160.0 + 262.0
    int_temperature = min(int_temperature, tmax)
    return max(0.0, table0(int_temperature + 0.1) - table0(int_temperature))


@gtscript.function
def lookup_0(temperature):
    int_temperature = temperature_index(temperature)
    return table0(int_temperature) + 10 * (
        temperature - int_temperature
    ) * table0_delta(int_temperature)


@gtscript.function
def table2_delta(int_temperature):
    tmax = constants.TICE0 - 160.0 + 262.0
    int_temperature = min(int_temperature, tmax)
    return max(0.0, table2(int_temperature + 0.1) - table2(int_temperature))


@gtscript.function
def lookup_2(temperature):
    int_temperature = temperature_index(temperature)
    return table2(int_temperature) + 10 * (
        temperature - int_temperature
    ) * table2_delta(int_temperature)


@gtscript.function
def wqs(temperature, density):
    tmin = constants.TICE0 - 160.0
    temp_limit = min(tmin + 262.1, max(tmin, temperature))
    qsat = lookup_0(temp_limit) / (constants.RVGAS * temperature * density)
    it = temperature_index(temp_limit - 0.05)
    dqdt = (
        10.0
        * (
            table0_delta(it)
            + 10 * (temp_limit - it) * (table0_delta(it + 0.1) - table0_delta(it))
        )
        / (constants.RVGAS * temperature * density)
    )
    return qsat, dqdt


@gtscript.function
def iqs(temperature, density):
    tmin = constants.TICE0 - 160.0
    temp_limit = min(tmin + 262.1, max(tmin, temperature))
    qsat = lookup_2(temp_limit) / (constants.RVGAS * temperature * density)
    it = temperature_index(temp_limit - 0.05)
    dqdt = (
        10.0
        * (
            table2_delta(it)
            + 10 * (temp_limit - it) * (table2_delta(it + 0.1) - table2_delta(it))
        )
        / (constants.RVGAS * temperature * density)
    )
    return qsat, dqdt


@gtscript.function
def moist_heat_capacity(qvapor, qliquid, qrain, qice, qsnow, qgraupel):

    from __externals__ import c1_ice, c1_liq, c1_vap

    q_liq = qliquid + qrain
    q_solid = qice + qsnow + qgraupel
    return 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice


@gtscript.function
def vent_coeff(qden, density_factor, c1, c2, blin, mu):
    """
    Ventilation coefficient, Lin et al. (1983)
    """

    return c1 + c2 * exp((3 + 2 * mu + blin) / (mu + 3) / 2 * log(6 * qden)) * sqrt(
        density_factor
    ) / exp((1 + mu) / (mu + 3) * log(6 * qden))


@gtscript.function
def melting_function(
    tc,
    dq,
    qden,
    pxacw,
    pxacr,
    density,
    density_factor,
    lcpk,
    icpk,
    cvm,
    blin,
    mu,
    c1,
    c2,
    c3,
    c4,
):
    """
    Melting function, Lin et al. (1983)
    Fortran name is pmlt
    """
    return (c1 / (icpk * cvm) * tc / density - c2 * lcpk / icpk * dq) * exp(
        (1 + mu) / (mu + 3) * log(6 * qden)
    ) * vent_coeff(qden, density_factor, c3, c4, blin, mu) + constants.C_LIQ0 / (
        icpk * cvm
    ) * tc * (
        pxacw + pxacr
    )


@gtscript.function
def sublimation_function(
    t2, dq, qden, qsat, density, density_factor, cpk, cvm, c1, c2, c3, c4, c5, blin, mu
):
    """
    Sublimation or evaporation function, Lin et al. (1983)
    Fortran name is psub
    """

    return (
        c1
        * t2
        * dq
        * exp((1 + mu) / (mu + 3) * log(6 * qden))
        * vent_coeff(qden, density_factor, c2, c3, blin, mu)
        / (c4 * t2 + c5 * (cpk * cvm) ** 2 * qsat * density)
    )


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
    """
    Fortran name is mte
    """
    from __externals__ import c1_ice, c1_liq, c1_vap, c_air

    q_liq = qliquid + qrain
    q_solid = qice + qsnow + qgraupel
    q_cond = q_liq + q_solid
    con = 1.0 - (qvapor + q_cond)
    if moist_q:
        cvm = con + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
    else:
        cvm = 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
    return constants.RGRAV * cvm * c_air * temp * delp
