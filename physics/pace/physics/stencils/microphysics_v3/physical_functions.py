from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import exp, log, sqrt

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
def accretion_2d(
    qden,
    denfac,
    c,
    blin,
    mu
):
    return denfac * c * exp((2 + mu + blin) / (mu + 3) * log(6*qden))


@gtscript.function
def accretion_3d(
    q1,
    q2,
    v1,
    v2,
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

    t1 = exp(1. / (acc1 + 3) * log(6 * q1 * density))
    t2 = exp(1. / (acc2 + 3) * log(6 * q2 * density))

    if vdiffflag == 1:
        vdiff = abs(v1-v2)
    elif vdiffflag == 2:
        sqrt((1.20 * v1 - 0.95 * v2) ** 2. + 0.08 * v1 * v2)
    else: # vdiffflag == 3:
        vdiff = sqrt((1.00 * v1 - 1.00 * v2) ** 2. + 0.04 * v1 * v2)

    accrete = c * vdiff / density
    tmp = acco1 * exp((6+acc1 - 1) * log(t1)) * exp((acc2 + 1 - 1) * log(t2))
    tmp += acco2 * exp((6+acc1 - 2) * log(t1)) * exp((acc2 + 2 - 1) * log(t2))
    tmp += acco3 * exp((6+acc1 - 3) * log(t1)) * exp((acc2 + 3 - 1) * log(t2))

    return accrete * tmp


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
        c1_vap,
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
    cvm = 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
    te = cvm + temp + lv00 * qvapor - li00 * q_solid
    lcpk = (lv00 + d1_vap * temp) / cvm
    icpk = (li00 + d1_ice * temp) / cvm
    tcpk = (li20 + (d1_vap + d1_ice) * temp) / cvm
    tcp3 = lcpk + icpk * min(1.0, basic.dim(tice, temp) / (tice - t_wfr))

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
    cvm = 1.0 + qvapor * c1_vap + q_l * c1_liq + q_solid * c1_ice

    tk = (te - lv00 + qvapor + li00 * (qice + qsnow + qgraupel)) / cvm
    lcpk = (lv00 + d1_vap * tk) / cvm
    icpk = (li00 + d1_ice * tk) / cvm
    tcpk = (li20 + (d1_vap + d1_ice) * tk) / cvm
    tcp3 = lcpk + icpk * min(1.0, basic.dim(tice, tk) / (tice - t_wfr))

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
            constants.DC_VAP * log(temp / constants.TICE)
            + constants.LV0 * (temp - constants.TICE) / (temp * constants.TICE)
        )
        / constants.RVGAS
    )


@gtscript.function
def table2(temp):
    if temp < constants.TICE:
        # Over ice between -160 degrees Celsius and 0 degrees Celsius
        return constants.E00 * exp(
            (
                constants.D2ICE * log(temp / constants.TICE)
                + constants.LI2 * (temp - constants.TICE) / (temp * constants.TICE)
            )
            / constants.RVGAS
        )

    else:
        # Over water between 0 degrees Celsius and 102 degrees Celsius
        return table0(temp)


@gtscript.function
def sat_spec_hum_water(temp, density):
    """
    qs_core with table 0 in microphysics
    compute the saturated specific humidity, core function
    """
    q = table0(temp) / (constants.RVGAS * temp * density)
    dqdt = q * (constants.DC_VAP + constants.LV0 / temp) / (constants.RVGAS * temp)
    return q, dqdt


@gtscript.function
def sat_spec_hum_water_ice(temp, density):
    if temp > constants.TICE + 102.0:
        temp = constants.TICE + 102.0
    q = table2(temp) / (constants.RVGAS * temp * density)
    if temp < constants.TICE:
        dqdt = q * (constants.D2ICE + constants.LI2 / temp) / (constants.RVGAS * temp)
    else:
        dqdt = q * (constants.DC_VAP + constants.LV0 / temp) / (constants.RVGAS * temp)
    return q, dqdt


@gtscript.function
def moist_heat_capacity(qvapor, qliquid, qrain, qice, qsnow, qgraupel):

    from __externals__ import c1_ice, c1_liq, c1_vap

    q_liq = qliquid + qrain
    q_solid = qice + qsnow + qgraupel
    return 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
