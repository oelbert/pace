from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import __INLINED, PARALLEL, computation, interval

import pace.util.constants as constants
from pace.dsl.typing import FloatField


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

    q_liq = qrain + qliquid
    q_solid = qice + qsnow + qgraupel
    q_cond = q_liq + q_solid
    con = 1.0 - (qvapor + q_cond)
    if moist_q:
        cvm = con + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
    else:
        cvm = 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice
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
        if __INLINED(hydrostatic):
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
