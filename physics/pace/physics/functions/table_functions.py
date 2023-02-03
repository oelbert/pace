from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import exp, log

import pace.util.constants as constants


@gtscript.stencil
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


@gtscript.stencil
def sat_spec_hum_water(temp, density):
    """
    qs_core with table 0 in microphysics
    compute the saturated specific humidity, core function
    """
    q = table0(temp) / (constants.RVGAS * temp * density)
    dqdt = q * (constants.DC_VAP + constants.LV0 / temp) / (constants.RVGAS * temp)
    return q, dqdt
