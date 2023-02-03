from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import exp, log

import pace.util.constants as constants


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
