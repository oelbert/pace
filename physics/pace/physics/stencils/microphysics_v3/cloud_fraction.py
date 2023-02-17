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
import pace.physics.stencils.microphysics_v3.physical_functions as physfun
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ

from ..._config import MicroPhysicsConfig


@gtscript.function
def cloud_scheme_1(
    qpz,
    qstar,
    q_cond,
    qa,
    pz,
    rh,
    h_var,
):
    from __externals__ import cld_min, do_cld_adj, f_dq_m, f_dq_p, icloud_f, rh_thres

    if (rh > rh_thres) and (qpz > constants.QCMIN):
        dq = h_var * qpz
        if __INLINED(do_cld_adj):
            q_plus = qpz + dq * f_dq_p * min(
                1.0, max(0.0, (pz - 200.0e2) / (1000.0e2 - 200.0e2))
            )
        else:
            q_plus = qpz + dq * f_dq_p
        q_minus = qpz - dq * f_dq_m
        if __INLINED(icloud_f == 2):
            if qstar < qpz:
                qa = 1.0
            else:
                qa = 0.0
        elif __INLINED(icloud_f == 3):
            if qstar < qpz:
                qa = 1.0
            else:
                if qstar < q_plus:
                    qa = (q_plus - qstar) / (dq * f_dq_p)
                else:
                    qa = 0.0
                if q_cond > constants.QCMIN:
                    qa = max(cld_min, qa)
                qa = min(1.0, qa)
        else:
            if qstar < q_minus:
                qa = 1.0
            else:
                if qstar < q_plus:
                    if __INLINED(icloud_f == 0):
                        qa = (q_plus - qstar) / (dq * f_dq_p + dq * f_dq_m)
                    else:  # icloud_f == 1:
                        qa = (q_plus - qstar) / (
                            (dq * f_dq_p + dq * f_dq_m) * (1.0 - q_cond)
                        )
                else:
                    qa = 0.0
                if q_cond > constants.QCMIN:
                    qa = max(cld_min, qa)
                qa = min(1.0, qa)
    else:
        qa = 0.0

    return qa


def cloud_fraction(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    qa: FloatField,
    temperature: FloatField,
    density: FloatField,
    pz: FloatField,
    te: FloatField,
    h_var: FloatFieldIJ,
):

    from __externals__ import (
        c1_ice,
        c1_liq,
        c1_vap,
        cfflag,
        li00,
        lv00,
        rad_graupel,
        rad_rain,
        rad_snow,
        t_wfr,
    )

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
        (qvapor, qliquid, qrain, qice, qsnow, qgraupel, temperature)

        # Combine water species
        ice = q_solid
        q_solid = qice
        if __INLINED(rad_snow is True):
            q_solid += qsnow
            if __INLINED(rad_graupel is True):
                q_solid += qgraupel

        liq = q_liq
        q_liq = qliquid
        if __INLINED(rad_rain is True):
            q_liq += qrain

        q_cond = q_liq + q_solid
        qpz = qvapor + q_cond

        # use the "liquid - frozen water temperature" (tin)
        # to compute saturated specific humidity
        ice = ice - q_solid
        liq = liq - q_liq
        tin = (te - lv00 * qpz + li00 * ice) / (
            1.0 + qpz * c1_vap + liq * c1_liq + ice * c1_ice
        )

        if tin <= t_wfr:
            qstar = physfun.sat_spec_hum_water_ice(tin, density)
        elif tin >= constants.TICE0:
            qstar = physfun.sat_spec_hum_water(tin, density)
        else:
            qsi = physfun.sat_spec_hum_water_ice(tin, density)
            qsw = physfun.sat_spec_hum_water(tin, density)
            if q_cond > constants.QCMIN:
                rqi = q_solid / q_cond
            else:
                (constants.TICE0 - tin) / (constants.TICE0 - t_wfr)
            qstar = rqi * qsi + (1.0 - rqi) * qsw

        # Cloud schemes
        rh = qpz / qstar

        if __INLINED(cfflag == 1):
            qa = cloud_scheme_1(
                qpz,
                qstar,
                q_cond,
                qa,
                pz,
                rh,
                h_var,
            )
    pass


class CloudFraction:
    def __init__(self):
        pass

    def __call__(
        self,
        qvapor: FloatField,
        qliquid: FloatField,
        qrain: FloatField,
        qice: FloatField,
        qsnow: FloatField,
        qgraupel: FloatField,
        qa: FloatField,
        temperature: FloatField,
        density: FloatField,
        pz: FloatField,
        h_var: FloatFieldIJ,
        gsize: FloatFieldIJ,
    ):
        pass
