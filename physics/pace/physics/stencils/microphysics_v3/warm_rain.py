import math

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import PARALLEL, computation, exp, interval, log, sqrt

import pace.physics.stencils.microphysics_v3.physical_functions as physfun
import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ

from ..._config import MicroPhysicsConfig


@gtscript.function
def vent_coeff(qden, density_factor, c1, c2):
    """
    Ventilation coefficient, Lin et al. (1983)
    """
    from __externals__ import blin, mu

    return c1 + c2 * exp((3 + 2 * mu + blin) / (mu + 3) / 2 * log(6 * qden)) * sqrt(
        density_factor
    ) / exp((1 + mu) / (mu + 3) * log(6 * qden))


@gtscript.function
def psub(t2, dq, qden, qsat, density, density_factor, cpk, cvm):
    """
    Sublimation or evaporation function, Lin et al. (1983)
    """
    from __externals__ import c1, c2, c3, c4, c5, mu

    return (
        c1
        * t2
        * dq
        * exp((1 + mu) / (mu + 3) * log(6 * qden))
        * vent_coeff(qden, density_factor, c2, c3)
        / (c4 * t2 + c5 * (cpk * cvm) ** 2 * qsat * density)
    )


def evaporate_rain(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    delp: FloatField,
    temperature: FloatField,
    density: FloatField,
    density_factor: FloatField,
    reevap: FloatFieldIJ,
    h_var: FloatFieldIJ,
):
    """
    Rain evaporation to form water vapor, Lin et al. (1983)
    Fortran name is prevp
    """
    from __externals__ import (
        c1_ice,
        c1_liq,
        c1_vap,
        fac_revap,
        lv00,
        rhc_revap,
        t_wfr,
        timestep,
        use_rhc_revap,
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
        ) = physfun.calc_heat_cap_and_latent_heat_coeff(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temperature,
        )

        tin = (temperature * cvm - lv00 * qliquid) / (
            1.0 + (qvapor + qliquid) * c1_vap + qrain * c1_liq + q_solid * c1_ice
        )

        # calculate supersaturation and subgrid variability of water
        qpz = qvapor + qliquid
        qsat, dqdt = physfun.sat_spec_hum_water(tin, density)
        dqv = qsat - qvapor

        dqh = max(qliquid, h_var * max(qpz, constants.QCMIN))
        dqh = min(dqh, 0.2 * qpz)

        q_minus = qpz - dqh
        q_plus = qpz + dqh

        rh_tem = qpz / qsat

        if (
            (temperature > t_wfr)
            and (qrain > constants.QCMIN)
            and (dqv > 0.0)
            and (qsat > q_minus)
        ):
            if qsat > q_plus:
                dq = qsat - qpz
            else:
                dq = 0.25 * (qsat - q_minus) ** 2 / dqh
            qden = qrain * density
            t2 = tin * tin
            sink = psub(t2, dq, qden, qsat, density, density_factor, lcpk, cvm)
            sink = min(qrain, timestep * fac_revap * sink, dqv / (1.0 + lcpk * dqdt))
            if (use_rhc_revap is True) and (rh_tem >= rhc_revap):
                sink = 0

            reevap += sink * delp

            (
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                cvm,
                temperature,
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
                0.0,
                -sink,
                0.0,
                0.0,
                0.0,
                te,
            )


class WarmRain:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: pace.util.QuantityFactory,
        config: MicroPhysicsConfig,
        timestep: float,
        convert_mm_day: float,
    ):

        self._idx: GridIndexing = stencil_factory.grid_indexing
        self._fac_revap = 1.0
        if config.tau_revp > 1.0e6:
            self._fac_revap = 1.0 - math.exp(-timestep / config.tau_revp)

        self._evaporate_rain = stencil_factory.from_origin_domain(
            func=evaporate_rain,
            externals={
                "timestep": timestep,
                "mu": config.mur,
                "c1_vap": config.c1_vap,
                "c1_liq": config.c1_liq,
                "c1_ice": config.c1_ice,
                "t_wfr": config.t_wfr,
                "fac_revap": self._fac_revap,
                "use_rhc_revap": config.use_rhc_revap,
                "rhc_revap": config.rhc_revap,
                "d1_ice": config.d1_ice,
                "d1_vap": config.d1_vap,
                "li00": config.li00,
                "li20": config.li20,
                "lv00": config.lv00,
                "tice": config.tice,
                "c1": config.crevp_1,
                "c2": config.crevp_2,
                "c3": config.crevp_3,
                "c4": config.crevp_4,
                "c5": config.crevp_5,
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
        temperature: FloatField,
        delp: FloatField,
        delz: FloatField,
        density: FloatField,
        density_factor: FloatField,
        vterminal_water: FloatField,
        vterminal_rain: FloatField,
        cloud_condensation_nuclei: FloatField,
        reevap: FloatFieldIJ,
        h_var: FloatFieldIJ,
    ):
        pass
