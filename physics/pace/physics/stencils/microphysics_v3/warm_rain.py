import math

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import __INLINED, FORWARD, PARALLEL, computation, exp, interval, log, sqrt

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


def accrete_rain(
    qliquid: FloatField,
    qrain: FloatField,
    temperature: FloatField,
    density: FloatField,
    density_factor: FloatField,
    vterminal_water: FloatField,
    vterminal_rain: FloatField,

):
    """
    Rain accretion with cloud water, Lin et al. (1983)
    Fortran name is pracw
    """
    from __externals__ import t_wfr, timestep, do_new_acc_water, cracw, acco0_4, acco1_4, acco2_4, acc8, acc9, blinr, mur

    with computation(PARALLEL), interval(...):
        if (temperature > t_wfr) and (qrain > constants.QCMIN) and (qliquid > constants.QCMIN):
            qden = qrain * density

            if __INLINED(do_new_acc_water is True):
                sink = timestep * physfun.accretion_3d(qliquid, qrain, vterminal_rain, vterminal_water, density, cracw, acco0_4, acco1_4, acco2_4, acc8, acc9)
            else:
                sink = timestep * physfun.accretion_2d(qden, cracw, density_factor, blinr, mur)
                sink = sink / (1. + sink) * qliquid

            qliquid -= sink
            qrain += sink


def autoconvert_water_rain(
    qliquid: FloatField,
    qrain: FloatField,
    cloud_condensation_nuclei: FloatField,
    temperature: FloatField,
    density: FloatField,
    h_var: FloatFieldIJ,
):
    """
    Cloud water to rain autoconversion, Hong et al. (2004)
    Fortran name is praut
    """
    from __externals__ import timestep, irain_f, z_slope_liq, t_wfr, do_psd_water_num, pcaw, pcbw, muw, fac_rc, cpaut
    
    with computation(FORWARD):
        with interval(0, 1):
            if __INLINED((irain_f == 0) and (z_slope_liq is True)):
                # linear_prof
                dl = 0.
        with interval(1, None):
            if __INLINED((irain_f == 0) and (z_slope_liq is True)):
                dq = 0.5 * (qliquid - qliquid[0, 0, -1])
        with interval(1,-1):
            if __INLINED((irain_f == 0) and (z_slope_liq is True)):
                # Use twice the strength of the
                # positive definiteness limiter (lin et al 1994)
                dl = 0.5 * min(abs(dq + dq[0, 0, +1]), 0.5 * qliquid[0, 0, 0])
                if dq * dq[0, 0, +1] <= 0.0:
                    if dq > 0.0:  # Local maximum
                        dl = min(dl, min(dq, -dq[0, 0, +1]))
                    else: # Local minimum
                        dl = 0.0
        with interval(-1, None):
            if __INLINED((irain_f == 0) and (z_slope_liq is True)):
                dl = 0.
    with computation(PARALLEL), interval(...):
        if __INLINED(irain_f == 0):
            if __INLINED(z_slope_liq is True):
                # Impose a presumed background horizontal variability that is
                # proportional to the value itself
                dl = max(dl, 0., h_var * qliquid)
            else:
                dl = max(0., h_var * qliquid)
            
            #rest of praut
            if (temperature > t_wfr) and (qliquid > constants.QCMIN):
                if __INLINED(do_psd_water_num):
                    cloud_condensation_nuclei = physfun.calc_particle_concentration(qliquid, density, pcaw, pcbw, muw)
                    cloud_condensation_nuclei /= density
                qc = fac_rc * cloud_condensation_nuclei
                dl = min(max(constants.QCMIN, dl), 0.5*qliquid)
                dq = 0.5 * (qliquid + dl - qc)

                if dq > 0.:
                    c_praut = cpaut * exp((-1./3.) * log(cloud_condensation_nuclei * constants.RHO_W))
                    sink = min (1., dq / dl) * timestep * c_praut * density * exp((7./3.) * log (qliquid))
                    sink = min(qliquid, sink)

                    qliquid -= sink
                    qrain += sink
        else: # if irain_f == 1:
            if (temperature > t_wfr) and (qliquid > constants.QCMIN):
                if __INLINED(do_psd_water_num):
                    cloud_condensation_nuclei = physfun.calc_particle_concentration(qliquid, density, pcaw, pcbw, muw)
                    cloud_condensation_nuclei /= density
                
                qc = fac_rc * cloud_condensation_nuclei
                dq = qliquid - qc

                if dq > 0:
                    c_praut = cpaut * exp((-1./3.) * log(cloud_condensation_nuclei * constants.RHO_W))
                    sink = min(dq, timestep * c_praut * density * exp((7./3.) * log(qliquid)))
                    sink = min(sink, qliquid)

                    qliquid -= sink
                    qrain += sink


class WarmRain:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        config: MicroPhysicsConfig,
        timestep: float,
    ):

        self._idx: GridIndexing = stencil_factory.grid_indexing
        self._fac_revap = 1.0
        if config.tau_revp > 1.0e6:
            self._fac_revap = 1.0 - math.exp(-timestep / config.tau_revp)

        fac_rc = (4. / 3.) * constants.PI * constants.RHO_W * config.rthresh ** 3
        aone = 2. / 9. * (3. / 4.) ** (4. / 3.) / constants.PI ** (1. / 3.)
        cpaut = config.c_paut * aone * constants.GRAV / constants.VISD

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

        self._accrete_rain = stencil_factory.from_origin_domain(
            func=accrete_rain,
            externals={
                "timestep": timestep,
                "t_wfr": config.t_wfr,
                "do_new_acc_water": config.do_new_acc_water,
                "cracw": config.cracw,
                "acco0_4": config.acco[0,4],
                "acco1_4": config.acco[1,4],
                "acco2_4": config.acco[2,4],
                "acc9": config.acc[8],
                "acc9": config.acc[9],
                "blinr": config.blinr,
                "mur": config.mur,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._autoconvert_water_rain = stencil_factory.from_origin_domain(
            func=autoconvert_water_rain,
            externals={
                "timestep": timestep,
                "t_wfr": config.t_wfr,
                "irain_f": config.irain_f,
                "z_slope_liq": config.z_slope_liq,
                "do_psd_water_num": config.do_psd_water_num,
                "pcaw": config.pcaw,
                "pcbw": config.pcbw,
                "muw": config.muw,
                "fac_rc": fac_rc,
                "cpaut": cpaut,
                
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
        density: FloatField,
        density_factor: FloatField,
        vterminal_water: FloatField,
        vterminal_rain: FloatField,
        cloud_condensation_nuclei: FloatField,
        reevap: FloatFieldIJ,
        h_var: FloatFieldIJ,
    ):
        """
        Warm rain cloud microphysics
        Args:
            qvapor (inout):
            qliquid (inout):
            qrain (inout):
            qice (inout):
            qsnow (inout):
            qgraupel (inout):
            temperature (inout):
            delp (in):
            density (in):
            density_factor (in):
            vterminal_water (in):
            vterminal_rain (in):
            cloud_condensation_nuclei (inout):
            reevap (inout):
            h_var (in):
        """
        self._evaporate_rain(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            delp,
            temperature,
            density,
            density_factor,
            reevap,
            h_var,
        )

        self._accrete_rain(
            qliquid,
            qrain,
            temperature,
            density,
            density_factor,
            vterminal_water,
            vterminal_rain,
        )

        self._autoconvert_water_rain(
            qliquid,
            qrain,
            cloud_condensation_nuclei,
            temperature,
            density,
            h_var,
        )
