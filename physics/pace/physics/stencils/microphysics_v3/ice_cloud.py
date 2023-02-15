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
)

import pace.fv3core.stencils.basic_operations as basic
import pace.physics.stencils.microphysics_v3.physical_functions as physfun
import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ
from pace.util import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM

from ..._config import MicroPhysicsConfig


@gtscript.function
def melt_cloud_ice(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temperature,
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
    tc = temperature - tice_mlt
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
        temperature,
        cvm,
        lcpk,
        icpk,
        tcpk,
        tcp3,
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
def melt_snow(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    temperature,
    density,
    density_factor,
    vterminal_w,
    vterminal_r,
    vterminal_s,
    cvm,
    te,
    lcpk,
    icpk,
    tcpk,
    tcp3,
):
    """
    Snow melting (includes snow accretion with cloud water and rain)
    to form cloud water and rain, Lin et al. (1983)
    Fortran name is psmlt
    """

    from __externals__ import timestep, do_new_acc_water, csacw, csacr, acco_0_0, acco_1_0, acco_2_0, acco_0_1, acco_1_1, acco_2_1, acco_0_6, acco_1_6, acco_2_6, acc0, acc1, acc2, acc3, acc12, acc13, blins, mus

    tc = temperature - constants.TICE0

    if (tc >= 0) and (qsnow > constants.QCMIN):
        psacw = 0.
        qden = qsnow * density
        if qliquid > constants.QCMIN:
            if __INLINED(do_new_acc_water is True):
                psacw = physfun.accretion_3d(qliquid, qsnow, vterminal_s, vterminal_w, density, csacw, acco_0_6, acco_1_6, acco_2_6, acc12, acc13)
            else:
                factor = physfun.accretion_2d(qden, density_factor, csacw, blins, mus)
                psacw = factor / (1. + timestep * factor) * qliquid

        psacr = 0.
        pracs = 0.
        if qrain > constants.QCMIN:
            psacr = min(qrain/timestep, physfun.accretion_3d(
                qrain, qsnow, vterminal_s, vterminal_r, density, csacr, acco_0_1, acco_1_1, acco_2_1, acc2, acc3
            ))
            pracs = physfun.accretion_3d(
                qsnow, qrain, vterminal_r, vterminal_s, density, csacr, acco_0_0, acco_1_0, acco_2_0, acc0, acc1
            )

        tin = temperature
        qsi, dqdt = physfun.sat_spec_hum_water_ice(tin, density)
        dq = qsi - qvapor



    pass


def ice_cloud(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    temperature: FloatField,
    density: FloatField,
    h_var: FloatFieldIJ,
):
    from __externals__ import z_slope_ice

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

        (
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temperature,
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
            temperature,
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
    
    with computation(FORWARD):
        with interval(0, 1):
            if __INLINED(z_slope_ice is True):
                # linear_prof
                di = 0.
        with interval(1, None):
            if __INLINED(z_slope_ice is True):
                dq = 0.5 * (qice - qice[0, 0, -1])
        with interval(1,-1):
            if __INLINED(z_slope_ice is True):
                # Use twice the strength of the
                # positive definiteness limiter (lin et al 1994)
                di = 0.5 * min(abs(dq + dq[0, 0, +1]), 0.5 * qice[0, 0, 0])
                if dq * dq[0, 0, +1] <= 0.0:
                    if dq > 0.0:  # Local maximum
                        di = min(di, min(dq, -dq[0, 0, +1]))
                    else: # Local minimum
                        di = 0.0
        with interval(-1, None):
            if __INLINED(z_slope_ice is True):
                di = 0.
    with computation(PARALLEL), interval(...):
        if __INLINED(z_slope_ice is True):
            # Impose a presumed background horizontal variability that is
            # proportional to the value itself
            di = max(di, 0., h_var * qice)
        else:
            di = max(0., h_var * qice)

    pass



class IceCloud:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: pace.util.QuantityFactory,
        config: MicroPhysicsConfig,
        timestep: float,
    ):
        self._idx: GridIndexing = stencil_factory.grid_indexing
        pass

    def __call__(
        self,
        qvapor: FloatField,
        qliquid: FloatField,
        qrain: FloatField,
        qice: FloatField,
        qsnow: FloatField,
        qgraupel: FloatField,
        temperature: FloatField,
        density: FloatField,
        density_factor: FloatField,
        vterminal_water: FloatField,
        vterminal_rain: FloatField,
        vterminal_ice: FloatField,
        vterminal_snow: FloatField,
        vterminal_graupel: FloatField,
        h_var: FloatFieldIJ,
    ):
        pass