import physical_functions as physfun  # noqa
from gt4py.cartesian import gtscript  # noqa
from gt4py.cartesian.gtscript import __INLINED, BACKWARD, FORWARD, computation, interval
from moist_total_energy import calc_moist_total_energy

import pace.util
import pace.util.constants as constants

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ
from pace.fv3core.stencils.map_single import MapSingle
from pace.util import X_DIM, Y_DIM, Z_DIM

from ..._config import MicroPhysicsConfig


def prep_terminal_fall(
    q_fall: FloatField,
    delp: FloatField,
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    temperature: FloatField,
    dm: FloatField,
    tot_e_initial: FloatField,
    no_fall: FloatFieldIJ,
):
    from __externals__ import do_sedi_w

    with computation(BACKWARD), interval(...):
        if q_fall < constants.QFMIN:
            no_fall = 0.0
    with computation(FORWARD), interval(...):
        if no_fall > 0.0:
            if __INLINED(do_sedi_w):
                dm = delp * (1.0 + qvapor + qliquid + qrain + qice + qsnow + qgraupel)
            tot_e_initial = calc_moist_total_energy(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                temperature,
                delp,
                False,
            )


def fall_implicit(
    q_fall: FloatField,
    v_terminal: FloatField,
    delp: FloatField,
    z_edge: FloatField,
    no_fall: FloatFieldIJ,
    m1: FloatField,
    precip: FloatFieldIJ,
):
    from __externals__ import timestep

    with computation(FORWARD):
        with interval(...):
            if no_fall > 0.0:
                delz = z_edge - z_edge[0, 0, 1]
                dd = timestep * v_terminal
                q_fall = q_fall * delp

        with interval(0, 1):
            if no_fall > 0.0:
                qm = q_fall / (delz + dd)

        with interval(1, None):
            if no_fall > 0.0:
                qm = (q_fall + qm[0, 0, -1] * dd[0, 0, -1]) / (delz + dd)

        with interval(...):
            if no_fall > 0.0:
                qm = qm * delz

        with interval(0, 1):
            if no_fall > 0.0:
                m1 = q_fall - qm

        with interval(1, None):
            if no_fall > 0.0:
                m1 = m1[0, 0, -1] + q_fall - qm

        with interval(-1, None):
            if no_fall > 0.0:
                precip += m1

        with interval(...):
            if no_fall > 0.0:
                q_fall = qm / delp


class TerminalFall:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: pace.util.QuantityFactory,
        config: MicroPhysicsConfig,
        timestep: float,
    ):
        self._sedflag = config.sedflag
        if self._sedflag not in [1, 2, 3, 4]:
            raise ValueError(
                f"sedimentation flag {self._sedflag} must be 1, 2, 3, or 4"
            )

        if self._sedflag == 2:
            raise NotImplementedError(
                f"sedimentation flag {self._sedflag}: "
                "Explicit Fall has not been implemented"
            )

        self._idx: GridIndexing = stencil_factory.grid_indexing

        self._timestep = timestep

        # allocate internal storages
        self._no_fall = quantity_factory.ones(dims=[X_DIM, Y_DIM], units="unknown")
        self._dm = quantity_factory.ones(dims=[X_DIM, Y_DIM, Z_DIM], units="Pa")
        self._prep_terminal_fall = stencil_factory.from_origin_domain(
            func=prep_terminal_fall,
            externals={
                "do_sedi_w": config.do_sedi_w,
                "c1_ice": config.c1_ice,
                "c1_liq": config.c1_liq,
                "c1_vap": config.c1_vap,
                "c_air": config.c_air,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        if (self._sedflag == 3) or (self._sedflag == 4):
            self._lagrangian = MapSingle(
                stencil_factory,
                quantity_factory,
                kord=9,
                mode=0,
                dims=[X_DIM, Y_DIM, Z_DIM],
            )

        pass

    def __call__(
        self,
        qvapor: FloatField,
        qliquid: FloatField,
        qrain: FloatField,
        qice: FloatField,
        qsnow: FloatField,
        qgraupel: FloatField,
        ua: FloatField,
        va: FloatField,
        wa: FloatField,
        temperature: FloatField,
        delp: FloatField,
        delz: FloatField,
        v_terminal: FloatField,
        z_surface: FloatFieldIJ,
        z_edge: FloatField,
        z_terminal: FloatField,
        column_energy_change: FloatFieldIJ,
    ):
        pass
