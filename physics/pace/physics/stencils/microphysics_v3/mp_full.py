from gt4py.cartesian.gtscript import FORWARD, computation, interval

import pace.util

# from pace.dsl.dace.orchestration import orchestrate
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ
from pace.physics.stencils.microphysics_v3.sedimentation import Sedimentation
from pace.util import X_DIM, Y_DIM, Z_DIM

from ..._config import MicroPhysicsConfig


def add_fluxes_and_surface_tracers(
    prefluxw: FloatField,
    prefluxr: FloatField,
    prefluxi: FloatField,
    prefluxs: FloatField,
    prefluxg: FloatField,
    pfw: FloatField,
    pfr: FloatField,
    pfi: FloatField,
    pfs: FloatField,
    pfg: FloatField,
    water: FloatFieldIJ,
    rain: FloatFieldIJ,
    ice: FloatFieldIJ,
    snow: FloatFieldIJ,
    graupel: FloatFieldIJ,
    w1: FloatFieldIJ,
    r1: FloatFieldIJ,
    i1: FloatFieldIJ,
    s1: FloatFieldIJ,
    g1: FloatFieldIJ,
):
    from __externals__ import convt

    with computation(FORWARD):
        with interval(...):
            prefluxw += pfw * convt
            prefluxr += pfr * convt
            prefluxi += pfi * convt
            prefluxs += pfs * convt
            prefluxg += pfg * convt
        with interval(-1, None):
            water += w1 * convt
            rain += r1 * convt
            ice += i1 * convt
            snow += s1 * convt
            graupel += g1 * convt


class FullMicrophysics:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: pace.util.QuantityFactory,
        config: MicroPhysicsConfig,
        timestep: float,
        convert_mm_day: float,
    ):

        self._idx: GridIndexing = stencil_factory.grid_indexing

        def make_quantity():
            return quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="unknown")

        def make_quantity_2D():
            return quantity_factory.zeros([X_DIM, Y_DIM], units="unknown")

        self._fluxw = make_quantity()
        self._fluxr = make_quantity()
        self._fluxi = make_quantity()
        self._fluxs = make_quantity()
        self._fluxg = make_quantity()

        self._vtw = make_quantity()
        self._vtr = make_quantity()
        self._vti = make_quantity()
        self._vts = make_quantity()
        self._vtg = make_quantity()

        self._w1 = make_quantity_2D()
        self._r1 = make_quantity_2D()
        self._i1 = make_quantity_2D()
        self._s1 = make_quantity_2D()
        self._g1 = make_quantity_2D()

        self._reevap = make_quantity_2D()

        self._sedimentation = Sedimentation(
            stencil_factory, quantity_factory, config, timestep, convert_mm_day
        )

        self._add_fluxes_and_surface_tracers = stencil_factory.from_origin_domain(
            func=add_fluxes_and_surface_tracers,
            externals={
                "convt": convert_mm_day,
            },
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
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
        density: FloatField,
        density_factor: FloatField,
        preflux_water: FloatField,
        preflux_rain: FloatField,
        preflux_ice: FloatField,
        preflux_snow: FloatField,
        preflux_graupel: FloatField,
        column_energy_change: FloatFieldIJ,
        surface_water: FloatFieldIJ,
        surface_rain: FloatFieldIJ,
        surface_ice: FloatFieldIJ,
        surface_snow: FloatFieldIJ,
        surface_graupel: FloatFieldIJ,
    ):
        """
        Full Microphysics Loop
        """
        self._sedimentation(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            ua,
            va,
            wa,
            temperature,
            delp,
            delz,
            density,
            density_factor,
            self._fluxw,
            self._fluxr,
            self._fluxi,
            self._fluxs,
            self._fluxg,
            self._vtw,
            self._vtr,
            self._vti,
            self._vts,
            self._vtg,
            column_energy_change,
            self._w1,
            self._r1,
            self._i1,
            self._s1,
            self._g1,
        )

        self._add_fluxes_and_surface_tracers(
            preflux_water,
            preflux_rain,
            preflux_ice,
            preflux_snow,
            preflux_graupel,
            self._fluxw,
            self._fluxr,
            self._fluxi,
            self._fluxs,
            self._fluxg,
            surface_water,
            surface_rain,
            surface_ice,
            surface_snow,
            surface_graupel,
            self._w1,
            self._r1,
            self._i1,
            self._s1,
            self._g1,
        )

        pass
