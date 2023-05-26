import math

from gt4py.cartesian.gtscript import (  # noqa
    __INLINED,
    FORWARD,
    PARALLEL,
    computation,
    interval,
)

import pace.dsl
import pace.util
import pace.util.constants as constants
from pace.dsl.typing import FloatField, FloatFieldIJ
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.warm_rain import (  # noqa
    accrete_rain,
    autoconvert_water_rain,
    evaporate_rain,
)
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class RainFunction:
    def __init__(
        self,
        stencil_factory,
        config,
        timestep: float,
    ):

        self._idx = stencil_factory.grid_indexing
        self._fac_revap = 1.0
        if config.tau_revp > 1.0e-6:
            self._fac_revap = 1.0 - math.exp(-timestep / config.tau_revp)

        fac_rc = (4.0 / 3.0) * constants.PI * constants.RHO_W * config.rthresh ** 3
        aone = 2.0 / 9.0 * (3.0 / 4.0) ** (4.0 / 3.0) / constants.PI ** (1.0 / 3.0)
        cpaut = config.c_paut * aone * constants.GRAV / constants.VISD

        self._evaporate_rain = stencil_factory.from_origin_domain(
            func=evaporate_rain,
            externals={
                "timestep": timestep,
                "mur": config.mur,
                "blinr": config.mur,
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
                "tice": constants.TICE0,
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
                "acco0_4": config.acco[0][4],
                "acco1_4": config.acco[1][4],
                "acco2_4": config.acco[2][4],
                "acc8": config.acc[8],
                "acc9": config.acc[9],
                "blinr": config.blinr,
                "mur": config.mur,
                "vdiffflag": config.vdiffflag,
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
        # self._evaporate_rain(
        #     qvapor,
        #     qliquid,
        #     qrain,
        #     qice,
        #     qsnow,
        #     qgraupel,
        #     delp,
        #     temperature,
        #     density,
        #     density_factor,
        #     reevap,
        #     h_var,
        # )

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


class TranslateWRainSubFunc(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "ws_qv", "mp3": True},
            "qliquid": {"serialname": "ws_ql", "mp3": True},
            "qrain": {"serialname": "ws_qr", "mp3": True},
            "qice": {"serialname": "ws_qi", "mp3": True},
            "qsnow": {"serialname": "ws_qs", "mp3": True},
            "qgraupel": {"serialname": "ws_qg", "mp3": True},
            "temperature": {"serialname": "ws_pt", "mp3": True},
            "density": {"serialname": "ws_den", "mp3": True},
            "density_factor": {"serialname": "ws_denfac", "mp3": True},
            "vterminal_water": {"serialname": "ws_vtw", "mp3": True},
            "vterminal_rain": {"serialname": "ws_vtr", "mp3": True},
            "delp": {"serialname": "ws_delp", "mp3": True},
            "h_var": {"serialname": "ws_h_var", "mp3": True},
            "cloud_condensation_nuclei": {"serialname": "ws_ccn", "mp3": True},
            "reevap": {"serialname": "ws_reevap", "mp3": True},
        }

        self.in_vars["parameters"] = [
            "dt",
        ]

        self.out_vars = {
            "qvapor": {"serialname": "ws_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "ws_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "ws_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "ws_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "ws_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "ws_qg", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "ws_pt", "kend": namelist.npz, "mp3": True},
            "cloud_condensation_nuclei": {
                "serialname": "ws_ccn",
                "kend": namelist.npz,
                "mp3": True,
            },
            "reevap": {"serialname": "ws_reevap", "mp3": True},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = RainFunction(
            self.stencil_factory,
            self.config,
            timestep=inputs.pop("dt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
