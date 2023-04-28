import pace.dsl
import pace.util
import pace.util.constants as constants
from pace.dsl.typing import FloatField, FloatFieldIJ
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.terminal_fall import (
    prep_terminal_fall,
    update_energy_wind_heat_post_fall,
)
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class StartFall:
    def __init__(
        self,
        stencil_factory: pace.dsl.StencilFactory,
        config,
    ):
        self._idx = stencil_factory.grid_indexing
        self.config = config
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

    def __call__(
        self,
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
        tot_e_initial: FloatFieldIJ,
        no_fall: FloatFieldIJ,
    ):
        self._prep_terminal_fall(
            q_fall,
            delp,
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temperature,
            dm,
            tot_e_initial,
            no_fall,
        )


class EndFall:
    def __init__(
        self,
        stencil_factory: pace.dsl.StencilFactory,
        config,
    ):
        self._idx = stencil_factory.grid_indexing
        self.config = config
        self._update_energy_wind_heat_post_fall = stencil_factory.from_origin_domain(
            func=update_energy_wind_heat_post_fall,
            externals={
                "do_sedi_uv": config.do_sedi_uv,
                "do_sedi_w": config.do_sedi_w,
                "do_sedi_heat": config.do_sedi_heat,
                "cw": constants.C_ICE,
                "c1_ice": config.c1_ice,
                "c1_liq": config.c1_liq,
                "c1_vap": config.c1_vap,
                "c_air": config.c_air,
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
        ua: FloatField,
        va: FloatField,
        wa: FloatField,
        temperature: FloatField,
        delp: FloatField,
        delz: FloatField,
        flux: FloatField,
        dm: FloatField,
        v_terminal: FloatField,
        tmp_energy1: FloatFieldIJ,
        tmp_energy2: FloatFieldIJ,
        column_energy_change: FloatFieldIJ,
        no_fall: FloatFieldIJ,
    ):
        self._update_energy_wind_heat_post_fall(
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
            flux,
            dm,
            v_terminal,
            tmp_energy1,
            tmp_energy2,
            column_energy_change,
            no_fall,
        )


class TranslateStartFall(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "q_fall": {"serialname": "sf_qf", "mp3": True},
            "qvapor": {"serialname": "sf_qv", "mp3": True},
            "qliquid": {"serialname": "sf_ql", "mp3": True},
            "qrain": {"serialname": "sf_qr", "mp3": True},
            "qice": {"serialname": "sf_qi", "mp3": True},
            "qsnow": {"serialname": "sf_qs", "mp3": True},
            "qgraupel": {"serialname": "sf_qg", "mp3": True},
            "temperature": {"serialname": "sf_pt", "mp3": True},
            "delp": {"serialname": "sf_delp", "mp3": True},
            "dm": {"serialname": "sf_dm", "mp3": True},
            "tot_e_initial": {"serialname": "sf_e1", "mp3": True},
            "no_fall": {"serialname": "sf_nf", "mp3": True},
        }

        self.out_vars = {
            "dm": {"serialname": "sf_dm", "kend": namelist.npz, "mp3": True},
            "tot_e_initial": {"serialname": "sf_e1", "mp3": True},
            "no_fall": {"serialname": "sf_nf", "mp3": True},
        }

        self.stencil_factory = stencil_factory
        pconf = PhysicsConfig.from_namelist(namelist)
        pconf.hydrostatic = True
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = StartFall(
            self.stencil_factory,
            self.config,
        )

        compute_func(**inputs)

        return self.slice_output(inputs)


class TranslateEndFall(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "ef_qv", "mp3": True},
            "qliquid": {"serialname": "ef_ql", "mp3": True},
            "qrain": {"serialname": "ef_qr", "mp3": True},
            "qice": {"serialname": "ef_qi", "mp3": True},
            "qsnow": {"serialname": "ef_qs", "mp3": True},
            "qgraupel": {"serialname": "ef_qg", "mp3": True},
            "ua": {"serialname": "ef_ua", "mp3": True},
            "va": {"serialname": "ef_va", "mp3": True},
            "wa": {"serialname": "ef_wa", "mp3": True},
            "temperature": {"serialname": "ef_pt", "mp3": True},
            "delp": {"serialname": "ef_dp", "mp3": True},
            "delz": {"serialname": "ef_dz", "mp3": True},
            "flux": {"serialname": "ef_pfi", "mp3": True},
            "dm": {"serialname": "ef_dm", "mp3": True},
            "v_terminal": {"serialname": "ef_vt", "mp3": True},
            "column_energy_change": {"serialname": "ef_dte", "mp3": True},
            "tmp_energy1": {"serialname": "ef_ie", "mp3": True},
            "tmp_energy2": {"serialname": "ef_fe", "mp3": True},
            "no_fall": {"serialname": "ef_nf", "mp3": True},
        }

        self.out_vars = {
            "ua": {"serialname": "ef_ua", "kend": namelist.npz, "mp3": True},
            "va": {"serialname": "ef_va", "kend": namelist.npz, "mp3": True},
            "wa": {"serialname": "ef_wa", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "ef_pt", "kend": namelist.npz, "mp3": True},
            "tmp_energy1": {
                "serialname": "ef_ie",
                "mp3": True,
            },
            "column_energy_change": {"serialname": "ef_dte", "mp3": True},
        }

        self.stencil_factory = stencil_factory
        pconf = PhysicsConfig.from_namelist(namelist)
        pconf.hydrostatic = True
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = EndFall(
            self.stencil_factory,
            self.config,
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
