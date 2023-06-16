import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.microphysics_v3 import (
    calc_sedimentation_energy_loss,
    calculate_total_energy_change_and_convert_temp,
    convert_mass_mixing_to_specific_ratios_and_update_temperatures,
    moist_total_energy_and_water,
    total_energy_check,
    update_temperature_pre_delp_q,
)
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class FinalCalcs:
    def __init__(
        self,
        stencil_factory,
        config,
        consv_te,
    ):
        self.config = config
        self._idx = stencil_factory.grid_indexing
        self.consv_te = consv_te
        self.do_qa = config.do_qa
        self.hydrostatic = config.hydrostatic
        self.consv_checker = config.consv_checker
        self.do_inline_mp = config.do_inline_mp
        self.fix_negative = config.fix_negative
        self._te_err = config.te_err
        self._tw_err = config.tw_err

        if self.config.consv_checker:
            self._moist_total_energy_and_water_mq = stencil_factory.from_origin_domain(
                func=moist_total_energy_and_water,
                externals={
                    "hydrostatic": self.config.hydrostatic,
                    "moist_q": True,
                    "c1_vap": self.config.c1_vap,
                    "c1_liq": self.config.c1_liq,
                    "c1_ice": self.config.c1_ice,
                    "lv00": self.config.lv00,
                    "li00": self.config.li00,
                    "c_air": self.config.c_air,
                    "timestep": self.config.dt_full,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

            self._moist_total_energy_and_water = stencil_factory.from_origin_domain(
                func=moist_total_energy_and_water,
                externals={
                    "hydrostatic": self.config.hydrostatic,
                    "moist_q": False,
                    "c1_vap": self.config.c1_vap,
                    "c1_liq": self.config.c1_liq,
                    "c1_ice": self.config.c1_ice,
                    "lv00": self.config.lv00,
                    "li00": self.config.li00,
                    "c_air": self.config.c_air,
                    "timestep": self.config.dt_full,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

            self._calc_sedimentation_energy_loss = stencil_factory.from_origin_domain(
                func=calc_sedimentation_energy_loss,
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        if (self.do_sedi_uv) or (self.do_sedi_w):
            self._update_temperature_pre_delp_q = stencil_factory.from_origin_domain(
                func=update_temperature_pre_delp_q,
                externals={
                    "do_sedi_uv": self.do_sedi_uv,
                    "do_sedi_w": self.do_sedi_w,
                    "c_air": config.c_air,
                    "c1_vap": config.c1_vap,
                    "c1_liq": config.c1_liq,
                    "c1_ice": config.c1_ice,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        self._convert_mass_mixing_to_specific_ratios_and_update_temperatures = (
            stencil_factory.from_origin_domain(
                func=convert_mass_mixing_to_specific_ratios_and_update_temperatures,
                externals={
                    "do_inline_mp": config.do_inline_mp,
                    "c_air": config.c_air,
                    "c1_vap": config.c1_vap,
                    "c1_liq": config.c1_liq,
                    "c1_ice": config.c1_ice,
                    "do_sedi_uv": config.do_sedi_uv,
                    "do_sedi_w": config.do_sedi_w,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )
        )

        self._calculate_total_energy_change_and_convert_temp = (
            stencil_factory.from_origin_domain(
                func=calculate_total_energy_change_and_convert_temp,
                externals={
                    "consv_te": consv_te,
                    "hydrostatic": self.hydrostatic,
                    "c_air": config.c_air,
                    "do_inline_mp": config.do_inline_mp,
                    "cp_heating": config.cp_heating,
                    "c1_vap": config.c1_vap,
                    "c1_liq": config.c1_liq,
                    "c1_ice": config.c1_ice,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )
        )

    def __call__(
        self,
        qvapor,
        qliquid,
        qrain,
        qice,
        qsnow,
        qgraupel,
        delp,
        delz,
        pt,
        ua,
        va,
        wa,
        total_energy,
        qvapor0,
        qliquid0,
        qrain0,
        qice0,
        qsnow0,
        qgraupel0,
        dp0,
        pt0,
        u0,
        v0,
        w0,
        tzuv,
        tzw,
        column_energy_change,
        adj_vmr,
        gsize,
        column_vapor,
        column_water,
        column_rain,
        column_ice,
        column_snow,
        column_graupel,
        total_energy_wet_begin,
        total_water_wet_begin,
        total_energy_bot_wet_begin,
        total_water_bot_wet_begin,
        total_energy_dry_begin,
        total_water_dry_begin,
        total_energy_bot_dry_begin,
        total_water_bot_dry_begin,
        total_energy_dry_end,
        total_energy_wet_end,
        total_water_dry_end,
        total_water_wet_end,
        total_energy_bot_dry_end,
        total_energy_bot_wet_end,
        total_water_bot_dry_end,
        total_water_bot_wet_end,
        column_energy_loss,
    ):
        if (self.do_sedi_uv) or (self.do_sedi_w):
            self._update_temperature_pre_delp_q(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                ua,
                va,
                wa,
                u0,
                v0,
                w0,
                pt,
                tzuv,
                tzw,
            )
        if self.consv_checker:
            self._moist_total_energy_and_water(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                pt,
                ua,
                va,
                wa,
                delp,
                gsize,
                column_energy_change,
                column_vapor,
                column_water,
                column_rain,
                column_ice,
                column_snow,
                column_graupel,
                0.0,
                0.0,
                total_energy_dry_end,
                total_water_dry_end,
                total_energy_bot_dry_end,
                total_water_bot_dry_end,
            )

            self._calc_sedimentation_energy_loss(
                column_energy_loss, column_energy_change, gsize
            )

        self._convert_mass_mixing_to_specific_ratios_and_update_temperatures(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            ua,
            va,
            wa,
            delp,
            qvapor0,
            qliquid0,
            qrain0,
            qice0,
            qsnow0,
            qgraupel0,
            u0,
            v0,
            w0,
            dp0,
            pt,
            tzuv,
            tzw,
            adj_vmr,
        )

        if self.consv_checker:
            self._moist_total_energy_and_water_mq(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                pt,
                ua,
                va,
                wa,
                delp,
                gsize,
                column_energy_change,
                column_vapor,
                column_water,
                column_rain,
                column_ice,
                column_snow,
                column_graupel,
                0.0,
                0.0,
                total_energy_wet_end,
                total_water_wet_end,
                total_energy_bot_wet_end,
                total_water_bot_wet_end,
            )

        self._calculate_total_energy_change_and_convert_temp(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            total_energy,
            pt,
            pt0,
            delp,
            delz,
        )

        if self.consv_checker:
            total_energy_check(
                total_energy_dry_end,
                total_energy_wet_end,
                total_water_dry_end,
                total_water_wet_end,
                total_energy_dry_begin,
                total_energy_wet_begin,
                total_water_dry_begin,
                total_water_wet_begin,
                total_energy_bot_dry_end,
                total_energy_bot_wet_end,
                total_water_bot_dry_end,
                total_water_bot_wet_end,
                total_energy_bot_dry_begin,
                total_energy_bot_wet_begin,
                total_water_bot_dry_begin,
                total_water_bot_wet_begin,
                self._idx.isc,
                self._idx.iec,
                self._idx.jsc,
                self._idx.jec,
                self._te_err,
                self._tw_err,
            )


class TranslateFinalCalculations(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "fin_qv", "mp3": True},
            "qliquid": {"serialname": "fin_ql", "mp3": True},
            "qrain": {"serialname": "fin_qr", "mp3": True},
            "qice": {"serialname": "fin_qi", "mp3": True},
            "qsnow": {"serialname": "fin_qs", "mp3": True},
            "qgraupel": {"serialname": "fin_qg", "mp3": True},
            "delp": {"serialname": "fin_delp", "mp3": True},
            "delz": {"serialname": "fin_delz", "mp3": True},
            "pt": {"serialname": "fin_pt", "mp3": True},
            "ua": {"serialname": "fin_ua", "mp3": True},
            "va": {"serialname": "fin_va", "mp3": True},
            "wa": {"serialname": "fin_wa", "mp3": True},
            "qvapor0": {"serialname": "fin_qv0", "mp3": True},
            "qliquid0": {"serialname": "fin_ql0", "mp3": True},
            "qrain0": {"serialname": "fin_qr0", "mp3": True},
            "qice0": {"serialname": "fin_qi0", "mp3": True},
            "qsnow0": {"serialname": "fin_qs0", "mp3": True},
            "qgraupel0": {"serialname": "fin_qg0", "mp3": True},
            "dp0": {"serialname": "fin_dp0", "mp3": True},
            "pt0": {"serialname": "fin_pt0", "mp3": True},
            "u0": {"serialname": "fin_u0", "mp3": True},
            "v0": {"serialname": "fin_v0", "mp3": True},
            "w0": {"serialname": "fin_w0", "mp3": True},
            "column_energy_change": {"serialname": "fin_dte", "mp3": True},
            "adj_vmr": {"serialname": "fin_adj_vmr", "mp3": True},
            "gsize": {"serialname": "fin_gsize", "mp3": True},
            "column_vapor": {"serialname": "fin_vapor", "mp3": True},
            "column_water": {"serialname": "fin_water", "mp3": True},
            "column_rain": {"serialname": "fin_rain", "mp3": True},
            "column_ice": {"serialname": "fin_ice", "mp3": True},
            "column_snow": {"serialname": "fin_snow", "mp3": True},
            "column_graupel": {"serialname": "fin_graupel", "mp3": True},
            "total_energy_wet_begin": {"serialname": "fin_ew0", "mp3": True},
            "total_water_wet_begin": {"serialname": "fin_ww0", "mp3": True},
            "total_energy_bot_wet_begin": {"serialname": "fin_bew0", "mp3": True},
            "total_water_bot_wet_begin": {"serialname": "fin_bww0", "mp3": True},
            "total_energy_dry_begin": {"serialname": "fin_ed0", "mp3": True},
            "total_water_dry_begin": {"serialname": "fin_wd0", "mp3": True},
            "total_energy_bot_dry_begin": {"serialname": "fin_bed0", "mp3": True},
            "total_water_bot_dry_begin": {"serialname": "fin_bwd0", "mp3": True},
            "total_energy_wet_end": {"serialname": "fin_ew", "mp3": True},
            "total_water_wet_end": {"serialname": "fin_ww", "mp3": True},
            "total_energy_bot_wet_end": {"serialname": "fin_bew", "mp3": True},
            "total_water_bot_wet_end": {"serialname": "fin_bww", "mp3": True},
            "total_energy_dry_end": {"serialname": "fin_ed", "mp3": True},
            "total_water_dry_end": {"serialname": "fin_wd", "mp3": True},
            "total_energy_bot_dry_end": {"serialname": "fin_bed", "mp3": True},
            "total_water_bot_dry_end": {"serialname": "fin_bwd", "mp3": True},
            "tzuv": {"serialname": "fin_tzuv", "mp3": True},
            "tzw": {"serialname": "fin_tzw", "mp3": True},
            "total_energy": {"serialname": "fin_te", "mp3": True},
            "column_energy_loss": {"serialname": "fin_te_loss", "mp3": True},
        }

        self.out_vars = {
            "qvapor": {"serialname": "fin_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "fin_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "fin_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "fin_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "fin_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "fin_qg", "kend": namelist.npz, "mp3": True},
            "delp": {"serialname": "fin_delp", "kend": namelist.npz, "mp3": True},
            "delz": {"serialname": "fin_delz", "kend": namelist.npz, "mp3": True},
            "pt": {"serialname": "fin_pt", "kend": namelist.npz, "mp3": True},
            "ua": {"serialname": "fin_ua", "kend": namelist.npz, "mp3": True},
            "va": {"serialname": "fin_va", "kend": namelist.npz, "mp3": True},
            "wa": {"serialname": "fin_wa", "kend": namelist.npz, "mp3": True},
            "column_energy_change": {"serialname": "fin_dte", "mp3": True},
            "adj_vmr": {"serialname": "fin_adj_vmr", "kend": namelist.npz, "mp3": True},
            "total_energy_wet_end": {
                "serialname": "fin_ew",
                "kend": namelist.npz,
                "mp3": True,
            },
            "total_water_wet_end": {
                "serialname": "fin_ww",
                "kend": namelist.npz,
                "mp3": True,
            },
            "total_energy_bot_wet_end": {"serialname": "fin_bew", "mp3": True},
            "total_water_bot_wet_end": {"serialname": "fin_bww", "mp3": True},
            "total_energy_dry_end": {
                "serialname": "fin_ed",
                "kend": namelist.npz,
                "mp3": True,
            },
            "total_water_dry_end": {
                "serialname": "fin_wd",
                "kend": namelist.npz,
                "mp3": True,
            },
            "total_energy_bot_dry_end": {"serialname": "fin_bed", "mp3": True},
            "total_water_bot_dry_end": {"serialname": "fin_bwd", "mp3": True},
            "tzuv": {"serialname": "fin_tzuv", "kend": namelist.npz, "mp3": True},
            "tzw": {"serialname": "fin_tzw", "kend": namelist.npz, "mp3": True},
            "total_energy": {"serialname": "fin_te", "kend": namelist.npz, "mp3": True},
            "column_energy_loss": {"serialname": "fin_te_loss", "mp3": True},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

        self.sizer = pace.util.SubtileGridSizer.from_tile_params(
            nx_tile=self.namelist.npx - 1,
            ny_tile=self.namelist.npy - 1,
            nz=self.namelist.npz,
            n_halo=3,
            extra_dim_lengths={},
            layout=self.namelist.layout,
        )

        self.quantity_factory = pace.util.QuantityFactory.from_backend(
            self.sizer, self.stencil_factory.backend
        )

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = FinalCalcs(self.stencil_factory, self.config, consv_te=False)

        compute_func(**inputs)

        return self.slice_output(inputs)
