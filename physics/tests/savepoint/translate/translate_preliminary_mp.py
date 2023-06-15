import pace.dsl
import pace.util
import pace.fv3core.stencils.basic_operations as basic
from pace.dsl.stencil import StencilFactory
from pace.physics._config import PhysicsConfig, MicroPhysicsConfig
from pace.physics.stencils.microphysics_v3.microphysics_v3 import reset_initial_values_and_make_copies, convert_virtual_to_true_temperature_and_calc_total_energy, moist_total_energy_and_water, convert_specific_to_mass_mixing_ratios_and_calculate_densities, cloud_nuclei_subgrid_and_relative_humidity
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py
from pace.dsl.typing import FloatField


class PrelimCalcs:
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

        self._copy_stencil = stencil_factory.from_origin_domain(
            basic.copy_defn,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._reset_initial_values_and_make_copies = stencil_factory.from_origin_domain(
            func=reset_initial_values_and_make_copies,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

        self._convert_virtual_to_true_temperature_and_calc_total_energy = (
            stencil_factory.from_origin_domain(
                func=convert_virtual_to_true_temperature_and_calc_total_energy,
                externals={
                    "consv_te": consv_te,
                    "do_inline_mp": self.do_inline_mp,
                    "hydrostatic": self.config.hydrostatic,
                    "c_air": self.config.c_air,
                    "c1_vap": self.config.c1_vap,
                    "c1_liq": self.config.c1_liq,
                    "c1_ice": self.config.c1_ice,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )
        )

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
                    "timestep": self.config.dt_atmos,
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
                    "timestep": self.config.dt_atmos,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        self._convert_specific_to_mass_mixing_ratios_and_calculate_densities = (
            stencil_factory.from_origin_domain(
                func=convert_specific_to_mass_mixing_ratios_and_calculate_densities,
                externals={
                    "do_inline_mp": self.config.do_inline_mp,
                },
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )
        )

        self._cloud_nuclei_subgrid_and_relative_humidity = (
            stencil_factory.from_origin_domain(
                func=cloud_nuclei_subgrid_and_relative_humidity,
                externals={
                    "prog_ccn": self.config.prog_ccn,
                    "ccn_l": self.config.ccn_l,
                    "ccn_o": self.config.ccn_o,
                    "dw_land": self.config.dw_land,
                    "dw_ocean": self.config.dw_ocean,
                    "rh_inc": self.config.rh_inc,
                    "rh_inr": self.config.rh_inr,
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
        density,
        pz,
        density_factor,
        pt,
        ua,
        va,
        wa,
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
        column_energy_change,
        cond,
        adj_vmr,
        gsize,
        hs,
        qcloud_cond_nuclei,
        qcloud_ice_nuclei,
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
        cloud_condensation_nuclei,
        cloud_ice_nuclei,
        h_var,
        rh_adj,
        rh_rain,
    ):
        self._reset_initial_values_and_make_copies(
            adj_vmr,
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            delp,
            pt,
            ua,
            va,
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
            column_energy_change,
            cond,
        )

        if not self.hydrostatic:
            self._copy_stencil(wa, w0)

        if self.do_inline_mp or self.consv_te:
            self._convert_virtual_to_true_temperature_and_calc_total_energy(
                qvapor,
                qliquid,
                qrain,
                qice,
                qsnow,
                qgraupel,
                pt,
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
                total_energy_wet_begin,
                total_water_wet_begin,
                total_energy_bot_wet_begin,
                total_water_bot_wet_begin,
            )

        self._convert_specific_to_mass_mixing_ratios_and_calculate_densities(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            delp,
            delz,
            pt,
            density,
            pz,
            density_factor,
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
                total_energy_dry_begin,
                total_water_dry_begin,
                total_energy_bot_dry_begin,
                total_water_bot_dry_begin,
            )

        self._cloud_nuclei_subgrid_and_relative_humidity(
            hs,
            qcloud_cond_nuclei,
            qcloud_ice_nuclei,
            density,
            cloud_condensation_nuclei,
            cloud_ice_nuclei,
            gsize,
            h_var,
            rh_adj,
            rh_rain,
        )


class TranslatePreliminaryCalculations(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "pp_qv", "mp3": True},
            "qliquid": {"serialname": "pp_ql", "mp3": True},
            "qrain": {"serialname": "pp_qr", "mp3": True},
            "qice": {"serialname": "pp_qi", "mp3": True},
            "qsnow": {"serialname": "pp_qs", "mp3": True},
            "qgraupel": {"serialname": "pp_qg", "mp3": True},
            "density": {"serialname": "pp_den", "mp3": True},
            "delp": {"serialname": "pp_delp", "mp3": True},
            "delz": {"serialname": "pp_delz", "mp3": True},
            "pz": {"serialname": "pp_pz", "mp3": True},
            "density_factor": {"serialname": "pp_denfac", "mp3": True},
            "pt": {"serialname": "pp_pt", "mp3": True},
            "ua": {"serialname": "pp_ua", "mp3": True},
            "va": {"serialname": "pp_va", "mp3": True},
            "wa": {"serialname": "pp_wa", "mp3": True},
            "qvapor0": {"serialname": "pp_qv0", "mp3": True},
            "qliquid0": {"serialname": "pp_ql0", "mp3": True},
            "qrain0": {"serialname": "pp_qr0", "mp3": True},
            "qice0": {"serialname": "pp_qi0", "mp3": True},
            "qsnow0": {"serialname": "pp_qs0", "mp3": True},
            "qgraupel0": {"serialname": "pp_qg0", "mp3": True},
            "dp0": {"serialname": "pp_dp0", "mp3": True},
            "pt0": {"serialname": "pp_pt0", "mp3": True},
            "u0": {"serialname": "pp_u0", "mp3": True},
            "v0": {"serialname": "pp_v0", "mp3": True},
            "w0": {"serialname": "pp_w0", "mp3": True},
            "column_energy_change": {"serialname": "pp_dte", "mp3": True},
            "cond": {"serialname": "pp_cond", "mp3": True},
            "adj_vmr": {"serialname": "pp_adj_vmr", "mp3": True},
            "gsize": {"serialname": "pp_gsize", "mp3": True},
            "hs": {"serialname": "pp_hs", "mp3": True},
            "qcloud_cond_nuclei": {"serialname": "pp_qnl", "mp3": True},
            "qcloud_ice_nuclei": {"serialname": "pp_qni", "mp3": True},
            "column_vapor": {"serialname": "pp_vapor", "mp3": True},
            "column_water": {"serialname": "pp_water", "mp3": True},
            "column_rain": {"serialname": "pp_rain", "mp3": True},
            "column_ice": {"serialname": "pp_ice", "mp3": True},
            "column_snow": {"serialname": "pp_snow", "mp3": True},
            "column_graupel": {"serialname": "pp_graupel", "mp3": True},
            "total_energy_wet_begin": {"serialname": "pp_ew0", "mp3": True},
            "total_water_wet_begin": {"serialname": "pp_ww0", "mp3": True},
            "total_energy_bot_wet_begin": {"serialname": "pp_bew0", "mp3": True},
            "total_water_bot_wet_begin": {"serialname": "pp_bww0", "mp3": True},
            "total_energy_dry_begin": {"serialname": "pp_ed0", "mp3": True},
            "total_water_dry_begin": {"serialname": "pp_wd0", "mp3": True},
            "total_energy_bot_dry_begin": {"serialname": "pp_bed0", "mp3": True},
            "total_water_bot_dry_begin": {"serialname": "pp_bwd0", "mp3": True},
            "cloud_condensation_nuclei": {"serialname": "pp_ccn", "mp3": True},
            "cloud_ice_nuclei": {"serialname": "pp_cin", "mp3": True},
            "h_var": {"serialname": "pp_h_var", "mp3": True},
            "rh_adj": {"serialname": "pp_rh_adj", "mp3": True},
            "rh_rain": {"serialname": "pp_rh_rain", "mp3": True},
        }

        self.out_vars = {
            "qvapor": {"serialname": "pp_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "pp_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "pp_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "pp_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "pp_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "pp_qg", "kend": namelist.npz, "mp3": True},
            "density": {"serialname": "pp_den", "kend": namelist.npz, "mp3": True},
            "delp": {"serialname": "pp_delp", "kend": namelist.npz, "mp3": True},
            "delz": {"serialname": "pp_delz", "kend": namelist.npz, "mp3": True},
            "pz": {"serialname": "pp_pz", "kend": namelist.npz, "mp3": True},
            "density_factor": {"serialname": "pp_denfac", "kend": namelist.npz, "mp3": True},
            "pt": {"serialname": "pp_pt", "kend": namelist.npz, "mp3": True},
            "ua": {"serialname": "pp_ua", "kend": namelist.npz, "mp3": True},
            "va": {"serialname": "pp_va", "kend": namelist.npz, "mp3": True},
            "wa": {"serialname": "pp_wa", "kend": namelist.npz, "mp3": True},
            "qvapor0": {"serialname": "pp_qv0", "kend": namelist.npz, "mp3": True},
            "qliquid0": {"serialname": "pp_ql0", "kend": namelist.npz, "mp3": True},
            "qrain0": {"serialname": "pp_qr0", "kend": namelist.npz, "mp3": True},
            "qice0": {"serialname": "pp_qi0", "kend": namelist.npz, "mp3": True},
            "qsnow0": {"serialname": "pp_qs0", "kend": namelist.npz, "mp3": True},
            "qgraupel0": {"serialname": "pp_qg0", "kend": namelist.npz, "mp3": True},
            "dp0": {"serialname": "pp_dp0", "kend": namelist.npz, "mp3": True},
            "pt0": {"serialname": "pp_pt0", "kend": namelist.npz, "mp3": True},
            "u0": {"serialname": "pp_u0", "kend": namelist.npz, "mp3": True},
            "v0": {"serialname": "pp_v0", "kend": namelist.npz, "mp3": True},
            "w0": {"serialname": "pp_w0", "kend": namelist.npz, "mp3": True},
            "column_energy_change": {"serialname": "pp_dte", "mp3": True},
            "cond": {"serialname": "pp_cond", "mp3": True},
            "adj_vmr": {"serialname": "pp_adj_vmr", "kend": namelist.npz, "mp3": True},
            "total_energy_wet_begin": {"serialname": "pp_ew0", "kend": namelist.npz, "mp3": True},
            "total_water_wet_begin": {"serialname": "pp_ww0", "kend": namelist.npz, "mp3": True},
            "total_energy_bot_wet_begin": {"serialname": "pp_bew0", "mp3": True},
            "total_water_bot_wet_begin": {"serialname": "pp_bww0", "mp3": True},
            "total_energy_dry_begin": {"serialname": "pp_ed0", "kend": namelist.npz, "mp3": True},
            "total_water_dry_begin": {"serialname": "pp_wd0", "kend": namelist.npz, "mp3": True},
            "total_energy_bot_dry_begin": {"serialname": "pp_bed0", "mp3": True},
            "total_water_bot_dry_begin": {"serialname": "pp_bwd0", "mp3": True},
            "cloud_condensation_nuclei": {"serialname": "pp_ccn", "kend": namelist.npz, "mp3": True},
            "cloud_ice_nuclei": {"serialname": "pp_cin", "kend": namelist.npz, "mp3": True},
            "h_var": {"serialname": "pp_h_var", "mp3": True},
            "rh_adj": {"serialname": "pp_rh_adj", "mp3": True},
            "rh_rain": {"serialname": "pp_rh_rain", "mp3": True},
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

        compute_func = PrelimCalcs(
            self.stencil_factory,
            self.config,
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
