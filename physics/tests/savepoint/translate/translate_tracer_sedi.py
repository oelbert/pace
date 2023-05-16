import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.sedimentation import (
    adjust_fluxes,
    calc_edge_and_terminal_height,
    calc_terminal_velocity_ice,
    calc_terminal_velocity_rsg,
    sedi_melt,
)
from pace.physics.stencils.microphysics_v3.terminal_fall import TerminalFall
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py
from pace.util import X_DIM, Y_DIM, Z_DIM


class TracerSedimentation:
    def __init__(
        self,
        stencil_factory: pace.dsl.StencilFactory,
        quantity_factory: pace.util.QuantityFactory,
        config,
        timestep,
        tracer,
    ):
        self._tracer = tracer

        self._idx = stencil_factory.grid_indexing
        self._is_ = self._idx.isc
        self._ie = self._idx.iec
        self._js = self._idx.jsc
        self._je = self._idx.jec
        self._ks = 0
        self._ke = config.npz - 1
        self.c1_vap = config.c1_vap
        self.c1_liq = config.c1_liq
        self.c1_ice = config.c1_ice

        if config.do_hail:
            self.tvag = config.tvah
            self.tvbg = config.tvbh
            self.mug = config.muh
            self.bling = config.blinh
        else:
            self.tvag = config.tvag
            self.tvbg = config.tvbg
            self.mug = config.mug
            self.bling = config.bling

        assert config.ifflag in [
            1,
            2,
        ], f"Ice Formation Flag must be 1 or 2 not {config.ifflag}"

        self._timestep = timestep
        self.config = config
        self.li00 = config.li00

        # allocate internal storages
        def make_quantity():
            return quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="unknown")

        self._icpk = make_quantity()
        self._cvm = make_quantity()

        if self._tracer == "qi":
            if config.do_psd_ice_fall is False:
                self._calc_terminal_ice_velocity = stencil_factory.from_origin_domain(
                    func=calc_terminal_velocity_ice,
                    externals={
                        "ifflag": config.ifflag,
                        "constant_v": config.const_vi,
                        "v_fac": config.vi_fac,
                        "v_max": config.vi_max,
                        "aa": -4.14122e-5,
                        "bb": -0.00538922,
                        "cc": -0.0516344,
                        "dd": 0.00216078,
                        "ee": 1.9714,
                    },
                    origin=self._idx.origin_compute(),
                    domain=self._idx.domain_compute(),
                )

        if False in [
            config.const_vr,
            config.const_vi,
            config.const_vs,
            config.const_vg,
        ]:
            self._calc_terminal_rsg_velocity = stencil_factory.from_origin_domain(
                func=calc_terminal_velocity_rsg,
                externals={"const_v": False},
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        if True in [config.const_vr, config.const_vi, config.const_vs, config.const_vg]:
            self._calc_terminal_rsg_velocity_const = stencil_factory.from_origin_domain(
                func=calc_terminal_velocity_rsg,
                externals={"const_v": True},
                origin=self._idx.origin_compute(),
                domain=self._idx.domain_compute(),
            )

        self._calc_edge_and_terminal_height = stencil_factory.from_origin_domain(
            func=calc_edge_and_terminal_height,
            externals={"timestep": self._timestep},
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(add=(0, 0, 1)),
        )

        self._terminal_fall = TerminalFall(
            stencil_factory, quantity_factory, config, timestep
        )

        self._adjust_fluxes = stencil_factory.from_origin_domain(
            func=adjust_fluxes,
            origin=self._idx.origin_compute(),
            domain=self._idx.domain_compute(),
        )

    def __call__(
        self,
        qvapor,
        qliquid,
        qrain,
        qice,
        qsnow,
        qgraupel,
        temperature,
        delp,
        delz,
        density,
        density_factor,
        ua,
        va,
        wa,
        column_energy_change,
        preflux_tracer,
        vterminal_tracer,
        column_water,
        column_rain,
        column_ice,
        column_snow,
        column_graupel,
        z_surface,
        z_edge,
        z_terminal,
    ):
        if self._tracer == "qi":
            # Terminal fall and melting of falling cloud ice into rain:
            if self.config.do_psd_ice_fall:
                if self.config.const_vi is False:
                    self._calc_terminal_rsg_velocity(
                        qice,
                        density,
                        density_factor,
                        vterminal_tracer,
                        self.config.vi_fac,
                        self.config.tvai,
                        self.config.tvbi,
                        self.config.mui,
                        self.config.blini,
                        self.config.vi_max,
                    )
                else:
                    self._calc_terminal_rsg_velocity_const(
                        qice,
                        density,
                        density_factor,
                        vterminal_tracer,
                        self.config.vi_fac,
                        self.config.tvai,
                        self.config.tvbi,
                        self.config.mui,
                        self.config.blini,
                        self.config.vi_max,
                    )
            else:
                self._calc_terminal_ice_velocity(
                    qice, temperature, density, vterminal_tracer
                )

            self._calc_edge_and_terminal_height(
                z_surface,
                z_edge,
                z_terminal,
                delz,
                vterminal_tracer,
            )

            if self.config.do_sedi_melt:
                sedi_melt(
                    qvapor,
                    qliquid,
                    qrain,
                    qice,
                    qsnow,
                    qgraupel,
                    self._cvm.data[:],
                    temperature,
                    delp,
                    z_edge.data[:],
                    z_terminal.data[:],
                    z_surface.data[:],
                    self._timestep,
                    vterminal_tracer,
                    column_rain,
                    self.config.tau_imlt,
                    self._icpk.data[:],
                    self.li00,
                    self.c1_vap,
                    self.c1_liq,
                    self.c1_ice,
                    self._ks,
                    self._ke,
                    self._is_,
                    self._ie,
                    self._js,
                    self._je,
                    "ice",
                )

            self._terminal_fall(
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
                vterminal_tracer,
                z_edge,
                z_terminal,
                preflux_tracer,
                column_ice,
                column_energy_change,
                "ice",
            )

        elif self._tracer == "qs":
            # Terminal fall and melting of falling snow into rain
            if self.config.const_vs is False:
                self._calc_terminal_rsg_velocity(
                    qsnow,
                    density,
                    density_factor,
                    vterminal_tracer,
                    self.config.vs_fac,
                    self.config.tvas,
                    self.config.tvbs,
                    self.config.mus,
                    self.config.blins,
                    self.config.vs_max,
                )
            else:
                self._calc_terminal_rsg_velocity_const(
                    qsnow,
                    density,
                    density_factor,
                    vterminal_tracer,
                    self.config.vs_fac,
                    self.config.tvas,
                    self.config.tvbs,
                    self.config.mus,
                    self.config.blins,
                    self.config.vs_max,
                )

            self._calc_edge_and_terminal_height(
                z_surface,
                z_edge,
                z_terminal,
                delz,
                vterminal_tracer,
            )

            if self.config.do_sedi_melt:
                sedi_melt(
                    qvapor,
                    qliquid,
                    qrain,
                    qice,
                    qsnow,
                    qgraupel,
                    self._cvm.data[:],
                    temperature,
                    delp,
                    z_edge.data[:],
                    z_terminal.data[:],
                    z_surface.data[:],
                    self._timestep,
                    vterminal_tracer,
                    column_rain,
                    self.config.tau_smlt,
                    self._icpk.data[:],
                    self.li00,
                    self.c1_vap,
                    self.c1_liq,
                    self.c1_ice,
                    self._ks,
                    self._ke,
                    self._is_,
                    self._ie,
                    self._js,
                    self._je,
                    "snow",
                )

            self._terminal_fall(
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
                vterminal_tracer,
                z_edge,
                z_terminal,
                preflux_tracer,
                column_snow,
                column_energy_change,
                "snow",
            )

        elif self._tracer == "qg":
            # Terminal fall and melting of falling graupel into rain
            if self.config.const_vg is False:
                self._calc_terminal_rsg_velocity(
                    qgraupel,
                    density,
                    density_factor,
                    vterminal_tracer,
                    self.config.vg_fac,
                    self.tvag,
                    self.tvbg,
                    self.mug,
                    self.bling,
                    self.config.vg_max,
                )
            else:
                self._calc_terminal_rsg_velocity_const(
                    qgraupel,
                    density,
                    density_factor,
                    vterminal_tracer,
                    self.config.vg_fac,
                    self.tvag,
                    self.tvbg,
                    self.mug,
                    self.bling,
                    self.config.vg_max,
                )

            self._calc_edge_and_terminal_height(
                z_surface,
                z_edge,
                z_terminal,
                delz,
                vterminal_tracer,
            )

            if self.config.do_sedi_melt:
                sedi_melt(
                    qvapor,
                    qliquid,
                    qrain,
                    qice,
                    qsnow,
                    qgraupel,
                    self._cvm.data[:],
                    temperature,
                    delp,
                    z_edge.data[:],
                    z_terminal.data[:],
                    z_surface.data[:],
                    self._timestep,
                    vterminal_tracer,
                    column_rain,
                    self.config.tau_gmlt,
                    self._icpk.data[:],
                    self.li00,
                    self.c1_vap,
                    self.c1_liq,
                    self.c1_ice,
                    self._ks,
                    self._ke,
                    self._is_,
                    self._ie,
                    self._js,
                    self._je,
                    "graupel",
                )

            self._terminal_fall(
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
                vterminal_tracer,
                z_edge,
                z_terminal,
                preflux_tracer,
                column_graupel,
                column_energy_change,
                "graupel",
            )

        elif self._tracer == "ql":
            # Terminal fall of cloud water
            if self.config.do_psd_water_fall:
                if self.config.const_vw is False:
                    self._calc_terminal_rsg_velocity(
                        qliquid,
                        density,
                        density_factor,
                        vterminal_tracer,
                        self.config.vw_fac,
                        self.config.tvaw,
                        self.config.tvbw,
                        self.config.muw,
                        self.config.blinw,
                        self.config.vw_max,
                    )
                else:
                    self._calc_terminal_rsg_velocity_const(
                        qliquid,
                        density,
                        density_factor,
                        vterminal_tracer,
                        self.config.vw_fac,
                        self.config.tvaw,
                        self.config.tvbw,
                        self.config.muw,
                        self.config.blinw,
                        self.config.vw_max,
                    )

                self._calc_edge_and_terminal_height(
                    z_surface,
                    z_edge,
                    z_terminal,
                    delz,
                    vterminal_tracer,
                )

                self._terminal_fall(
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
                    vterminal_tracer,
                    z_edge,
                    z_terminal,
                    preflux_tracer,
                    column_water,
                    column_energy_change,
                    "liquid",
                )

        elif self._tracer == "qr":
            # terminal fall of rain
            if self.config.const_vs is False:
                self._calc_terminal_rsg_velocity(
                    qrain,
                    density,
                    density_factor,
                    vterminal_tracer,
                    self.config.vr_fac,
                    self.config.tvar,
                    self.config.tvbr,
                    self.config.mur,
                    self.config.blinr,
                    self.config.vr_max,
                )
            else:
                self._calc_terminal_rsg_velocity_const(
                    qrain,
                    density,
                    density_factor,
                    vterminal_tracer,
                    self.config.vr_fac,
                    self.config.tvar,
                    self.config.tvbr,
                    self.config.mur,
                    self.config.blinr,
                    self.config.vr_max,
                )

            self._calc_edge_and_terminal_height(
                z_surface,
                z_edge,
                z_terminal,
                delz,
                vterminal_tracer,
            )

            self._terminal_fall(
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
                vterminal_tracer,
                z_edge,
                z_terminal,
                preflux_tracer,
                column_rain,
                column_energy_change,
                "rain",
            )

        else:
            raise ValueError(f"Tracer {self._tracer} not recognized")

        self._adjust_fluxes(preflux_tracer)


class TranslateTracerSed(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)

        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "ts_qv", "mp3": True},
            "qliquid": {"serialname": "ts_ql", "mp3": True},
            "qrain": {"serialname": "ts_qr", "mp3": True},
            "qice": {"serialname": "ts_qi", "mp3": True},
            "qsnow": {"serialname": "ts_qs", "mp3": True},
            "qgraupel": {"serialname": "ts_qg", "mp3": True},
            "temperature": {"serialname": "ts_pt", "mp3": True},
            "delp": {"serialname": "ts_delp", "mp3": True},
            "delz": {"serialname": "ts_delz", "mp3": True},
            "density": {"serialname": "ts_den", "mp3": True},
            "density_factor": {"serialname": "ts_denfac", "mp3": True},
            "ua": {"serialname": "ts_ua", "mp3": True},
            "va": {"serialname": "ts_va", "mp3": True},
            "wa": {"serialname": "ts_wa", "mp3": True},
            "column_energy_change": {"serialname": "ts_dte", "mp3": True},
            "preflux_tracer": {"serialname": "ts_pf", "mp3": True},
            "vterminal_tracer": {"serialname": "ts_vt", "mp3": True},
            "column_water": {"serialname": "ts_w1", "mp3": True},
            "column_rain": {"serialname": "ts_r1", "mp3": True},
            "column_ice": {"serialname": "ts_i1", "mp3": True},
            "column_snow": {"serialname": "ts_s1", "mp3": True},
            "column_graupel": {"serialname": "ts_g1", "mp3": True},
            "z_edge": {"serialname": "ts_ze", "mp3": True},
            "z_terminal": {"serialname": "ts_zt", "mp3": True},
            "z_surface": {"serialname": "ts_zs", "mp3": True},
        }
        self.in_vars["parameters"] = ["dt"]
        self.out_vars = {
            "qvapor": {"serialname": "ts_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "ts_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "ts_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "ts_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "ts_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "ts_qg", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "ts_pt", "kend": namelist.npz, "mp3": True},
            "ua": {"serialname": "ts_ua", "kend": namelist.npz, "mp3": True},
            "va": {"serialname": "ts_va", "kend": namelist.npz, "mp3": True},
            "wa": {"serialname": "ts_wa", "kend": namelist.npz, "mp3": True},
            "preflux_tracer": {
                "serialname": "ts_pf",
                "kend": namelist.npz,
                "mp3": True,
            },
            "vterminal_tracer": {
                "serialname": "ts_vt",
                "kend": namelist.npz,
                "mp3": True,
            },
            "column_water": {"serialname": "ts_w1", "mp3": True},
            "column_rain": {"serialname": "ts_r1", "mp3": True},
            "column_ice": {"serialname": "ts_i1", "mp3": True},
            "column_snow": {"serialname": "ts_s1", "mp3": True},
            "column_graupel": {"serialname": "ts_g1", "mp3": True},
            "column_energy_change": {
                "serialname": "ts_dte",
                "kend": namelist.npz,
                "mp3": True,
            },
            "z_edge": {"serialname": "ts_ze", "kend": namelist.npz + 1, "mp3": True},
            "z_terminal": {
                "serialname": "ts_zt",
                "kend": namelist.npz + 1,
                "mp3": True,
            },
            "z_surface": {"serialname": "ts_zs", "kend": namelist.npz + 1, "mp3": True},
        }

        self.stencil_factory = stencil_factory
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

        sizer = pace.util.SubtileGridSizer.from_tile_params(
            nx_tile=self.namelist.npx - 1,
            ny_tile=self.namelist.npy - 1,
            nz=self.namelist.npz,
            n_halo=3,
            extra_dim_lengths={},
            layout=self.namelist.layout,
        )

        self.quantity_factory = pace.util.QuantityFactory.from_backend(
            sizer, self.stencil_factory.backend
        )

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        ts_q = "qs"

        compute_func = TracerSedimentation(
            self.stencil_factory,
            self.quantity_factory,
            self.config,
            timestep=inputs.pop("dt"),
            tracer=ts_q,
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
