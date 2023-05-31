import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
import pace.physics.stencils.microphysics_v3.physical_functions as physfun
import pace.util.constants as constants
from gt4py.cartesian.gtscript import __INLINED, FORWARD, computation, interval
from pace.physics.stencils.microphysics_v3.subgrid_z_proc import ( # noqa
    perform_instant_processes,
    cloud_condensation_evaporation,
    complete_freeze,
    wegener_bergeron_findeisen,
    freeze_bigg,
    deposit_and_sublimate_graupel,
    deposit_and_sublimate_ice,
    deposit_and_sublimate_snow,
)
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py
from pace.dsl.stencil import GridIndexing, StencilFactory
from pace.dsl.typing import FloatField, FloatFieldIJ


def vertical_subgrid_processes(
    qvapor: FloatField,
    qliquid: FloatField,
    qrain: FloatField,
    qice: FloatField,
    qsnow: FloatField,
    qgraupel: FloatField,
    temperature: FloatField,
    density: FloatField,
    density_factor: FloatField,
    delp: FloatField,
    cloud_condensation_nuclei: FloatField,
    cloud_ice_nuclei: FloatField,
    te: FloatField,
    cond: FloatFieldIJ,
    dep: FloatFieldIJ,
    reevap: FloatFieldIJ,
    sub: FloatFieldIJ,
    rh_adj: FloatFieldIJ,
):
    """"""
    from __externals__ import do_warm_rain_mp, do_wbf # noqa

    # with computation(FORWARD):
    #     with interval(-1, None):
    #         cond = 0
    #         dep = 0
    #         reevap = 0
    #         sub = 0

    with computation(FORWARD):
        with interval(...):
            # (
            #     q_liq,
            #     q_solid,
            #     cvm,
            #     te,
            #     lcpk,
            #     icpk,
            #     tcpk,
            #     tcp3,
            # ) = physfun.calc_heat_cap_and_latent_heat_coeff(
            #     qvapor, qliquid, qrain, qice, qsnow, qgraupel, temperature
            # )

            if __INLINED(not do_warm_rain_mp):
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
                    dep,
                    reevap,
                    sub,
                ) = perform_instant_processes(
                    qvapor,
                    qliquid,
                    qrain,
                    qice,
                    qsnow,
                    qgraupel,
                    temperature,
                    density,
                    delp,
                    te,
                    rh_adj,
                    dep,
                    reevap,
                    sub,
                )

            # (
            #     qvapor,
            #     qliquid,
            #     qrain,
            #     qice,
            #     qsnow,
            #     qgraupel,
            #     temperature,
            #     cvm,
            #     lcpk,
            #     icpk,
            #     tcpk,
            #     tcp3,
            #     cond,
            #     reevap,
            # ) = cloud_condensation_evaporation(
            #     qvapor,
            #     qliquid,
            #     qrain,
            #     qice,
            #     qsnow,
            #     qgraupel,
            #     temperature,
            #     delp,
            #     density,
            #     te,
            #     tcp3,
            #     cond,
            #     reevap,
            # )

            # if __INLINED(not do_warm_rain_mp):
            #     (
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         temperature,
            #         cvm,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #     ) = complete_freeze(
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         temperature,
            #         cvm,
            #         te,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #     )

            #     if __INLINED(do_wbf):
            #         (
            #             qvapor,
            #             qliquid,
            #             qrain,
            #             qice,
            #             qsnow,
            #             qgraupel,
            #             temperature,
            #             cvm,
            #             lcpk,
            #             icpk,
            #             tcpk,
            #             tcp3,
            #         ) = wegener_bergeron_findeisen(
            #             qvapor,
            #             qliquid,
            #             qrain,
            #             qice,
            #             qsnow,
            #             qgraupel,
            #             temperature,
            #             density,
            #             cvm,
            #             te,
            #             lcpk,
            #             icpk,
            #             tcpk,
            #             tcp3,
            #         )

            #     (
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         cloud_condensation_nuclei,
            #         temperature,
            #         cvm,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #     ) = freeze_bigg(
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         cloud_condensation_nuclei,
            #         temperature,
            #         density,
            #         cvm,
            #         te,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #     )

            #     (
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         cloud_ice_nuclei,
            #         temperature,
            #         cvm,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #         dep,
            #         sub,
            #     ) = deposit_and_sublimate_ice(
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         cloud_ice_nuclei,
            #         temperature,
            #         delp,
            #         density,
            #         cvm,
            #         te,
            #         dep,
            #         sub,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #     )

            #     (
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         temperature,
            #         cvm,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #         dep,
            #         sub,
            #     ) = deposit_and_sublimate_snow(
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         temperature,
            #         delp,
            #         density,
            #         density_factor,
            #         cvm,
            #         te,
            #         dep,
            #         sub,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #     )

            #     (
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         temperature,
            #         cvm,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #         dep,
            #         sub,
            #     ) = deposit_and_sublimate_graupel(
            #         qvapor,
            #         qliquid,
            #         qrain,
            #         qice,
            #         qsnow,
            #         qgraupel,
            #         temperature,
            #         delp,
            #         density,
            #         density_factor,
            #         cvm,
            #         te,
            #         dep,
            #         sub,
            #         lcpk,
            #         icpk,
            #         tcpk,
            #         tcp3,
            #     )


class SubSubgridProcesses:
    """
    subgrid_z_proc in Fortran
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        config,
        timestep: float,
    ):
        self._idx: GridIndexing = stencil_factory.grid_indexing

        if config.do_hail:
            mu_g = config.muh
            blin_g = config.blinh
        else:
            mu_g = config.mug
            blin_g = config.bling

        if config.inflag not in [1, 2, 3, 4, 5]:
            raise ValueError(
                f"Ice Nucleation Flag must be an integer from 1 to 5"
                f"not {config.inflag}"
            )

        if config.igflag not in [1, 2, 3, 4]:
            raise ValueError(
                f"Ice Generation Flag must be an int from 1 to 4, got {config.igflag}"
            )

        self._vertical_subgrid_processes = stencil_factory.from_origin_domain(
            func=vertical_subgrid_processes,
            externals={
                "timestep": timestep,
                "c1_ice": config.c1_ice,
                "c1_liq": config.c1_liq,
                "c1_vap": config.c1_vap,
                "d1_ice": config.d1_ice,
                "d1_vap": config.d1_vap,
                "li00": config.li00,
                "li20": config.li20,
                "lv00": config.lv00,
                "t_wfr": config.t_wfr,
                "tice": constants.TICE0,
                "t_min": config.t_min,
                "t_sub": config.t_sub,
                "do_cond_timescale": config.do_cond_timescale,
                "rh_fac": config.rh_fac,
                "rhc_cevap": config.rhc_cevap,
                "tau_l2v": config.tau_l2v,
                "tau_v2l": config.tau_v2l,
                "use_rhc_cevap": config.use_rhc_cevap,
                "qi0_crt": config.qi0_crt,
                "tau_wbf": config.tau_wbf,
                "do_psd_water_num": config.do_psd_water_num,
                "muw": config.muw,
                "pcaw": config.pcaw,
                "pcbw": config.pcbw,
                "do_psd_ice_num": config.do_psd_ice_num,
                "igflag": config.igflag,
                "inflag": config.inflag,
                "is_fac": config.is_fac,
                "mui": config.mui,
                "pcai": config.pcai,
                "pcbi": config.pcbi,
                "prog_ccn": config.prog_ccn,
                "qi_lim": config.qi_lim,
                "blins": config.blins,
                "mus": config.mus,
                "cssub_1": config.cssub_1,
                "cssub_2": config.cssub_2,
                "cssub_3": config.cssub_3,
                "cssub_4": config.cssub_4,
                "cssub_5": config.cssub_5,
                "ss_fac": config.ss_fac,
                "cgsub_1": config.cgsub_1,
                "cgsub_2": config.cgsub_2,
                "cgsub_3": config.cgsub_3,
                "cgsub_4": config.cgsub_4,
                "cgsub_5": config.cgsub_5,
                "bling": config.bling,
                "mug": config.mug,
                "gs_fac": config.gs_fac,
                "do_warm_rain_mp": config.do_warm_rain_mp,
                "do_wbf": config.do_wbf,
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
        density: FloatField,
        density_factor: FloatField,
        delp: FloatField,
        cloud_condensation_nuclei: FloatField,
        cloud_ice_nuclei: FloatField,
        te: FloatField,
        cond: FloatFieldIJ,
        dep: FloatFieldIJ,
        reevap: FloatFieldIJ,
        sub: FloatFieldIJ,
        rh_adj: FloatFieldIJ,
    ):
        """
        Temperature sentive high vertical resolution processes
        Args:
            qvapor (inout):
            qliquid (inout):
            qrain (inout):
            qice (inout):
            qsnow (inout):
            qgraupel (inout):
            temperature (inout):
            density (in):
            density_factor (in):
            delp (in):
            cloud_condensation_nuclei (inout):
            cloud_ice_nuclei (inout):
            te (in):
            cond (out):
            dep (out):
            reevap (out):
            sub (out):
            rh_adj (in):
        """

        self._vertical_subgrid_processes(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temperature,
            density,
            density_factor,
            delp,
            cloud_condensation_nuclei,
            cloud_ice_nuclei,
            te,
            cond,
            dep,
            reevap,
            sub,
            rh_adj,
        )


class TranslateSubgridZSubs(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "szs_qv", "mp3": True},
            "qliquid": {"serialname": "szs_ql", "mp3": True},
            "qrain": {"serialname": "szs_qr", "mp3": True},
            "qice": {"serialname": "szs_qi", "mp3": True},
            "qsnow": {"serialname": "szs_qs", "mp3": True},
            "qgraupel": {"serialname": "szs_qg", "mp3": True},
            "temperature": {"serialname": "szs_pt", "mp3": True},
            "density": {"serialname": "szs_den", "mp3": True},
            "density_factor": {"serialname": "szs_denfac", "mp3": True},
            "delp": {"serialname": "szs_delp", "mp3": True},
            "rh_adj": {"serialname": "szs_rh_adj", "mp3": True},
            "cloud_condensation_nuclei": {"serialname": "szs_ccn", "mp3": True},
            "cloud_ice_nuclei": {"serialname": "szs_cin", "mp3": True},
            "te": {"serialname": "szs_te", "mp3": True},
            "cond": {"serialname": "szs_cond", "mp3": True},
            "dep": {"serialname": "szs_dep", "mp3": True},
            "reevap": {"serialname": "szs_reevap", "mp3": True},
            "sub": {"serialname": "szs_sub", "mp3": True},
        }

        self.in_vars["parameters"] = [
            "dt",
        ]

        self.out_vars = {
            "qvapor": {"serialname": "szs_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "szs_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "szs_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "szs_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "szs_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "szs_qg", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "szs_pt", "kend": namelist.npz, "mp3": True},
            "cloud_condensation_nuclei": {
                "serialname": "szs_ccn",
                "kend": namelist.npz,
                "mp3": True,
            },
            "cloud_ice_nuclei": {
                "serialname": "szs_cin",
                "kend": namelist.npz,
                "mp3": True,
            },
            "cond": {"serialname": "szs_cond", "kend": namelist.npz, "mp3": True},
            "dep": {"serialname": "szs_dep", "kend": namelist.npz, "mp3": True},
            "reevap": {"serialname": "szs_reevap", "kend": namelist.npz, "mp3": True},
            "sub": {"serialname": "szs_sub", "kend": namelist.npz, "mp3": True},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = SubSubgridProcesses(
            self.stencil_factory,
            self.config,
            timestep=inputs.pop("dt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
