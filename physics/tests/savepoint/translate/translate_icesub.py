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
from pace.physics.stencils.microphysics_v3.ice_cloud import (  # noqa
    accrete_graupel_with_cloud_water_and_rain,
    accrete_graupel_with_ice,
    accrete_graupel_with_snow,
    accrete_snow_with_ice,
    accrete_snow_with_rain_and_freeze_to_graupel,
    autoconvert_ice_to_snow,
    autoconvert_snow_to_graupel,
    melt_snow,
)
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


def test_func_stencil(
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
    cvm: FloatField,
    te: FloatField,
    lcpk: FloatField,
    icpk: FloatField,
    tcpk: FloatField,
    tcp3: FloatField,
    di: FloatField,
):
    from __externals__ import z_slope_ice  # noqa

    # with computation(PARALLEL), interval(...):
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
    #     ) = melt_cloud_ice(
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
    #     ) = freeze_cloud_water(
    #         qvapor,
    #         qliquid,
    #         qrain,
    #         qice,
    #         qsnow,
    #         qgraupel,
    #         temperature,
    #         density,
    #         cvm,
    #         te,
    #         lcpk,
    #         icpk,
    #         tcpk,
    #         tcp3,
    #     )
    # with computation(FORWARD):
    #     with interval(0, 1):
    #         if __INLINED(z_slope_ice):
    #             # linear_prof
    #             di = 0.0
    #     with interval(1, None):
    #         if __INLINED(z_slope_ice):
    #             dq = 0.5 * (qice - qice[0, 0, -1])
    # with computation(FORWARD):
    #     with interval(1, -1):
    #         if __INLINED(z_slope_ice):
    #             # Use twice the strength of the
    #             # positive definiteness limiter (lin et al 1994)
    #             di = 0.5 * min(abs(dq + dq[0, 0, +1]), 0.5 * qice[0, 0, 0])
    #             if dq * dq[0, 0, +1] <= 0.0:
    #                 if dq > 0.0:  # Local maximum
    #                     di = min(di, min(dq, -dq[0, 0, +1]))
    #                 else:  # Local minimum
    #                     di = 0.0
    #     with interval(-1, None):
    #         if __INLINED(z_slope_ice):
    #             di = 0.0
    # with computation(PARALLEL), interval(...):
    #     if __INLINED(z_slope_ice):
    #         # Impose a presumed background horizontal variability that is
    #         # proportional to the value itself
    #         di = max(di, max(0.0, h_var * qice))
    #     else:
    #         di = max(0.0, h_var * qice)

    with computation(PARALLEL), interval(...):
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
        ) = melt_snow(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temperature,
            density,
            density_factor,
            vterminal_water,
            vterminal_rain,
            vterminal_snow,
            cvm,
            te,
            lcpk,
            icpk,
            tcpk,
            tcp3,
        )
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
    #     ) = melt_graupel(
    #         qvapor,
    #         qliquid,
    #         qrain,
    #         qice,
    #         qsnow,
    #         qgraupel,
    #         temperature,
    #         density,
    #         density_factor,
    #         vterminal_water,
    #         vterminal_rain,
    #         vterminal_graupel,
    #         cvm,
    #         te,
    #         lcpk,
    #         icpk,
    #         tcpk,
    #         tcp3,
    #     )
    #     qice, qsnow = accrete_snow_with_ice(
    #         qice,
    #         qsnow,
    #         temperature,
    #         density,
    #         density_factor,
    #         vterminal_ice,
    #         vterminal_snow,
    #     )
    #     (qice, qsnow, di, temperature) = autoconvert_ice_to_snow(
    #         qice, qsnow, temperature, density, di
    #     )
    #     qice, qgraupel = accrete_graupel_with_ice(
    #         qice,
    #         qgraupel,
    #         temperature,
    #         density,
    #         density_factor,
    #         vterminal_ice,
    #         vterminal_graupel,
    #     )

    # with computation(PARALLEL), interval(...):
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
    #     ) = accrete_snow_with_rain_and_freeze_to_graupel(
    #         qvapor,
    #         qliquid,
    #         qrain,
    #         qice,
    #         qsnow,
    #         qgraupel,
    #         temperature,
    #         density,
    #         vterminal_rain,
    #         vterminal_snow,
    #         te,
    #         cvm,
    #         lcpk,
    #         icpk,
    #         tcpk,
    #         tcp3,
    #     )

    # qsnow, qgraupel = accrete_graupel_with_snow(
    #     qsnow,
    #     qgraupel,
    #     temperature,
    #     density,
    #     vterminal_snow,
    #     vterminal_graupel,
    # )

    # qsnow, qgraupel = autoconvert_snow_to_graupel(
    #     qsnow, qgraupel, temperature, density
    # )

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
    # ) = accrete_graupel_with_cloud_water_and_rain(
    #     qvapor,
    #     qliquid,
    #     qrain,
    #     qice,
    #     qsnow,
    #     qgraupel,
    #     temperature,
    #     density,
    #     density_factor,
    #     vterminal_rain,
    #     vterminal_graupel,
    #     te,
    #     cvm,
    #     lcpk,
    #     icpk,
    #     tcpk,
    #     tcp3,
    # )


class IceFunction:
    def __init__(
        self,
        stencil_factory: pace.dsl.StencilFactory,
        config,
        timestep,
    ):
        self._idx = stencil_factory.grid_indexing
        if config.do_hail:
            mu_g = config.muh
            blin_g = config.blinh
        else:
            mu_g = config.mug
            blin_g = config.bling

        self._test_func_stencil = stencil_factory.from_origin_domain(
            func=test_func_stencil,
            externals={
                "timestep": timestep,
                "z_slope_ice": config.z_slope_ice,
                "vdiffflag": config.vdiffflag,
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
                "ql_mlt": config.ql_mlt,
                "tau_imlt": config.tau_imlt,
                "tice_mlt": config.tice_mlt,
                "tau_i2s": config.tau_i2s,
                "qi0_crt": config.qi0_crt,
                "qs0_crt": config.qs0_crt,
                "acc0": config.acc[0],
                "acc1": config.acc[1],
                "acc2": config.acc[2],
                "acc3": config.acc[3],
                "acc4": config.acc[4],
                "acc5": config.acc[5],
                "acc6": config.acc[6],
                "acc7": config.acc[7],
                "acc12": config.acc[12],
                "acc13": config.acc[13],
                "acc14": config.acc[14],
                "acc15": config.acc[15],
                "acc16": config.acc[16],
                "acc17": config.acc[17],
                "acc18": config.acc[18],
                "acc19": config.acc[19],
                "acco_0_0": config.acco[0][0],
                "acco_1_0": config.acco[1][0],
                "acco_2_0": config.acco[2][0],
                "acco_0_1": config.acco[0][1],
                "acco_1_1": config.acco[1][1],
                "acco_2_1": config.acco[2][1],
                "acco_0_2": config.acco[0][2],
                "acco_1_2": config.acco[1][2],
                "acco_2_2": config.acco[2][2],
                "acco_0_3": config.acco[0][3],
                "acco_1_3": config.acco[1][3],
                "acco_2_3": config.acco[2][3],
                "acco_0_6": config.acco[0][6],
                "acco_1_6": config.acco[1][6],
                "acco_2_6": config.acco[2][6],
                "acco_0_7": config.acco[0][7],
                "acco_1_7": config.acco[1][7],
                "acco_2_7": config.acco[2][7],
                "acco_0_8": config.acco[0][8],
                "acco_1_8": config.acco[1][8],
                "acco_2_8": config.acco[2][8],
                "acco_0_9": config.acco[0][9],
                "acco_1_9": config.acco[1][9],
                "acco_2_9": config.acco[2][9],
                "cracs": config.cracs,
                "csacr": config.csacr,
                "csacw": config.csacw,
                "cgacw": config.cgacw,
                "cgacr": config.cgacr,
                "cgacs": config.cgacs,
                "csaci": config.csaci,
                "cgaci": config.cgaci,
                "csmlt_1": config.csmlt_1,
                "csmlt_2": config.csmlt_2,
                "csmlt_3": config.csmlt_3,
                "csmlt_4": config.csmlt_4,
                "cgmlt_1": config.cgmlt_1,
                "cgmlt_2": config.cgmlt_2,
                "cgmlt_3": config.cgmlt_3,
                "cgmlt_4": config.cgmlt_4,
                "cgfr_1": config.cgfr_1,
                "cgfr_2": config.cgfr_2,
                "do_new_acc_water": config.do_new_acc_water,
                "do_new_acc_ice": config.do_new_acc_ice,
                "mur": config.mur,
                "blinr": config.blinr,
                "mus": config.mus,
                "blins": config.blins,
                "qs_mlt": config.qs_mlt,
                "mug": mu_g,
                "bling": blin_g,
                "fi2s_fac": config.fi2s_fac,
                "fi2g_fac": config.fi2g_fac,
                "fs2g_fac": config.fs2g_fac,
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
        vterminal_water: FloatField,
        vterminal_rain: FloatField,
        vterminal_ice: FloatField,
        vterminal_snow: FloatField,
        vterminal_graupel: FloatField,
        h_var: FloatFieldIJ,
        cvm: FloatField,
        te8: FloatField,
        lcpk: FloatField,
        icpk: FloatField,
        tcpk: FloatField,
        tcp3: FloatField,
        di: FloatField,
    ):
        self._test_func_stencil(
            qvapor,
            qliquid,
            qrain,
            qice,
            qsnow,
            qgraupel,
            temperature,
            density,
            density_factor,
            vterminal_water,
            vterminal_rain,
            vterminal_ice,
            vterminal_snow,
            vterminal_graupel,
            h_var,
            cvm,
            te8,
            lcpk,
            icpk,
            tcpk,
            tcp3,
            di,
        )


class TranslateIceSubFunc(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "isub_qv", "mp3": True},
            "qliquid": {"serialname": "isub_ql", "mp3": True},
            "qrain": {"serialname": "isub_qr", "mp3": True},
            "qice": {"serialname": "isub_qi", "mp3": True},
            "qsnow": {"serialname": "isub_qs", "mp3": True},
            "qgraupel": {"serialname": "isub_qg", "mp3": True},
            "temperature": {"serialname": "isub_pt", "mp3": True},
            "density": {"serialname": "isub_den", "mp3": True},
            "density_factor": {"serialname": "isub_denfac", "mp3": True},
            "vterminal_water": {"serialname": "isub_vtw", "mp3": True},
            "vterminal_rain": {"serialname": "isub_vtr", "mp3": True},
            "vterminal_ice": {"serialname": "isub_vti", "mp3": True},
            "vterminal_snow": {"serialname": "isub_vts", "mp3": True},
            "vterminal_graupel": {"serialname": "isub_vtg", "mp3": True},
            "h_var": {"serialname": "isub_h_var", "mp3": True},
            "cvm": {"serialname": "isub_cvm", "mp3": True},
            "te8": {"serialname": "isub_te8", "mp3": True},
            "lcpk": {"serialname": "isub_lcpk", "mp3": True},
            "icpk": {"serialname": "isub_icpk", "mp3": True},
            "tcpk": {"serialname": "isub_tcpk", "mp3": True},
            "tcp3": {"serialname": "isub_tcp3", "mp3": True},
            "di": {"serialname": "isub_di", "mp3": True},
        }

        self.in_vars["parameters"] = [
            "dt",
        ]

        self.out_vars = {
            "qvapor": {"serialname": "isub_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "isub_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "isub_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "isub_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "isub_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "isub_qg", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "isub_pt", "kend": namelist.npz, "mp3": True},
            "cvm": {"serialname": "isub_cvm", "kend": namelist.npz, "mp3": True},
            "te8": {"serialname": "isub_te8", "kend": namelist.npz, "mp3": True},
            "lcpk": {"serialname": "isub_lcpk", "kend": namelist.npz, "mp3": True},
            "icpk": {"serialname": "isub_icpk", "kend": namelist.npz, "mp3": True},
            "tcpk": {"serialname": "isub_tcpk", "kend": namelist.npz, "mp3": True},
            "tcp3": {"serialname": "isub_tcp3", "kend": namelist.npz, "mp3": True},
            "di": {"serialname": "isub_di", "kend": namelist.npz, "mp3": True},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        compute_func = IceFunction(
            self.stencil_factory,
            self.config,
            timestep=inputs.pop("dt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
