import copy  # noqa

import numpy as np  # noqa

import pace.dsl
import pace.dsl.gt4py_utils as utils  # noqa
import pace.util
from pace.physics.stencils.microphysics_v3.microphysics_v3 import Microphysics  # noqa
from pace.physics.stencils.physics import PhysicsState  # noqa
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class TranslateMicrophysics3(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "mp_qv"},
            "qliquid": {"serialname": "mp_ql"},
            "qrain": {"serialname": "mp_qr"},
            "qice": {"serialname": "mp_qi"},
            "qsnow": {"serialname": "mp_qs"},
            "qgraupel": {"serialname": "mp_qg"},
            "qcld": {"serialname": "mp_qa"},
            "qcloud_cond_nuclei": {"serialname": "mp_qnl"},
            "qcloud_ice_nuclei": {"serialname": "mp_qni"},
            "pt": {"serialname": "mp_pt"},
            "ua": {"serialname": "mp_ua"},
            "va": {"serialname": "mp_va"},
            "wa": {"serialname": "mp_wa"},
            "delz": {"serialname": "mp_delz"},
            "delp": {"serialname": "mp_delp"},
            "gsize": {"serialname": "mp_gsize"},
            "geopotential_surface_height": {"serialname": "mp_hs"},
            "column_water": {"serialname": "mp_water"},
            "column_rain": {"serialname": "mp_rain"},
            "column_ice": {"serialname": "mp_ice"},
            "column_snow": {"serialname": "mp_snow"},
            "column_graupel": {"serialname": "mp_graupel"},
            "q_cond": {"serialname": "mp_q_con"},
            "cappa": {"serialname": "mp_cappa"},
            "total_energy": {"serialname": "mp_te"},
            "preflux_water": {"serialname": "mp_prefluxw"},
            "preflux_rain": {"serialname": "mp_prefluxr"},
            "preflux_ice": {"serialname": "mp_prefluxi"},
            "preflux_snow": {"serialname": "mp_prefluxs"},
            "preflux_graupel": {"serialname": "mp_prefluxg"},
            "condensation": {"serialname": "mp_cond"},
            "deposition": {"serialname": "mp_dep"},
            "evaporation": {"serialname": "mp_reevap"},
            "sublimation": {"serialname": "mp_sub"},
        }
        self.in_vars["parameters"] = {
            "timestep": {"serialname": "mp_dt"},
            "hydrostatic": False,
            "consv_te": False,
            "last_step": True,
            "do_inline_mp": False,
        }

        self.out_vars = {
            "qvapor": {"serialname": "mp_qv"},
            "qliquid": {"serialname": "mp_ql"},
            "qrain": {"serialname": "mp_qr"},
            "qice": {"serialname": "mp_qi"},
            "qsnow": {"serialname": "mp_qs"},
            "qgraupel": {"serialname": "mp_qg"},
            "qcld": {"serialname": "mp_qa"},
            "pt": {"serialname": "mp_pt"},
            "ua": {"serialname": "mp_ua"},
            "va": {"serialname": "mp_va"},
            "wa": {"serialname": "mp_wa"},
            "delz": {"serialname": "mp_delz"},
            "delp": {"serialname": "mp_delp"},
            "column_water": {"serialname": "mp_water"},
            "column_rain": {"serialname": "mp_rain"},
            "column_ice": {"serialname": "mp_ice"},
            "column_snow": {"serialname": "mp_snow"},
            "column_graupel": {"serialname": "mp_graupel"},
            "q_cond": {"serialname": "mp_q_con"},
            "cappa": {"serialname": "mp_cappa"},
            "total_energy": {"serialname": "mp_te"},
            "column_energy_change": {"serialname": "mp_dte"},
            "preflux_water": {"serialname": "mp_prefluxw"},
            "preflux_rain": {"serialname": "mp_prefluxr"},
            "preflux_ice": {"serialname": "mp_prefluxi"},
            "preflux_snow": {"serialname": "mp_prefluxs"},
            "preflux_graupel": {"serialname": "mp_prefluxg"},
            "condensation": {"serialname": "mp_cond"},
            "deposition": {"serialname": "mp_dep"},
            "evaporation": {"serialname": "mp_reevap"},
            "sublimation": {"serialname": "mp_sub"},
            "adj_vmr": {"serialname": "mp_adj_vmr"},
            "particle_concentration_w": {"serialname": "mp_pcw"},
            "effective_diameter_w": {"serialname": "mp_edw"},
            "optical_extinction_w": {"serialname": "mp_oew"},
            "radar_reflectivity_w": {"serialname": "mp_rrw"},
            "terminal_velocity_w": {"serialname": "mp_tvw"},
            "particle_concentration_r": {"serialname": "mp_pcr"},
            "effective_diameter_r": {"serialname": "mp_edr"},
            "optical_extinction_r": {"serialname": "mp_oer"},
            "radar_reflectivity_r": {"serialname": "mp_rrr"},
            "terminal_velocity_r": {"serialname": "mp_tvr"},
            "particle_concentration_i": {"serialname": "mp_pci"},
            "effective_diameter_i": {"serialname": "mp_edi"},
            "optical_extinction_i": {"serialname": "mp_oei"},
            "radar_reflectivity_i": {"serialname": "mp_rri"},
            "terminal_velocity_i": {"serialname": "mp_tvi"},
            "particle_concentration_s": {"serialname": "mp_pcs"},
            "effective_diameter_s": {"serialname": "mp_eds"},
            "optical_extinction_s": {"serialname": "mp_oes"},
            "radar_reflectivity_s": {"serialname": "mp_rrs"},
            "terminal_velocity_s": {"serialname": "mp_tvs"},
            "particle_concentration_g": {"serialname": "mp_pcg"},
            "effective_diameter_g": {"serialname": "mp_edg"},
            "optical_extinction_g": {"serialname": "mp_oeg"},
            "radar_reflectivity_g": {"serialname": "mp_rrg"},
            "terminal_velocity_g": {"serialname": "mp_tvg"},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)
