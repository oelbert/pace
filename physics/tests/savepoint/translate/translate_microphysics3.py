import copy  # noqa

import numpy as np  # noqa

import pace.dsl
import pace.dsl.gt4py_utils as utils
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.microphysics_state import MicrophysicsState
from pace.physics.stencils.microphysics_v3.microphysics_v3 import Microphysics
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
        self.in_vars["parameters"] = [
            "timestep",
            "consv_te",
            "last_step",
        ]

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
            "adj_vmr": {"serialname": "mp_adj_vmr"},
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
        pconf = PhysicsConfig.from_namelist(namelist)
        self.config = pconf.microphysics

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)

        storage = utils.make_storage_from_shape(
            self.grid_indexing.domain_full(add=(1, 1, 1)),
            origin=self.grid_indexing.origin_compute(),
            backend=self.stencil_factory.backend,
        )

        sizer = pace.util.SubtileGridSizer.from_tile_params(
            nx_tile=self.namelist.npx - 1,
            ny_tile=self.namelist.npy - 1,
            nz=self.namelist.npz,
            n_halo=3,
            extra_dim_lengths={},
            layout=self.namelist.layout,
        )

        quantity_factory = pace.util.QuantityFactory.from_backend(
            sizer, self.stencil_factory.backend
        )

        microphysics_state = MicrophysicsState.init_from_storages(
            inputs,
            sizer=sizer,
            quantity_factory=quantity_factory,
        )

        microphysics = Microphysics(
            self.stencil_factory,
            quantity_factory,
            self.grid.grid_data,
            self.config,
            consv_te=inputs["consv_te"],
        )

        microphysics(
            microphysics_state,
            timestep=inputs["timestep"],
            last_step=inputs["last_step"],
        )

        # copy microphysics state back to inputs
        inputs["column_water"] = microphysics_state.column_water
        inputs["column_rain"] = microphysics_state.column_rain
        inputs["column_ice"] = microphysics_state.column_ice
        inputs["column_snow"] = microphysics_state.column_snow
        inputs["column_graupel"] = microphysics_state.column_graupel
        inputs["condensation"] = microphysics_state.condensation
        inputs["deposition"] = microphysics_state.deposition
        inputs["sublimation"] = microphysics_state.sublimation
        inputs["evaporation"] = microphysics_state.evaporation
        inputs["column_energy_change"] = microphysics_state.column_energy_change
        inputs["adj_vmr"] = microphysics_state.adj_vmr
        inputs["particle_concentration_w"] = microphysics_state.particle_concentration_w
        inputs["effective_diameter_w"] = microphysics_state.effective_diameter_w
        inputs["optical_extinction_w"] = microphysics_state.optical_extinction_w
        inputs["radar_reflectivity_w"] = microphysics_state.radar_reflectivity_w
        inputs["terminal_velocity_w"] = microphysics_state.terminal_velocity_w
        inputs["particle_concentration_r"] = microphysics_state.particle_concentration_r
        inputs["effective_diameter_r"] = microphysics_state.effective_diameter_r
        inputs["optical_extinction_r"] = microphysics_state.optical_extinction_r
        inputs["radar_reflectivity_r"] = microphysics_state.radar_reflectivity_r
        inputs["terminal_velocity_r"] = microphysics_state.terminal_velocity_r
        inputs["particle_concentration_i"] = microphysics_state.particle_concentration_i
        inputs["effective_diameter_i"] = microphysics_state.effective_diameter_i
        inputs["optical_extinction_i"] = microphysics_state.optical_extinction_i
        inputs["radar_reflectivity_i"] = microphysics_state.radar_reflectivity_i
        inputs["terminal_velocity_i"] = microphysics_state.terminal_velocity_i
        inputs["particle_concentration_s"] = microphysics_state.particle_concentration_s
        inputs["effective_diameter_s"] = microphysics_state.effective_diameter_s
        inputs["optical_extinction_s"] = microphysics_state.optical_extinction_s
        inputs["radar_reflectivity_s"] = microphysics_state.radar_reflectivity_s
        inputs["terminal_velocity_s"] = microphysics_state.terminal_velocity_s
        inputs["particle_concentration_g"] = microphysics_state.particle_concentration_g
        inputs["effective_diameter_g"] = microphysics_state.effective_diameter_g
        inputs["optical_extinction_g"] = microphysics_state.optical_extinction_g
        inputs["radar_reflectivity_g"] = microphysics_state.radar_reflectivity_g
        inputs["terminal_velocity_g"] = microphysics_state.terminal_velocity_g

        out = self.slice_output(inputs)
