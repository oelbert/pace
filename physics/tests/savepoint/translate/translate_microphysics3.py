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
            "qvapor": {"serialname": "mp_qv", "microph": True},
            "qliquid": {"serialname": "mp_ql", "microph": True},
            "qrain": {"serialname": "mp_qr", "microph": True},
            "qice": {"serialname": "mp_qi", "microph": True},
            "qsnow": {"serialname": "mp_qs", "microph": True},
            "qgraupel": {"serialname": "mp_qg", "microph": True},
            "qcld": {"serialname": "mp_qa", "microph": True},
            "qcloud_cond_nuclei": {"serialname": "mp_qnl", "microph": True},
            "qcloud_ice_nuclei": {"serialname": "mp_qni", "microph": True},
            "pt": {"serialname": "mp_pt", "microph": True},
            "ua": {"serialname": "mp_ua", "microph": True},
            "va": {"serialname": "mp_va", "microph": True},
            "wa": {"serialname": "mp_wa", "microph": True},
            "delz": {"serialname": "mp_delz", "microph": True},
            "delp": {"serialname": "mp_delp", "microph": True},
            "gsize": {"serialname": "mp_gsize", "microph": True},
            "geopotential_surface_height": {"serialname": "mp_hs", "microph": True},
            "column_water": {"serialname": "mp_water", "microph": True},
            "column_rain": {"serialname": "mp_rain", "microph": True},
            "column_ice": {"serialname": "mp_ice", "microph": True},
            "column_snow": {"serialname": "mp_snow", "microph": True},
            "column_graupel": {"serialname": "mp_graupel", "microph": True},
            "qcon": {"serialname": "mp_q_con", "microph": True},
            "cappa": {"serialname": "mp_cappa", "microph": True},
            "total_energy": {"serialname": "mp_te", "microph": True},
            "preflux_water": {"serialname": "mp_prefluxw", "microph": True},
            "preflux_rain": {"serialname": "mp_prefluxr", "microph": True},
            "preflux_ice": {"serialname": "mp_prefluxi", "microph": True},
            "preflux_snow": {"serialname": "mp_prefluxs", "microph": True},
            "preflux_graupel": {"serialname": "mp_prefluxg", "microph": True},
            "condensation": {"serialname": "mp_cond", "microph": True},
            "deposition": {"serialname": "mp_dep", "microph": True},
            "evaporation": {"serialname": "mp_reevap", "microph": True},
            "sublimation": {"serialname": "mp_sub", "microph": True},
        }
        self.in_vars["parameters"] = [
            "timestep",
            "consv_te",
            "last_step",
        ]

        self.out_vars = {
            "qvapor": {"serialname": "mp_qv", "kend": namelist.npz - 1},
            "qliquid": {"serialname": "mp_ql", "kend": namelist.npz - 1},
            "qrain": {"serialname": "mp_qr", "kend": namelist.npz - 1},
            "qice": {"serialname": "mp_qi", "kend": namelist.npz - 1},
            "qsnow": {"serialname": "mp_qs", "kend": namelist.npz - 1},
            "qgraupel": {"serialname": "mp_qg", "kend": namelist.npz - 1},
            "qcld": {"serialname": "mp_qa", "kend": namelist.npz - 1},
            "pt": {"serialname": "mp_pt", "kend": namelist.npz - 1},
            "ua": {"serialname": "mp_ua", "kend": namelist.npz - 1},
            "va": {"serialname": "mp_va", "kend": namelist.npz - 1},
            "wa": {"serialname": "mp_wa", "kend": namelist.npz - 1},
            "delz": {"serialname": "mp_delz", "kend": namelist.npz - 1},
            "delp": {"serialname": "mp_delp", "kend": namelist.npz - 1},
            "column_water": {"serialname": "mp_water"},
            "column_rain": {"serialname": "mp_rain"},
            "column_ice": {"serialname": "mp_ice"},
            "column_snow": {"serialname": "mp_snow"},
            "column_graupel": {"serialname": "mp_graupel"},
            "qcon": {"serialname": "mp_q_con", "kend": namelist.npz - 1},
            "cappa": {"serialname": "mp_cappa", "kend": namelist.npz - 1},
            "adj_vmr": {"serialname": "mp_adj_vmr", "kend": namelist.npz - 1},
            "total_energy": {"serialname": "mp_te", "kend": namelist.npz - 1},
            "column_energy_change": {"serialname": "mp_dte"},
            "preflux_water": {"serialname": "mp_prefluxw", "kend": namelist.npz - 1},
            "preflux_rain": {"serialname": "mp_prefluxr", "kend": namelist.npz - 1},
            "preflux_ice": {"serialname": "mp_prefluxi", "kend": namelist.npz - 1},
            "preflux_snow": {"serialname": "mp_prefluxs", "kend": namelist.npz - 1},
            "preflux_graupel": {"serialname": "mp_prefluxg", "kend": namelist.npz - 1},
            "condensation": {"serialname": "mp_cond"},
            "deposition": {"serialname": "mp_dep"},
            "evaporation": {"serialname": "mp_reevap"},
            "sublimation": {"serialname": "mp_sub"},
            "particle_concentration_w": {
                "serialname": "mp_pcw",
                "kend": namelist.npz - 1,
            },
            "effective_diameter_w": {"serialname": "mp_edw", "kend": namelist.npz - 1},
            "optical_extinction_w": {"serialname": "mp_oew", "kend": namelist.npz - 1},
            "radar_reflectivity_w": {"serialname": "mp_rrw", "kend": namelist.npz - 1},
            "terminal_velocity_w": {"serialname": "mp_tvw", "kend": namelist.npz - 1},
            "particle_concentration_r": {
                "serialname": "mp_pcr",
                "kend": namelist.npz - 1,
            },
            "effective_diameter_r": {"serialname": "mp_edr", "kend": namelist.npz - 1},
            "optical_extinction_r": {"serialname": "mp_oer", "kend": namelist.npz - 1},
            "radar_reflectivity_r": {"serialname": "mp_rrr", "kend": namelist.npz - 1},
            "terminal_velocity_r": {"serialname": "mp_tvr", "kend": namelist.npz - 1},
            "particle_concentration_i": {
                "serialname": "mp_pci",
                "kend": namelist.npz - 1,
            },
            "effective_diameter_i": {"serialname": "mp_edi", "kend": namelist.npz - 1},
            "optical_extinction_i": {"serialname": "mp_oei", "kend": namelist.npz - 1},
            "radar_reflectivity_i": {"serialname": "mp_rri", "kend": namelist.npz - 1},
            "terminal_velocity_i": {"serialname": "mp_tvi", "kend": namelist.npz - 1},
            "particle_concentration_s": {
                "serialname": "mp_pcs",
                "kend": namelist.npz - 1,
            },
            "effective_diameter_s": {"serialname": "mp_eds", "kend": namelist.npz - 1},
            "optical_extinction_s": {"serialname": "mp_oes", "kend": namelist.npz - 1},
            "radar_reflectivity_s": {"serialname": "mp_rrs", "kend": namelist.npz - 1},
            "terminal_velocity_s": {"serialname": "mp_tvs", "kend": namelist.npz - 1},
            "particle_concentration_g": {
                "serialname": "mp_pcg",
                "kend": namelist.npz - 1,
            },
            "effective_diameter_g": {"serialname": "mp_edg", "kend": namelist.npz - 1},
            "optical_extinction_g": {"serialname": "mp_oeg", "kend": namelist.npz - 1},
            "radar_reflectivity_g": {"serialname": "mp_rrg", "kend": namelist.npz - 1},
            "terminal_velocity_g": {"serialname": "mp_tvg", "kend": namelist.npz - 1},
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

        microphysics_state = MicrophysicsState.init_from_storages(
            inputs,
            sizer=self.sizer,
            quantity_factory=self.quantity_factory,
        )

        microphysics = Microphysics(
            self.stencil_factory,
            self.quantity_factory,
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

        return self.slice_output(inputs)
