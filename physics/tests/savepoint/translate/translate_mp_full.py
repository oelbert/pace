import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.mp_full import FullMicrophysics
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class TranslateMPFull(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "mpf_qv", "mp3": True},
            "qliquid": {"serialname": "mpf_ql", "mp3": True},
            "qrain": {"serialname": "mpf_qr", "mp3": True},
            "qice": {"serialname": "mpf_qi", "mp3": True},
            "qsnow": {"serialname": "mpf_qs", "mp3": True},
            "qgraupel": {"serialname": "mpf_qg", "mp3": True},
            "qa": {"serialname": "mpf_qa", "mp3": True},
            "ua": {"serialname": "mpf_u", "mp3": True},
            "va": {"serialname": "mpf_v", "mp3": True},
            "wa": {"serialname": "mpf_w", "mp3": True},
            "temperature": {"serialname": "mpf_pt", "mp3": True},
            "delp": {"serialname": "mpf_delp", "mp3": True},
            "delz": {"serialname": "mpf_delz", "mp3": True},
            "density": {"serialname": "mpf_den", "mp3": True},
            "density_factor": {"serialname": "mpf_denfac", "mp3": True},
            "cloud_condensation_nuclei": {"serialname": "mpf_ccn", "mp3": True},
            "cloud_ice_nuclei": {"serialname": "mpf_cin", "mp3": True},
            "preflux_water": {"serialname": "mpf_pfw", "mp3": True},
            "preflux_rain": {"serialname": "mpf_pfr", "mp3": True},
            "preflux_ice": {"serialname": "mpf_pfi", "mp3": True},
            "preflux_snow": {"serialname": "mpf_pfs", "mp3": True},
            "preflux_graupel": {"serialname": "mpf_pfg", "mp3": True},
            "h_var": {"serialname": "mpf_h_var", "mp3": True},
            "rh_adj": {"serialname": "mpf_rh_adj", "mp3": True},
            "column_energy_change": {"serialname": "mpf_dte", "mp3": True},
            "surface_water": {"serialname": "mpf_water", "mp3": True},
            "surface_rain": {"serialname": "mpf_rain", "mp3": True},
            "surface_ice": {"serialname": "mpf_ice", "mp3": True},
            "surface_snow": {"serialname": "mpf_snow", "mp3": True},
            "surface_graupel": {"serialname": "mpf_graupel", "mp3": True},
            "condensation": {"serialname": "mpf_cond", "mp3": True},
            "deposition": {"serialname": "mpf_dep", "mp3": True},
            "evaporation": {"serialname": "mpf_evap", "mp3": True},
            "sublimation": {"serialname": "mpf_sub", "mp3": True},
        }

        self.in_vars["parameters"] = ["convt", "dt"]

        self.out_vars = {
            "qvapor": {"serialname": "mpf_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "mpf_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "mpf_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "mpf_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "mpf_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "mpf_qg", "kend": namelist.npz, "mp3": True},
            "qa": {"serialname": "mpf_qa", "kend": namelist.npz, "mp3": True},
            "ua": {"serialname": "mpf_u", "kend": namelist.npz, "mp3": True},
            "va": {"serialname": "mpf_v", "kend": namelist.npz, "mp3": True},
            "wa": {"serialname": "mpf_w", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "mpf_pt", "kend": namelist.npz, "mp3": True},
            "delp": {"serialname": "mpf_delp", "kend": namelist.npz, "mp3": True},
            "delz": {"serialname": "mpf_delz", "kend": namelist.npz, "mp3": True},
            "density": {"serialname": "mpf_den", "kend": namelist.npz, "mp3": True},
            "density_factor": {
                "serialname": "mpf_denfac",
                "kend": namelist.npz,
                "mp3": True,
            },
            "cloud_condensation_nuclei": {
                "serialname": "mpf_ccn",
                "kend": namelist.npz,
                "mp3": True,
            },
            "cloud_ice_nuclei": {
                "serialname": "mpf_cin",
                "kend": namelist.npz,
                "mp3": True,
            },
            "preflux_water": {
                "serialname": "mpf_pfw",
                "kend": namelist.npz,
                "mp3": True,
            },
            "preflux_rain": {
                "serialname": "mpf_pfr",
                "kend": namelist.npz,
                "mp3": True,
            },
            "preflux_ice": {"serialname": "mpf_pfi", "kend": namelist.npz, "mp3": True},
            "preflux_snow": {
                "serialname": "mpf_pfs",
                "kend": namelist.npz,
                "mp3": True,
            },
            "preflux_graupel": {
                "serialname": "mpf_pfg",
                "kend": namelist.npz,
                "mp3": True,
            },
            "h_var": {"serialname": "mpf_h_var", "mp3": True},
            "rh_adj": {"serialname": "mpf_rh_adj", "mp3": True},
            "column_energy_change": {"serialname": "mpf_dte", "mp3": True},
            "surface_water": {"serialname": "mpf_water", "mp3": True},
            "surface_rain": {"serialname": "mpf_rain", "mp3": True},
            "surface_ice": {"serialname": "mpf_ice", "mp3": True},
            "surface_snow": {"serialname": "mpf_snow", "mp3": True},
            "surface_graupel": {"serialname": "mpf_graupel", "mp3": True},
            "condensation": {"serialname": "mpf_cond", "mp3": True},
            "deposition": {"serialname": "mpf_dep", "mp3": True},
            "evaporation": {"serialname": "mpf_evap", "mp3": True},
            "sublimation": {"serialname": "mpf_sub", "mp3": True},
        }

        self.stencil_factory = stencil_factory
        self.grid_indexing = self.stencil_factory.grid_indexing
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

        compute_func = FullMicrophysics(
            self.stencil_factory,
            self.quantity_factory,
            self.config,
            inputs.pop("dt"),
            self.config.ntimes,
            inputs.pop("convt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
