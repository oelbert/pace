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
            "qvapor": {"serialname": "mpf_qv"},
            "qliquid": {"serialname": "mpf_ql"},
            "qrain": {"serialname": "mpf_qr"},
            "qice": {"serialname": "mpf_qi"},
            "qsnow": {"serialname": "mpf_qs"},
            "qgraupel": {"serialname": "mpf_qg"},
            "qa": {"serialname": "mpf_qa"},
            "ua": {"serialname": "mpf_u"},
            "va": {"serialname": "mpf_v"},
            "wa": {"serialname": "mpf_w"},
            "temperature": {"serialname": "mpf_pt"},
            "delp": {"serialname": "mpf_delp"},
            "delz": {"serialname": "mpf_delz"},
            "density": {"serialname": "mpf_den"},
            "density_factor": {"serialname": "mpf_denfac"},
            "cloud_condensation_nuclei": {"serialname": "mpf_ccn"},
            "cloud_ice_nuclei": {"serialname": "mpf_cin"},
            "preflux_water": {"serialname": "mpf_pfw"},
            "preflux_rain": {"serialname": "mpf_pfr"},
            "preflux_ice": {"serialname": "mpf_pfi"},
            "preflux_snow": {"serialname": "mpf_pfs"},
            "preflux_graupel": {"serialname": "mpf_pfg"},
            "h_var": {"serialname": "mpf_h_var"},
            "rh_adj": {"serialname": "mpf_rh_adj"},
            "column_energy_change": {"serialname": "mpf_dte"},
            "surface_water": {"serialname": "mpf_water"},
            "surface_rain": {"serialname": "mpf_rain"},
            "surface_ice": {"serialname": "mpf_ice"},
            "surface_snow": {"serialname": "mpf_snow"},
            "surface_graupel": {"serialname": "mpf_graupel"},
            "condensation": {"serialname": "mpf_cond"},
            "deposition": {"serialname": "mpf_dep"},
            "evaporation": {"serialname": "mpf_evap"},
            "sublimation": {"serialname": "mpf_sub"},
        }

        self.in_vars["parameters"] = [
            "convt",
            "dt"
        ]

        self.out_vars = {
            "qvapor": {"serialname": "mpf_qv", "kend": namelist.npz - 1},
            "qliquid": {"serialname": "mpf_ql", "kend": namelist.npz - 1},
            "qrain": {"serialname": "mpf_qr", "kend": namelist.npz - 1},
            "qice": {"serialname": "mpf_qi", "kend": namelist.npz - 1},
            "qsnow": {"serialname": "mpf_qs", "kend": namelist.npz - 1},
            "qgraupel": {"serialname": "mpf_qg", "kend": namelist.npz - 1},
            "qa": {"serialname": "mpf_qa", "kend": namelist.npz - 1},
            "ua": {"serialname": "mpf_u", "kend": namelist.npz - 1},
            "va": {"serialname": "mpf_v", "kend": namelist.npz - 1},
            "wa": {"serialname": "mpf_w", "kend": namelist.npz - 1},
            "temperature": {"serialname": "mpf_pt", "kend": namelist.npz - 1},
            "delp": {"serialname": "mpf_delp", "kend": namelist.npz - 1},
            "delz": {"serialname": "mpf_delz", "kend": namelist.npz - 1},
            "density": {"serialname": "mpf_den", "kend": namelist.npz - 1},
            "density_factor": {"serialname": "mpf_denfac", "kend": namelist.npz - 1},
            "cloud_condensation_nuclei": {"serialname": "mpf_ccn", "kend": namelist.npz - 1},
            "cloud_ice_nuclei": {"serialname": "mpf_cin", "kend": namelist.npz - 1},
            "preflux_water": {"serialname": "mpf_pfw", "kend": namelist.npz - 1},
            "preflux_rain": {"serialname": "mpf_pfr", "kend": namelist.npz - 1},
            "preflux_ice": {"serialname": "mpf_pfi", "kend": namelist.npz - 1},
            "preflux_snow": {"serialname": "mpf_pfs", "kend": namelist.npz - 1},
            "preflux_graupel": {"serialname": "mpf_pfg", "kend": namelist.npz - 1},
            "h_var": {"serialname": "mpf_h_var"},
            "rh_adj": {"serialname": "mpf_rh_adj"},
            "column_energy_change": {"serialname": "mpf_dte"},
            "surface_water": {"serialname": "mpf_water"},
            "surface_rain": {"serialname": "mpf_rain"},
            "surface_ice": {"serialname": "mpf_ice"},
            "surface_snow": {"serialname": "mpf_snow"},
            "surface_graupel": {"serialname": "mpf_graupel"},
            "condensation": {"serialname": "mpf_cond"},
            "deposition": {"serialname": "mpf_dep"},
            "evaporation": {"serialname": "mpf_evap"},
            "sublimation": {"serialname": "mpf_sub"},
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
