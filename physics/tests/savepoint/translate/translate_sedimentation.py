import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.sedimentation import Sedimentation
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class TranslateSedimentation(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "sd_qv"},
            "qliquid": {"serialname": "sd_ql"},
            "qrain": {"serialname": "sd_qr"},
            "qice": {"serialname": "sd_qi"},
            "qsnow": {"serialname": "sd_qs"},
            "qgraupel": {"serialname": "sd_qg"},
            "temperature": {"serialname": "sd_pt"},
            "delp": {"serialname": "sd_delp"},
            "delz": {"serialname": "sd_delz"},
            "density": {"serialname": "sd_den"},
            "density_factor": {"serialname": "sd_denfac"},
            "ua": {"serialname": "sd_u"},
            "va": {"serialname": "sd_v"},
            "wa": {"serialname": "sd_w"},
            "column_energy_change": {"serialname": "sd_dte"},
            "preflux_water": {"serialname": "sd_pfw"},
            "preflux_rain": {"serialname": "sd_pfr"},
            "preflux_ice": {"serialname": "sd_pfi"},
            "preflux_snow": {"serialname": "sd_pfs"},
            "preflux_graupel": {"serialname": "sd_pfg"},
            "vterminal_water": {"serialname": "sd_vtw"},
            "vterminal_rain": {"serialname": "sd_vtr"},
            "vterminal_ice": {"serialname": "sd_vti"},
            "vterminal_snow": {"serialname": "sd_vts"},
            "vterminal_graupel": {"serialname": "sd_vtg"},
            "column_water": {"serialname": "sd_w1"},
            "column_rain": {"serialname": "sd_r1"},
            "column_ice": {"serialname": "sd_i1"},
            "column_snow": {"serialname": "sd_s1"},
            "column_graupel": {"serialname": "sd_g1"},
        }

        self.in_vars["parameters"] = ["dt", "convt"]

        self.out_vars = {
            "qvapor": {"serialname": "sd_qv", "kend": namelist.npz - 1},
            "qliquid": {"serialname": "sd_ql", "kend": namelist.npz - 1},
            "qrain": {"serialname": "sd_qr", "kend": namelist.npz - 1},
            "qice": {"serialname": "sd_qi", "kend": namelist.npz - 1},
            "qsnow": {"serialname": "sd_qs", "kend": namelist.npz - 1},
            "qgraupel": {"serialname": "sd_qg", "kend": namelist.npz - 1},
            "temperature": {"serialname": "sd_pt", "kend": namelist.npz - 1},
            "ua": {"serialname": "sd_u", "kend": namelist.npz - 1},
            "va": {"serialname": "sd_v", "kend": namelist.npz - 1},
            "wa": {"serialname": "sd_w", "kend": namelist.npz - 1},
            "preflux_water": {"serialname": "sd_pfw", "kend": namelist.npz - 1},
            "preflux_rain": {"serialname": "sd_pfr", "kend": namelist.npz - 1},
            "preflux_ice": {"serialname": "sd_pfi", "kend": namelist.npz - 1},
            "preflux_snow": {"serialname": "sd_pfs", "kend": namelist.npz - 1},
            "preflux_graupel": {"serialname": "sd_pfg", "kend": namelist.npz - 1},
            "vterminal_water": {"serialname": "sd_vtw", "kend": namelist.npz - 1},
            "vterminal_rain": {"serialname": "sd_vtr", "kend": namelist.npz - 1},
            "vterminal_ice": {"serialname": "sd_vti", "kend": namelist.npz - 1},
            "vterminal_snow": {"serialname": "sd_vts", "kend": namelist.npz - 1},
            "vterminal_graupel": {"serialname": "sd_vtg", "kend": namelist.npz - 1},
            "column_water": {"serialname": "sd_w1"},
            "column_rain": {"serialname": "sd_r1"},
            "column_ice": {"serialname": "sd_i1"},
            "column_snow": {"serialname": "sd_s1"},
            "column_graupel": {"serialname": "sd_g1"},
            "column_energy_change": {"serialname": "sd_dte", "kend": namelist.npz - 1},
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

        compute_func = Sedimentation(
            self.stencil_factory,
            self.quantity_factory,
            self.config,
            timestep=inputs.pop("dt"),
            convert_mm_day=inputs.pop("convt"),
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
