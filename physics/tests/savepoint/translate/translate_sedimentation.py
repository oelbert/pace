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
            "qvapor": {"serialname": "sd_qv", "mp3": True},
            "qliquid": {"serialname": "sd_ql", "mp3": True},
            "qrain": {"serialname": "sd_qr", "mp3": True},
            "qice": {"serialname": "sd_qi", "mp3": True},
            "qsnow": {"serialname": "sd_qs", "mp3": True},
            "qgraupel": {"serialname": "sd_qg", "mp3": True},
            "temperature": {"serialname": "sd_pt", "mp3": True},
            "delp": {"serialname": "sd_delp", "mp3": True},
            "delz": {"serialname": "sd_delz", "mp3": True},
            "density": {"serialname": "sd_den", "mp3": True},
            "density_factor": {"serialname": "sd_denfac", "mp3": True},
            "ua": {"serialname": "sd_u", "mp3": True},
            "va": {"serialname": "sd_v", "mp3": True},
            "wa": {"serialname": "sd_w", "mp3": True},
            "column_energy_change": {"serialname": "sd_dte", "mp3": True},
            "preflux_water": {"serialname": "sd_pfw", "mp3": True},
            "preflux_rain": {"serialname": "sd_pfr", "mp3": True},
            "preflux_ice": {"serialname": "sd_pfi", "mp3": True},
            "preflux_snow": {"serialname": "sd_pfs", "mp3": True},
            "preflux_graupel": {"serialname": "sd_pfg", "mp3": True},
            "vterminal_water": {"serialname": "sd_vtw", "mp3": True},
            "vterminal_rain": {"serialname": "sd_vtr", "mp3": True},
            "vterminal_ice": {"serialname": "sd_vti", "mp3": True},
            "vterminal_snow": {"serialname": "sd_vts", "mp3": True},
            "vterminal_graupel": {"serialname": "sd_vtg", "mp3": True},
            "column_water": {"serialname": "sd_w1", "mp3": True},
            "column_rain": {"serialname": "sd_r1", "mp3": True},
            "column_ice": {"serialname": "sd_i1", "mp3": True},
            "column_snow": {"serialname": "sd_s1", "mp3": True},
            "column_graupel": {"serialname": "sd_g1", "mp3": True},
        }

        self.in_vars["parameters"] = ["dt", "convt"]

        self.out_vars = {
            "qvapor": {"serialname": "sd_qv", "kend": namelist.npz, "mp3": True},
            "qliquid": {"serialname": "sd_ql", "kend": namelist.npz, "mp3": True},
            "qrain": {"serialname": "sd_qr", "kend": namelist.npz, "mp3": True},
            "qice": {"serialname": "sd_qi", "kend": namelist.npz, "mp3": True},
            "qsnow": {"serialname": "sd_qs", "kend": namelist.npz, "mp3": True},
            "qgraupel": {"serialname": "sd_qg", "kend": namelist.npz, "mp3": True},
            "temperature": {"serialname": "sd_pt", "kend": namelist.npz, "mp3": True},
            "ua": {"serialname": "sd_u", "kend": namelist.npz, "mp3": True},
            "va": {"serialname": "sd_v", "kend": namelist.npz, "mp3": True},
            "wa": {"serialname": "sd_w", "kend": namelist.npz, "mp3": True},
            "preflux_water": {"serialname": "sd_pfw", "kend": namelist.npz, "mp3": True},
            "preflux_rain": {"serialname": "sd_pfr", "kend": namelist.npz, "mp3": True},
            "preflux_ice": {"serialname": "sd_pfi", "kend": namelist.npz, "mp3": True},
            "preflux_snow": {"serialname": "sd_pfs", "kend": namelist.npz, "mp3": True},
            "preflux_graupel": {"serialname": "sd_pfg", "kend": namelist.npz, "mp3": True},
            "vterminal_water": {"serialname": "sd_vtw", "kend": namelist.npz, "mp3": True},
            "vterminal_rain": {"serialname": "sd_vtr", "kend": namelist.npz, "mp3": True},
            "vterminal_ice": {"serialname": "sd_vti", "kend": namelist.npz, "mp3": True},
            "vterminal_snow": {"serialname": "sd_vts", "kend": namelist.npz, "mp3": True},
            "vterminal_graupel": {"serialname": "sd_vtg", "kend": namelist.npz, "mp3": True},
            "column_water": {"serialname": "sd_w1", "mp3": True},
            "column_rain": {"serialname": "sd_r1", "mp3": True},
            "column_ice": {"serialname": "sd_i1", "mp3": True},
            "column_snow": {"serialname": "sd_s1", "mp3": True},
            "column_graupel": {"serialname": "sd_g1", "mp3": True},
            "column_energy_change": {"serialname": "sd_dte", "kend": namelist.npz, "mp3": True},
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

        for var in inputs.keys():
            if len(inputs[var].shape) == 3:
                inputs[var] = pace.util.Quantity(
                    inputs[var],
                    dims=[pace.util.X_DIM, pace.util.Y_DIM, pace.util.Z_DIM],
                    units = "unknown"
                )
            elif len(inputs[var].shape) == 2:
                inputs[var] = pace.util.Quantity(
                    inputs[var],
                    dims=[pace.util.X_DIM, pace.util.Y_DIM],
                    units = "unknown"
                )
            else:
                raise TypeError(
                    f"input data with strange len: {len(inputs[var].shape)}"
                )

        compute_func(**inputs)

        return self.slice_output(inputs)
