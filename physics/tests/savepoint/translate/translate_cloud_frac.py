import pace.dsl
import pace.util
from pace.physics._config import PhysicsConfig
from pace.physics.stencils.microphysics_v3.cloud_fraction import CloudFraction
from pace.stencils.testing.translate_physics import TranslatePhysicsFortranData2Py


class TranslateCloudFrac(TranslatePhysicsFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: pace.util.Namelist,
        stencil_factory: pace.dsl.StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.in_vars["data_vars"] = {
            "qvapor": {"serialname": "cf_qv", "mp3": True},
            "qliquid": {"serialname": "cf_ql", "mp3": True},
            "qrain": {"serialname": "cf_qr", "mp3": True},
            "qice": {"serialname": "cf_qi", "mp3": True},
            "qsnow": {"serialname": "cf_qs", "mp3": True},
            "qgraupel": {"serialname": "cf_qg", "mp3": True},
            "qa": {"serialname": "cf_qa", "mp3": True},
            "temperature": {"serialname": "cf_pt", "mp3": True},
            "density": {"serialname": "cf_den", "mp3": True},
            "pz": {"serialname": "cf_pz", "mp3": True},
            "h_var": {"serialname": "cf_h_var", "mp3": True},
            "gsize": {"serialname": "cf_gsize", "mp3": True},
        }

        self.out_vars = {
            "qa": {"serialname": "ne_qa", "kend": namelist.npz, "mp3": True},
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

        compute_func = CloudFraction(
            self.stencil_factory,
            self.quantity_factory,
            self.config,
        )

        compute_func(**inputs)

        return self.slice_output(inputs)
