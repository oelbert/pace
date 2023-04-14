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
            "qvapor": {"serialname": "cf_qv", "microph": True},
            "qliquid": {"serialname": "cf_ql", "microph": True},
            "qrain": {"serialname": "cf_qr", "microph": True},
            "qice": {"serialname": "cf_qi", "microph": True},
            "qsnow": {"serialname": "cf_qs", "microph": True},
            "qgraupel": {"serialname": "cf_qg", "microph": True},
            "qa": {"serialname": "cf_qa", "microph": True},
            "temperature": {"serialname": "cf_pt", "microph": True},
            "density": {"serialname": "cf_den", "microph": True},
            "pz": {"serialname": "cf_pz", "microph": True},
            "h_var": {"serialname": "cf_h_var", "microph": True},
            "gsize": {"serialname": "cf_gsize", "microph": True},
        }

        self.out_vars = {
            "qa": {"serialname": "ne_qa", "kend": namelist.npz - 1},
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
