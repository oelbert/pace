import pace.util


class MicrophysicsState:
    """
    A state that contains everything that goes into or comes out of microphysics
    """

    def __init__(
        self,
        qvapor: pace.util.Quantity,
        qliquid: pace.util.Quantity,
        qrain: pace.util.Quantity,
        qice: pace.util.Quantity,
        qsnow: pace.util.Quantity,
        qgraupel: pace.util.Quantity,
        qcld: pace.util.Quantity,
        qcloud_cond_nuclei: pace.util.Quantity,
        qcloud_ice_nuclei: pace.util.Quantity,
        ua: pace.util.Quantity,
        va: pace.util.Quantity,
        wa: pace.util.Quantity,
        delp: pace.util.Quantity,
        delz: pace.util.Quantity,
        pt: pace.util.Quantity,
        geopotential_surface_height: pace.util.Quantity,
        preflux_water: pace.util.Quantity,
        preflux_rain: pace.util.Quantity,
        preflux_ice: pace.util.Quantity,
        preflux_snow: pace.util.Quantity,
        preflux_graupel: pace.util.Quantity,
        column_water: pace.util.Quantity,
        column_rain: pace.util.Quantity,
        column_ice: pace.util.Quantity,
        column_snow: pace.util.Quantity,
        column_graupel: pace.util.Quantity,
        condensation: pace.util.Quantity,
        deposition: pace.util.Quantity,
        evaporation: pace.util.Quantity,
        sublimation: pace.util.Quantity,
        total_energy: pace.util.Quantity,
        column_energy_change: pace.util.Quantity,
        adj_vmr: pace.util.Quantity,
        particle_concentration_w: pace.util.Quantity,
        effective_diameter_w: pace.util.Quantity,
        optical_extinction_w: pace.util.Quantity,
        radar_reflectivity_w: pace.util.Quantity,
        terminal_velocity_w: pace.util.Quantity,
        particle_concentration_r: pace.util.Quantity,
        effective_diameter_r: pace.util.Quantity,
        optical_extinction_r: pace.util.Quantity,
        radar_reflectivity_r: pace.util.Quantity,
        terminal_velocity_r: pace.util.Quantity,
        particle_concentration_i: pace.util.Quantity,
        effective_diameter_i: pace.util.Quantity,
        optical_extinction_i: pace.util.Quantity,
        radar_reflectivity_i: pace.util.Quantity,
        terminal_velocity_i: pace.util.Quantity,
        particle_concentration_s: pace.util.Quantity,
        effective_diameter_s: pace.util.Quantity,
        optical_extinction_s: pace.util.Quantity,
        radar_reflectivity_s: pace.util.Quantity,
        terminal_velocity_s: pace.util.Quantity,
        particle_concentration_g: pace.util.Quantity,
        effective_diameter_g: pace.util.Quantity,
        optical_extinction_g: pace.util.Quantity,
        radar_reflectivity_g: pace.util.Quantity,
        terminal_velocity_g: pace.util.Quantity,
    ):
        self.qvapor = qvapor
        self.qliquid = qliquid
        self.qrain = qrain
        self.qice = qice
        self.qsnow = qsnow
        self.qgraupel = qgraupel
        self.qcld = qcld
        self.qcloud_cond_nuclei = qcloud_cond_nuclei
        self.qcloud_ice_nuclei = qcloud_ice_nuclei
        self.ua = ua
        self.va = va
        self.wa = wa
        self.delp = delp
        self.delz = delz
        self.pt = pt
        self.geopotential_surface_height = geopotential_surface_height
        self.preflux_water = preflux_water
        self.preflux_rain = preflux_rain
        self.preflux_ice = preflux_ice
        self.preflux_snow = preflux_snow
        self.preflux_graupel = preflux_graupel
        self.column_water = column_water
        self.column_rain = column_rain
        self.column_ice = column_ice
        self.column_snow = column_snow
        self.column_graupel = column_graupel
        self.condensation = condensation
        self.deposition = deposition
        self.evaporation = evaporation
        self.sublimation = sublimation
        self.total_energy = total_energy
        self.column_energy_change = column_energy_change
        self.adj_vmr = adj_vmr
        self.particle_concentration_w = particle_concentration_w
        self.effective_diameter_w = effective_diameter_w
        self.optical_extinction_w = optical_extinction_w
        self.radar_reflectivity_w = radar_reflectivity_w
        self.terminal_velocity_w = terminal_velocity_w
        self.particle_concentration_r = particle_concentration_r
        self.effective_diameter_r = effective_diameter_r
        self.optical_extinction_r = optical_extinction_r
        self.radar_reflectivity_r = radar_reflectivity_r
        self.terminal_velocity_r = terminal_velocity_r
        self.particle_concentration_i = particle_concentration_i
        self.effective_diameter_i = effective_diameter_i
        self.optical_extinction_i = optical_extinction_i
        self.radar_reflectivity_i = radar_reflectivity_i
        self.terminal_velocity_i = terminal_velocity_i
        self.particle_concentration_s = particle_concentration_s
        self.effective_diameter_s = effective_diameter_s
        self.optical_extinction_s = optical_extinction_s
        self.radar_reflectivity_s = radar_reflectivity_s
        self.terminal_velocity_s = terminal_velocity_s
        self.particle_concentration_g = particle_concentration_g
        self.effective_diameter_g = effective_diameter_g
        self.optical_extinction_g = optical_extinction_g
        self.radar_reflectivity_g = radar_reflectivity_g
        self.terminal_velocity_g = terminal_velocity_g
