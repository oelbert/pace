import dataclasses
from typing import Tuple

import f90nml


DEFAULT_INT = 0
DEFAULT_STR = ""
DEFAULT_FLOAT = 0.0
DEFAULT_BOOL = False
# Global set of namelist defaults, attached to class for namespacing and static typing
class NamelistDefaults:
    layout = (1, 1)
    grid_type = 0
    do_f3d = False
    inline_q = False
    do_skeb = False  # save dissipation estimate
    use_logp = False
    moist_phys = True
    check_negative = False
    # gfdl_cloud_mucrophys.F90
    tau_r2g = 900.0  # rain freezing during fast_sat
    tau_smlt = 900.0  # snow melting timescale
    tau_gmlt = 600.0  # snow melting timescale
    tau_g2r = 600.0  # graupel melting to rain
    tau_imlt = 600.0  # cloud ice melting
    tau_i2s = 1000.0  # cloud ice to snow auto - conversion
    tau_l2r = 900.0  # cloud water to rain auto - conversion
    tau_g2v = 1200.0  # graupel sublimation
    tau_v2g = 21600.0  # graupel deposition -- make it a slow process
    tau_wbf = 300.0  # wegener bergeron findeisen timescale
    sat_adj0 = 0.90  # adjustment factor (0: no, 1: full) during fast_sat_adj
    ql_gen = 1.0e-3  # max new cloud water during remapping step if fast_sat_adj = .t.
    ql_mlt = 2.0e-3  # max value of cloud water allowed from melted cloud ice
    qs_mlt = 1.0e-6  # max cloud water due to snow melt
    ql0_max = 2.0e-3  # max cloud water value (auto converted to rain)
    t_sub = 184.0  # min temp for sublimation of cloud ice
    qi_gen = 1.82e-6  # max cloud ice generation during remapping step
    qi_lim = 1.0  # cloud ice limiter to prevent large ice build up
    qi0_max = 1.0e-4  # max cloud ice value (by other sources)
    rad_snow = True  # consider snow in cloud fraciton calculation
    rad_rain = True  # consider rain in cloud fraction calculation
    rad_graupel = True  # consider graupel in cloud fraction calculation
    tintqs = False  # use temperature in the saturation mixing in PDF
    dw_ocean = 0.10  # base value for ocean
    dw_land = 0.20  # base value for subgrid deviation / variability over land
    # cloud scheme 0 - ?
    # 1: old fvgfs gfdl) mp implementation
    # 2: binary cloud scheme (0 / 1)
    icloud_f = 0
    cld_min = 0.05  # !< minimum cloud fraction
    tau_l2v = 300.0  # cloud water to water vapor (evaporation)
    tau_v2l = 150.0  # water vapor to cloud water (condensation)
    c2l_ord = 4
    regional = False
    m_split = 0
    convert_ke = False
    breed_vortex_inline = False
    use_old_omega = True
    use_logp = False
    rf_fast = False
    p_ref = 1e5  # Surface pressure used to construct a horizontally-uniform reference
    adiabatic = False
    nf_omega = 1
    fv_sg_adj = -1
    n_sponge = 1
    fast_sat_adj = True
    qc_crt = 5.0e-8  # Minimum condensate mixing ratio to allow partial cloudiness
    c_cracw = 0.8  # Rain accretion efficiency
    c_paut = (
        0.5  # Autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
    )
    c_pgacs = 0.01  # Snow to graupel "accretion" eff. (was 0.1 in zetac)
    c_psaci = 0.05  # Accretion: cloud ice to snow (was 0.1 in zetac)
    ccn_l = 270.0  # CCN over land (cm^-3)
    ccn_o = 90.0  # CCN over ocean (cm^-3)
    use_rhc_cevap = False  # cap of rh for cloud water evaporation
    const_vg = False  # Fall velocity tuning constant of graupel
    const_vi = False  # Fall velocity tuning constant of ice
    const_vr = False  # Fall velocity tuning constant of rain water
    const_vs = False  # Fall velocity tuning constant of snow
    is_fac = 0.2  # Cloud ice sublimation temperature factor
    rh_fac = 10.0  # cloud water condensation / evaporation relative humidity factor
    sed_fac = 1.0  # coefficient for sedimentation fall,
    # Scale from 1.0 (implicit) to 0.0 (lagrangian)
    rhc_cevap = 0.85  # maximum relative humidity for cloud water evaporation
    vw_fac = 1.0
    vi_fac = 1.0  # if const_vi: 1/3
    vs_fac = 1.0  # if const_vs: 1.
    vg_fac = 1.0  # if const_vg: 2.
    vr_fac = 1.0  # if const_vr: 4.
    de_ice = False  # To prevent excessive build-up of cloud ice from external sources
    do_qa = True  # Do inline cloud fraction
    do_sedi_heat = False  # Transport of heat in sedimentation
    do_sedi_melt = True  # Melt cloud ice, snow, and graupel during sedimentation
    do_sedi_uv = True  # Transport of horizontal momentum in sedimentation
    do_sedi_w = True  # Transport of vertical motion in sedimentation
    fix_negative = True  # Fix negative water species
    do_cond_timescale = False  # Whether to apply a timescale to condensation
    do_hail = False  # Use hail parameters instead of graupel
    consv_checker = False  # Turn on energy and water conservation check in microphysics
    do_warm_rain = False  # Do only warm rain microphysics
    do_wbf = False  # Do Wegener Bergeron Findeisen process
    do_psd_water_fall = False  # Calculate cloud water terminal velocity based on PSD
    do_psd_ice_fall = False  # Calculate cloud ice terminal velocity based on PSD
    do_psd_water_num = False  # Calculate cloud water number concentration based on PSD
    do_psd_ice_num = False  # Calculate cloud ice number concentration based on PSD
    irain_f = 0  # Cloud water to rain auto conversion scheme
    mono_prof = False  # Perform terminal fall with mono ppm scheme
    mp_time = 225.0  # Maximum microphysics timestep (sec)
    prog_ccn = False  # Do prognostic ccn (yi ming's method)
    qi0_crt = 1.0e-04  # Cloud ice to snow autoconversion threshold
    qs0_crt = 0.003  # Snow to graupel density threshold (0.6e-3 in purdue lin scheme)
    rh_inc = 0.25  # RH increment for complete evaporation of cloud water and cloud ice
    rh_inr = 0.25  # RH increment for minimum evaporation of rain
    rthresh = 1e-05  # Critical cloud drop radius (micrometers)
    sedi_transport = True  # Transport of momentum in sedimentation
    use_ppm = False  # Use ppm fall scheme
    vw_max = 0.01  # Maximum fall speed for cloud water
    vg_max = 16.0  # Maximum fall speed for graupel
    vi_max = 1.0  # Maximum fall speed for ice
    vr_max = 16.0  # Maximum fall speed for rain
    vs_max = 2.0  # Maximum fall speed for snow
    z_slope_ice = True  # Use linear mono slope for autoconversions
    z_slope_liq = True  # Use linear mono slope for autoconversions
    tice = 273.16  # set tice = 165. to turn off ice - phase phys (kessler emulator)
    alin = 842.0  # "a" in lin1983
    alinw = 3.0e7  # "a" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
    alini = 7.0e2  # "a" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
    alinr = 842.0  # "a" in Lin et al. (1983) for rain (Liu and Orville 1969)
    alins = 4.8  # "a" in Lin et al. (1983) for snow (straka 2009)
    aling = 1.0  # "a" in Lin et al. (1983) for graupel (Pruppacher and Klett 2010)
    alinh = 1.0  # "a" in Lin et al. (1983) for hail (Pruppacher and Klett 2010)
    blinw = 2.0  # "b" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
    blini = 1.0  # "b" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
    blinr = 0.8  # "b" in Lin et al. (1983) for rain (Liu and Orville 1969)
    blins = 0.25  # "b" in Lin et al. (1983) for snow (straka 2009)
    bling = 0.5  # "b" in Lin et al. (1983) for graupel (Pruppacher and Klett 2010)
    blinh = 0.5  # "b" in Lin et al. (1983) for hail (Pruppacher and Klett 2010)
    clin = 4.8  # "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)
    ntimes = 1  # Number of cloud microphysics sub cycles
    do_inline_mp = False  # Whether the microphsyics is called inside of the dycore
    n0w_sig = 1.1  # cwater significand (Lin et al. 1983) (m^-4) (Martin et al. 1994)
    n0i_sig = 1.3  # cice significand (Lin et al. 1983) (m^-4) (McFarquhar et al. 2015)
    n0r_sig = 8.0
    # rain significand (Lin et al. 1983) (m^-4) (Marshall and Palmer 1948)
    n0s_sig = 3.0  # snow significand (Lin et al. 1983) (m^-4) (Gunn and Marshall 1958)
    n0g_sig = 4.0
    # graupel significand (Rutledge and Hobbs 1984) (m^-4) (Houze et al. 1979)
    n0h_sig = 4.0
    # hail significand (Lin et al. 1983) (m^-4) (Federer and Waldvogel 1975)
    n0w_exp = 41  # cwater exponent (Lin et al. 1983) (m^-4) (Martin et al. 1994)
    n0i_exp = 18  # cice exponent (Lin et al. 1983) (m^-4) (McFarquhar et al. 2015)
    n0r_exp = 6  # rain exponent (Lin et al. 1983) (m^-4) (Marshall and Palmer 1948)
    n0s_exp = 6  # snow exponent (Lin et al. 1983) (m^-4) (Gunn and Marshall 1958)
    n0g_exp = 6
    # graupel exponent (Rutledge and Hobbs 1984) (m^-4) (Houze et al. 1979)
    n0h_exp = 4  # hail exponent (Lin et al. 1983) (m^-4) (Federer and Waldvogel 1975)
    muw = (
        6.0  # shape parameter of cloud water in Gamma distribution (Martin et al. 1994)
    )
    mui = 3.35  # Gamma shape parameter of cloud ice (McFarquhar et al. 2015)
    mur = 1.0
    # shape parameter of rain in Gamma distribution (Marshall and Palmer 1948)
    mus = 1.0  # shape parameter of snow in Gamma distribution (Gunn and Marshall 1958)
    mug = 1.0  # shape parameter of graupel in Gamma distribution (Houze et al. 1979)
    muh = 1.0
    # shape parameter of hail in Gamma distribution (Federer and Waldvogel 1975)
    inflag = 1  # Ice nucleation scheme:
    # 1: Hong et al. (2004)
    # 2: Meyers et al. (1992)
    # 3: Meyers et al. (1992)
    # 4: Cooper (1986)
    # 5: Fletcher (1962)
    igflag = 3  # Ice generation scheme
    # 1: WSM6
    # 2: WSM6 with 0 at 0 C
    # 3: WSM6 with 0 at 0 C and fixed value at - 10 C
    # 4: combination of 1 and 3
    ifflag = 1  # Ice fall scheme
    # 1: Deng and Mace (2008)
    # 2: Heymsfield and Donner (1990)
    sedflag = 1  # sedimentation scheme
    # 1: implicit scheme
    # 2: explicit scheme
    # 3: lagrangian scheme
    # 4: combined implicit and lagrangian scheme

    @classmethod
    def as_dict(cls):
        return {
            name: default
            for name, default in cls.__dict__.items()
            if not name.startswith("_")
        }


@dataclasses.dataclass
class Namelist:
    # data_set: Any
    # date_out_of_range: str
    # do_sst_pert: bool
    # interp_oi_sst: bool
    # no_anom_sst: bool
    # sst_pert: float
    # sst_pert_type: str
    # use_daily: bool
    # use_ncep_ice: bool
    # use_ncep_sst: bool
    # blocksize: int
    # chksum_debug: bool
    """
    note: dycore_only may not be used in this model
     the same way it is in the Fortran version, watch for
     consequences of these inconsistencies, or more closely
     parallel the Fortran structure
    """
    dycore_only: bool = DEFAULT_BOOL
    # fdiag: float
    # knob_ugwp_azdir: Tuple[int, int, int, int]
    # knob_ugwp_doaxyz: int
    # knob_ugwp_doheat: int
    # knob_ugwp_dokdis: int
    # knob_ugwp_effac: Tuple[int, int, int, int]
    # knob_ugwp_ndx4lh: int
    # knob_ugwp_solver: int
    # knob_ugwp_source: Tuple[int, int, int, int]
    # knob_ugwp_stoch: Tuple[int, int, int, int]
    # knob_ugwp_version: int
    # knob_ugwp_wvspec: Tuple[int, int, int, int]
    # launch_level: int
    # reiflag: int
    # reimax: float
    # reimin: float
    # rewmax: float
    # rewmin: float
    # atmos_nthreads: int
    # calendar: Any
    # current_date: Any
    days: int = 0
    dt_atmos: int = DEFAULT_INT
    # dt_ocean: Any
    hours: int = 0
    # memuse_verbose: Any
    minutes: int = 0
    # months: Any
    # ncores_per_node: Any
    seconds: int = 0
    # use_hyper_thread: Any
    # max_axes: Any
    # max_files: Any
    # max_num_axis_sets: Any
    # prepend_date: Any
    # checker_tr: Any
    # filtered_terrain: Any
    # gfs_dwinds: Any
    # levp: Any
    # nt_checker: Any
    # checksum_required: Any
    # max_files_r: Any
    # max_files_w: Any
    # clock_grain: Any
    # domains_stack_size: Any
    # print_memory_usage: Any
    a_imp: float = DEFAULT_FLOAT
    # adjust_dry_mass: Any
    beta: float = DEFAULT_FLOAT
    # consv_am: Any
    consv_te: float = DEFAULT_FLOAT
    d2_bg: float = DEFAULT_FLOAT
    d2_bg_k1: float = DEFAULT_FLOAT
    d2_bg_k2: float = DEFAULT_FLOAT
    d4_bg: float = DEFAULT_FLOAT
    d_con: float = DEFAULT_FLOAT
    d_ext: float = DEFAULT_FLOAT
    dddmp: float = DEFAULT_FLOAT
    delt_max: float = DEFAULT_FLOAT
    # dnats: int
    do_sat_adj: bool = DEFAULT_BOOL
    do_vort_damp: bool = DEFAULT_BOOL
    # dwind_2d: Any
    # external_ic: Any
    fill: bool = DEFAULT_BOOL
    # fill_dp: bool
    # fv_debug: Any
    # gfs_phil: Any
    hord_dp: int = DEFAULT_INT
    hord_mt: int = DEFAULT_INT
    hord_tm: int = DEFAULT_INT
    hord_tr: int = DEFAULT_INT
    hord_vt: int = DEFAULT_INT
    hydrostatic: bool = DEFAULT_BOOL
    # io_layout: Any
    k_split: int = DEFAULT_INT
    ke_bg: float = DEFAULT_FLOAT
    kord_mt: int = DEFAULT_INT
    kord_tm: int = DEFAULT_INT
    kord_tr: int = DEFAULT_INT
    kord_wz: int = DEFAULT_INT
    layout: Tuple[int, int] = (1, 1)
    # make_nh: bool
    # mountain: bool
    n_split: int = DEFAULT_INT
    # na_init: Any
    # ncep_ic: Any
    # nggps_ic: Any
    nord: int = DEFAULT_INT
    npx: int = DEFAULT_INT
    npy: int = DEFAULT_INT
    npz: int = DEFAULT_INT
    ntiles: int = DEFAULT_INT
    # nudge: Any
    # nudge_qv: Any
    nwat: int = DEFAULT_INT
    p_fac: float = DEFAULT_FLOAT
    # phys_hydrostatic: Any
    # print_freq: Any
    # range_warn: Any
    # reset_eta: Any
    rf_cutoff: float = DEFAULT_FLOAT
    tau: float = DEFAULT_FLOAT
    # tau_h2o: Any
    # use_hydro_pressure: Any
    vtdm4: float = DEFAULT_FLOAT
    # warm_start: bool
    z_tracer: bool = DEFAULT_BOOL
    c_cracw: float = NamelistDefaults.c_cracw
    c_paut: float = NamelistDefaults.c_paut
    c_pgacs: float = NamelistDefaults.c_pgacs
    c_psaci: float = NamelistDefaults.c_psaci
    ccn_l: float = NamelistDefaults.ccn_l
    ccn_o: float = NamelistDefaults.ccn_o
    const_vg: bool = NamelistDefaults.const_vg
    const_vi: bool = NamelistDefaults.const_vi
    const_vr: bool = NamelistDefaults.const_vr
    const_vs: bool = NamelistDefaults.const_vs
    qc_crt: float = NamelistDefaults.qc_crt
    vw_fac: float = NamelistDefaults.vw_fac
    vs_fac: float = NamelistDefaults.vs_fac
    vg_fac: float = NamelistDefaults.vg_fac
    vi_fac: float = NamelistDefaults.vi_fac
    vr_fac: float = NamelistDefaults.vr_fac
    de_ice: bool = NamelistDefaults.de_ice
    do_qa: bool = NamelistDefaults.do_qa
    do_sedi_heat: bool = NamelistDefaults.do_sedi_heat
    do_sedi_melt: bool = NamelistDefaults.do_sedi_melt
    do_sedi_uv: bool = NamelistDefaults.do_sedi_uv
    do_sedi_w: bool = NamelistDefaults.do_sedi_w
    fast_sat_adj: bool = NamelistDefaults.fast_sat_adj
    fix_negative: bool = NamelistDefaults.fix_negative
    do_cond_timescale: bool = NamelistDefaults.do_cond_timescale
    do_hail: bool = NamelistDefaults.do_hail
    consv_checker: bool = NamelistDefaults.consv_checker
    do_warm_rain: bool = NamelistDefaults.do_warm_rain
    do_wbf: bool = NamelistDefaults.do_wbf
    do_psd_water_fall: bool = NamelistDefaults.do_psd_water_fall
    do_psd_ice_fall: bool = NamelistDefaults.do_psd_ice_fall
    do_psd_water_num: bool = NamelistDefaults.do_psd_water_num
    do_psd_ice_num: bool = NamelistDefaults.do_psd_ice_num
    irain_f: int = NamelistDefaults.irain_f
    mono_prof: bool = NamelistDefaults.mono_prof
    mp_time: float = NamelistDefaults.mp_time
    prog_ccn: bool = NamelistDefaults.prog_ccn
    qi0_crt: float = NamelistDefaults.qi0_crt
    qs0_crt: float = NamelistDefaults.qs0_crt
    rhc_cevap: float = NamelistDefaults.rhc_cevap
    is_fac: float = NamelistDefaults.is_fac
    rh_fac: float = NamelistDefaults.rh_fac
    sed_fac: float = NamelistDefaults.sed_fac
    rh_inc: float = NamelistDefaults.rh_inc
    rh_inr: float = NamelistDefaults.rh_inr
    # rh_ins: Any
    rthresh: float = NamelistDefaults.rthresh
    sedi_transport: bool = NamelistDefaults.sedi_transport
    # use_ccn: Any
    use_ppm: bool = NamelistDefaults.use_ppm
    vw_max: float = NamelistDefaults.vw_max
    vg_max: float = NamelistDefaults.vg_max
    vi_max: float = NamelistDefaults.vi_max
    vr_max: float = NamelistDefaults.vr_max
    vs_max: float = NamelistDefaults.vs_max
    z_slope_ice: bool = NamelistDefaults.z_slope_ice
    z_slope_liq: bool = NamelistDefaults.z_slope_liq
    tice: float = NamelistDefaults.tice
    alin: float = NamelistDefaults.alin
    alinw: float = NamelistDefaults.alinw
    alini: float = NamelistDefaults.alini
    alinr: float = NamelistDefaults.alinr
    alins: float = NamelistDefaults.alins
    aling: float = NamelistDefaults.aling
    alinh: float = NamelistDefaults.alinh
    blinw: float = NamelistDefaults.blinw
    blini: float = NamelistDefaults.blini
    blinr: float = NamelistDefaults.blinr
    blins: float = NamelistDefaults.blins
    bling: float = NamelistDefaults.bling
    blinh: float = NamelistDefaults.blinh
    clin: float = NamelistDefaults.clin
    ntimes: int = NamelistDefaults.ntimes
    do_inline_mp: bool = NamelistDefaults.do_inline_mp
    use_rhc_cevap: bool = NamelistDefaults.use_rhc_cevap
    # c0s_shal: Any
    # c1_shal: Any
    # cal_pre: Any
    # cdmbgwd: Any
    # cnvcld: Any
    # cnvgwd: Any
    # debug: Any
    # do_deep: Any
    # dspheat: Any
    # fhcyc: Any
    # fhlwr: Any
    # fhswr: Any
    # fhzero: Any
    # hybedmf: Any
    # iaer: Any
    # ialb: Any
    # ico2: Any
    # iems: Any
    # imfdeepcnv: Any
    # imfshalcnv: Any
    # imp_physics: Any
    # isol: Any
    # isot: Any
    # isubc_lw: Any
    # isubc_sw: Any
    # ivegsrc: Any
    # ldiag3d: Any
    # lwhtr: Any
    # ncld: int
    # nst_anl: Any
    # pdfcld: Any
    # pre_rad: Any
    # prslrd0: Any
    # random_clds: Any
    # redrag: Any
    # satmedmf: Any
    # shal_cnv: Any
    # swhtr: Any
    # trans_trac: Any
    # use_ufo: Any
    # xkzm_h: Any
    # xkzm_m: Any
    # xkzminv: Any
    # interp_method: Any
    # lat_s: Any
    # lon_s: Any
    # ntrunc: Any
    # fabsl: Any
    # faisl: Any
    # faiss: Any
    # fnabsc: Any
    # fnacna: Any
    # fnaisc: Any
    # fnalbc: Any
    # fnalbc2: Any
    # fnglac: Any
    # fnmskh: Any
    # fnmxic: Any
    # fnslpc: Any
    # fnsmcc: Any
    # fnsnoa: Any
    # fnsnoc: Any
    # fnsotc: Any
    # fntg3c: Any
    # fntsfa: Any
    # fntsfc: Any
    # fnvegc: Any
    # fnvetc: Any
    # fnvmnc: Any
    # fnvmxc: Any
    # fnzorc: Any
    # fsicl: Any
    # fsics: Any
    # fslpl: Any
    # fsmcl: Any
    # fsnol: Any
    # fsnos: Any
    # fsotl: Any
    # ftsfl: Any
    # ftsfs: Any
    # fvetl: Any
    # fvmnl: Any
    # fvmxl: Any
    # ldebug: Any
    grid_type: int = NamelistDefaults.grid_type
    do_f3d: bool = NamelistDefaults.do_f3d
    inline_q: bool = NamelistDefaults.inline_q
    do_skeb: bool = NamelistDefaults.do_skeb  # save dissipation estimate
    use_logp: bool = NamelistDefaults.use_logp
    moist_phys: bool = NamelistDefaults.moist_phys
    check_negative: bool = NamelistDefaults.check_negative
    # gfdl_cloud_microphys.F90
    tau_r2g: float = NamelistDefaults.tau_r2g  # rain freezing during fast_sat
    tau_smlt: float = NamelistDefaults.tau_smlt  # snow melting
    tau_gmlt: float = NamelistDefaults.tau_gmlt  # snow melting
    tau_g2r: float = NamelistDefaults.tau_g2r  # graupel melting to rain
    tau_imlt: float = NamelistDefaults.tau_imlt  # cloud ice melting
    tau_i2s: float = NamelistDefaults.tau_i2s  # cloud ice to snow auto - conversion
    tau_l2r: float = NamelistDefaults.tau_l2r  # cloud water to rain auto - conversion
    tau_g2v: float = NamelistDefaults.tau_g2v  # graupel sublimation
    tau_v2g: float = (
        NamelistDefaults.tau_v2g
    )  # graupel deposition -- make it a slow process
    tau_wbf: float = NamelistDefaults.tau_wbf
    sat_adj0: float = (
        NamelistDefaults.sat_adj0
    )  # adjustment factor (0: no 1: full) during fast_sat_adj
    ql_gen: float = (
        1.0e-3  # max new cloud water during remapping step if fast_sat_adj = .t.
    )
    ql_mlt: float = (
        NamelistDefaults.ql_mlt
    )  # max value of cloud water allowed from melted cloud ice
    qs_mlt: float = NamelistDefaults.qs_mlt  # max cloud water due to snow melt
    ql0_max: float = (
        NamelistDefaults.ql0_max
    )  # max cloud water value (auto converted to rain)
    t_sub: float = NamelistDefaults.t_sub  # min temp for sublimation of cloud ice
    qi_gen: float = (
        NamelistDefaults.qi_gen
    )  # max cloud ice generation during remapping step
    qi_lim: float = (
        NamelistDefaults.qi_lim
    )  # cloud ice limiter to prevent large ice build up
    qi0_max: float = NamelistDefaults.qi0_max  # max cloud ice value (by other sources)
    rad_snow: bool = (
        NamelistDefaults.rad_snow
    )  # consider snow in cloud fraction calculation
    rad_rain: bool = (
        NamelistDefaults.rad_rain
    )  # consider rain in cloud fraction calculation
    rad_graupel: bool = (
        NamelistDefaults.rad_graupel
    )  # consider graupel in cloud fraction calculation
    tintqs: bool = (
        NamelistDefaults.tintqs
    )  # use temperature in the saturation mixing in PDF
    dw_ocean: float = NamelistDefaults.dw_ocean  # base value for ocean
    dw_land: float = (
        NamelistDefaults.dw_land
    )  # base value for subgrid deviation / variability over land
    # cloud scheme 0 - ?
    # 1: old fvgfs gfdl) mp implementation
    # 2: binary cloud scheme (0 / 1)
    icloud_f: int = NamelistDefaults.icloud_f
    cld_min: float = NamelistDefaults.cld_min  # !< minimum cloud fraction
    tau_l2v: float = (
        NamelistDefaults.tau_l2v
    )  # cloud water to water vapor (evaporation)
    tau_v2l: float = (
        NamelistDefaults.tau_v2l
    )  # water vapor to cloud water (condensation)
    n0w_sig: float = NamelistDefaults.n0w_sig
    n0i_sig: float = NamelistDefaults.n0i_sig
    n0r_sig: float = NamelistDefaults.n0r_sig
    n0s_sig: float = NamelistDefaults.n0s_sig
    n0g_sig: float = NamelistDefaults.n0g_sig
    n0h_sig: float = NamelistDefaults.n0h_sig
    n0w_exp: float = NamelistDefaults.n0w_exp
    n0i_exp: float = NamelistDefaults.n0i_exp
    n0r_exp: float = NamelistDefaults.n0r_exp
    n0s_exp: float = NamelistDefaults.n0s_exp
    n0g_exp: float = NamelistDefaults.n0g_exp
    n0h_exp: float = NamelistDefaults.n0h_exp
    muw: float = NamelistDefaults.muw
    mui: float = NamelistDefaults.mui
    mur: float = NamelistDefaults.mur
    mus: float = NamelistDefaults.mus
    muh: float = NamelistDefaults.muh
    mug: float = NamelistDefaults.mug
    inflag: int = NamelistDefaults.inflag
    igflag: int = NamelistDefaults.igflag
    ifflag: int = NamelistDefaults.ifflag
    sedflag: int = NamelistDefaults.sedflag
    c2l_ord: int = NamelistDefaults.c2l_ord
    regional: bool = NamelistDefaults.regional
    m_split: int = NamelistDefaults.m_split
    convert_ke: bool = NamelistDefaults.convert_ke
    breed_vortex_inline: bool = NamelistDefaults.breed_vortex_inline
    use_old_omega: bool = NamelistDefaults.use_old_omega
    rf_fast: bool = NamelistDefaults.rf_fast
    adiabatic: bool = NamelistDefaults.adiabatic
    nf_omega: int = NamelistDefaults.nf_omega
    fv_sg_adj: int = NamelistDefaults.fv_sg_adj
    n_sponge: int = NamelistDefaults.n_sponge

    @classmethod
    def from_f90nml(cls, namelist: f90nml.Namelist):
        namelist_dict = namelist_to_flatish_dict(namelist.items())
        namelist_dict = {
            key: value
            for key, value in namelist_dict.items()
            if key in cls.__dataclass_fields__  # type: ignore
        }
        return cls(**namelist_dict)


def namelist_to_flatish_dict(nml_input):
    nml = dict(nml_input)
    for name, value in nml.items():
        if isinstance(value, f90nml.Namelist):
            nml[name] = namelist_to_flatish_dict(value)
    flatter_namelist = {}
    for key, value in nml.items():
        if isinstance(value, dict):
            for subkey, subvalue in value.items():
                if subkey in flatter_namelist:
                    raise ValueError(
                        "Cannot flatten this namelist, duplicate keys: " + subkey
                    )
                flatter_namelist[subkey] = subvalue
        else:
            flatter_namelist[key] = value
    return flatter_namelist
