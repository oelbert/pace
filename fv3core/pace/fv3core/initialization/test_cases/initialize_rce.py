import math

import numpy as np

import pace.dsl.gt4py_utils as utils
import pace.fv3core.initialization.init_utils as init_utils
import pace.util as fv3util
import pace.util.constants as constants
from pace.fv3core.dycore_state import DycoreState
from pace.util.grid import GridData, great_circle_distance_lon_lat, lon_lat_midpoint

# Based on initialization procedure for RCEMIP simulations
# See https://gmd.copernicus.org/articles/11/793/2018/

# Surface pressure (Pa)
p0 = 1014.8e2
# Surface temperature (K)
T0 = 300.0
# Surface specific humidity (kg/kg)
q0 = 18.65e-3
# Stratospheric specific humidity (kg/kg)
qt = 1e-14
# Tropopause height (m)
zt = 15e3 
# Specific humidity scale heights (m)
zq1 = 4e3 
zq2 = 7.5e3
# Virtual temperature lapse rate (K/m)
Gamma = 0.0067

def _compute_rce_thermodynamic_profiles(p):

    Rd = constants.RDGAS
    g = constants.GRAV
    zvir = constants.ZVIR

    # Compute surface virtual temperature
    Tv0 = T0*(1 + zvir*q0)

    # Compute pressure and virtual temperature at tropopause
    Tvt = Tv0 - Gamma*zt
    pt = p0*(Tvt/Tv0)**(g/(Rd*Gamma))

    # Invert (5) to get height as a function of pressure 
    z = np.zeros_like(p)
    trop = p >= pt
    strat = p < pt
    z[trop] = Tv0*(1 - (p[trop]/p0)**(Rd*Gamma/g))/Gamma
    z[strat] = zt - Rd*Tvt/g*np.log(p[strat]/pt)

    # Compute temperature and specific humidity
    Tv = np.zeros_like(p)
    qv = np.zeros_like(p)
    Tv[trop] = Tv0 - Gamma*z[trop]
    Tv[strat] = Tvt 
    qv[trop] = q0*np.exp(-z[trop]/zq1)*np.exp(-(z[trop]/zq2)**2)
    qv[strat] = qt
    T = Tv/(1 + zvir*qv)

    return T, qv

def init_rce_state(
    grid_data: GridData,
    quantity_factory: fv3util.QuantityFactory,
    hydrostatic: bool,
    comm: fv3util.TileCommunicator,
) -> DycoreState:
    """
    Create a DycoreState object with quantities initialized following
    the RCEMIP protocol for 300 K SST.
    """

    sample_quantity = grid_data.lat
    shape = (*sample_quantity.data.shape[:2], grid_data.ak.data.shape[0])
    numpy_state = init_utils.empty_numpy_dycore_state(shape)

    # Initialize vertical grid
    delp = np.zeros(shape)
    ps = np.full(shape[:2], p0)
    ak = grid_data.ak.data
    bk = grid_data.bk.data
    ptop = ak[0] + p0*bk[0]
    delp[:,:,:-1] = init_utils.initialize_delp(ps, ak, bk)
    pe = init_utils.initialize_edge_pressure(delp, ptop)
    peln = np.log(pe)
    pk, pkz = init_utils.initialize_kappa_pressures(pe, peln, ptop)

    # Initialize temperature and specific humidity
    pt = np.zeros(shape)
    qvapor = np.zeros(shape)
    p = 0.5*(ak[1:] + ak[:-1] + p0*(bk[1:] + bk[:-1]))
    pt_profile, qvapor_profile = _compute_rce_thermodynamic_profiles(p)
    pt[:,:,:-1] = pt_profile[np.newaxis,np.newaxis,:]
    qvapor[:,:,:-1] = qvapor_profile[np.newaxis,np.newaxis,:]

    # Initialize layer thicknesses
    delz = np.zeros(shape)
    print(pe.shape)
    print(delz.shape)
    delz[:, :, :-1] = (
        constants.RDGAS
        * pt[:, :, :-1]
        * (1 + constants.ZVIR * qvapor[:, :, :-1])
        / constants.GRAV
        * np.log(pe[:, :, :-1] / pe[:, :, 1:])
    )

    # Break symmetry by adding thermal noise to lowest five levels 
    # TODO: best way to add noise?
    # np.random.seed(comm.rank)
    # for amp, k in zip([0.1, 0.08, 0.06, 0.04, 0.02], [-2, -3, -4, -5, -6]):
    #     pt[:,:,k] += amp*np.random.rand(pt.shape[0], pt.shape[1])

    # Set state
    # Note that fields not explicitly set are 0
    numpy_state.delp[:] = delp
    numpy_state.delz[:] = delz
    numpy_state.pe[:] = pe
    numpy_state.peln[:] = peln
    numpy_state.pk[:] = pk
    numpy_state.pkz[:] = pkz
    numpy_state.ps[:] = pe[:, :, -1]
    numpy_state.pt[:] = pt
    numpy_state.qvapor[:] = qvapor
    state = DycoreState.init_from_numpy_arrays(
        numpy_state.__dict__,
        sizer=quantity_factory.sizer,
        backend=sample_quantity.metadata.gt4py_backend,
    )

    return state