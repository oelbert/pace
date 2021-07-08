#!/usr/bin/env python3

import os
import sys
import numpy as np

SERIALBOX_DIR = "/usr/local/serialbox"
sys.path.append(SERIALBOX_DIR + "/python")
import serialbox as ser
import physics_driver


# Serialized Variables from GFSPhysicsDriver-In Savepoint
IN_VARS_GFSPD = ["IPD_area", "IPD_gq0", "IPD_gt0", "IPD_gu0", "IPD_gv0", "IPD_kdt", 
                 "IPD_levs", "IPD_lradar", "IPD_ntrac", "IPD_phii", 
                 "IPD_prsi", "IPD_qgrs", "IPD_refl_10cm", "IPD_tgrs", "IPD_xlon"]

# Serialized Variables from GFSPhysicsDriver-Out Savepoint
OUT_VARS_GFSPD = ["IPD_gq0", "IPD_gt0", "IPD_gu0", "IPD_gv0"]

# Serialized Variables for inputs into get_prs_fv3
IN_VARS_PRS = ["prs_ix", "prs_levs", "prs_ntrac", "prs_phii", "prs_prsi", "prs_tgrs", "prs_qgrs", "prs_del", "prs_del_gz"]

# Serialized Variables for outputs from get_prs_fv3
OUT_VARS_PRS = ["prs_del", "prs_del_gz"]

# Serialized Variables for inputs into get_phi_fv3
IN_VARS_PHI = ["phi_del_gz", "phi_gq0", "phi_gt0", "phi_ix", "phi_levs", "phi_ntrac", "phi_phii", "phi_phil"]

# Serialized Variables for outputs from get_phi_fv3
OUT_VARS_PHI = ["phi_del_gz", "phi_phii", "phi_phil"]

SELECT_SP = None
# SELECT_SP = {"tile": 2, "savepoint": "satmedmfvdif-in-iter2-000000"}


def data_dict_from_var_list(var_list, serializer, savepoint):
    d = {}
    for var in var_list:
        data = serializer.read(var, savepoint)
        # convert single element numpy arrays to scalars
        if data.size == 1:
            data = data.item()
        d[var] = data
    return d


def compare_data(exp_data, ref_data):
    assert set(exp_data.keys()) == set(
        ref_data.keys()
    ), "Entries of exp and ref dictionaries don't match"
    for key in ref_data:
        ind = np.array(
            np.nonzero(~np.isclose(exp_data[key], ref_data[key], equal_nan=True))
        )
        if ind.size > 0:
            i = tuple(ind[:, 0])
            print("FAIL at ", key, i, exp_data[key][i], ref_data[key][i])
        # assert np.allclose(exp_data[key], ref_data[key], equal_nan=True), (
        #     "Data does not match for field " + key
        # )


for tile in range(6):

    if SELECT_SP is not None:
        if tile != SELECT_SP["tile"]:
            continue

    serializer = ser.Serializer(
        ser.OpenModeKind.Read, "c12_6ranks_baroclinic_dycore_microphysics", "Generator_rank" + str(tile)
    )

    savepoints = serializer.savepoint_list()

    isready = False
    for sp in savepoints:

        if SELECT_SP is not None:
            if sp.name != SELECT_SP["savepoint"] and sp.name != SELECT_SP[
                "savepoint"
            ].replace("-in-", "-out-"):
                continue

        if sp.name.startswith("GFSPhysicsDriver-IN"):

            if isready:
                raise Exception("out-of-order data enountered: " + sp.name)

            print("> running ", f"tile-{tile}", sp)

            # read serialized input data
            in_data = data_dict_from_var_list(IN_VARS_GFSPD, serializer, sp)

            # run Python version
            out_data_postphi = physics_driver.run(in_data)

            isready = True

        # if sp.name.startswith("PrsFV3-In"):
        #     print("> running ", f"tile-{tile}", sp)

        #     # read serialized input data
        #     ref_data = data_dict_from_var_list(IN_VARS_PRS, serializer, sp)
   
        #     compare_data(out_data_preprs, ref_data)

        # if sp.name.startswith("PrsFV3-Out"):
        #     print("> running ", f"tile-{tile}", sp)

        #     # read serialized input data
        #     ref_data = data_dict_from_var_list(OUT_VARS_PRS, serializer, sp)
            
        #     compare_data(out_data_postprs, ref_data)

        # if sp.name.startswith("PhiFV3-In"):
        #     print("> running ", f"tile-{tile}", sp)

        #     # read serialized input data
        #     ref_data = data_dict_from_var_list(IN_VARS_PHI, serializer, sp)
   
        #     compare_data(out_data_prephi, ref_data)

        if sp.name.startswith("PhiFV3-Out"):
            print("> running ", f"tile-{tile}", sp)

            # read serialized input data
            ref_data = data_dict_from_var_list(OUT_VARS_PHI, serializer, sp)
   
            compare_data(out_data_postphi, ref_data)