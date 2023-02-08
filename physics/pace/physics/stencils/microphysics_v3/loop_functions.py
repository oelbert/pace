import pace.util.constants as constants


QCMIN = 1.0e-15


def moist_heat_capacity(qvapor, qliquid, qrain, qice, qsnow, qgraupel):
    pass


def column_sedi_melt(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    cvm,
    temperature,
    delz,
    delp,
    ze,
    zt,
    zs,
    timestep,
    v_terminal,
    r1,
    tau_mlt,
    icpk,
    ks,
    ke,
    mode,
):
    li00 = -5
    if mode == "ice":
        q_melt = qice
    elif mode == "snow":
        q_melt = qsnow
    elif mode == "graupel":
        q_melt = qgraupel
    else:
        raise ValueError(f"sedi_melt mode {mode} not ice, snow, or graupel")
    for k in range(ke, ks - 1, -1):
        if v_terminal[k] >= 1.0e-10:
            continue
        if q_melt > constants.QCMIN:
            for m in range(k + 1, ke + 1):
                if zt[k + 1] >= ze[m]:
                    break
                if (zt[k] < ze[m + 1]) and (temperature[m] > constants.TICE):
                    cvm[k] = moist_heat_capacity(
                        qvapor[k], qliquid[k], qrain[k], qice[k], qsnow[k], qgraupel[k]
                    )
                    cvm[m] = moist_heat_capacity(
                        qvapor[m], qliquid[m], qrain[m], qice[m], qsnow[m], qgraupel[m]
                    )
                    dtime = min(timestep, (ze[m] - ze[m + 1]) / v_terminal[k])
                    dtime = min(1.0, dtime / tau_mlt)
                    sink = min(
                        q_melt[k] * delp[k] / delp[m],
                        dtime * (temperature[m] - constants.TICE) / icpk[m],
                    )
                    q_melt -= sink * delp[m] / delp[k]
                    if zt[k] < zs:
                        r1 += sink * delp[m]
                    else:
                        qrain[m] += sink
                    if mode == "ice":
                        qice = q_melt
                    elif mode == "snow":
                        qsnow = q_melt
                    else:
                        qgraupel = q_melt
                    temperature[k] = (
                        temperature[k] * cvm[k] - li00 * sink * delp[m] / delp[k]
                    ) / moist_heat_capacity(
                        qvapor[k], qliquid[k], qrain[k], qice[k], qsnow[k], qgraupel[k]
                    )
                    temperature[m] = (temperature[m] * cvm[m]) / moist_heat_capacity(
                        qvapor[m], qliquid[m], qrain[m], qice[m], qsnow[m], qgraupel[m]
                    )
                if q_melt[k] < QCMIN:
                    break

    return q_melt, qrain, r1, temperature, cvm
