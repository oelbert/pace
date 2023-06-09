import numbers
from typing import Tuple

import numpy as np
import scipy.special

class Field:
    def __init__(self, field, offsets: Tuple[int, ...], dimensions: Tuple[bool, bool, bool]):
        ii = iter(range(3))
        self.idx_to_data = tuple(
            [next(ii) if has_dim else None for has_dim in dimensions]
            + list(range(sum(dimensions), len(field.shape)))
        )

        shape = [field.shape[i] if i is not None else 1 for i in self.idx_to_data]
        self.field_view = np.reshape(field.data, shape).view(np.ndarray)

        self.offsets = offsets

    @classmethod
    def empty(cls, shape, dtype, offset):
        return cls(np.empty(shape, dtype=dtype), offset, (True, True, True))

    def shim_key(self, key):
        new_args = []
        if not isinstance(key, tuple):
            key = (key, )
        for index in self.idx_to_data:
            if index is None:
                new_args.append(slice(None, None))
            else:
                idx = key[index]
                offset = self.offsets[index]
                if isinstance(idx, slice):
                    new_args.append(
                        slice(idx.start + offset, idx.stop + offset, idx.step) if offset else idx
                    )
                else:
                    new_args.append(idx + offset)
        if not isinstance(new_args[2], (numbers.Integral, slice)):
            new_args = self.broadcast_and_clip_variable_k(new_args)
        return tuple(new_args)

    def broadcast_and_clip_variable_k(self, new_args: tuple):
        assert isinstance(new_args[0], slice) and isinstance(new_args[1], slice)
        if np.max(new_args[2]) >= self.field_view.shape[2] or np.min(new_args[2]) < 0:
            new_args[2] = np.clip(new_args[2].copy(), 0, self.field_view.shape[2]-1)
        new_args[:2] = np.broadcast_arrays(
            np.expand_dims(
                np.arange(new_args[0].start, new_args[0].stop),
                axis=tuple(i for i in range(self.field_view.ndim) if i != 0)
            ),
            np.expand_dims(
                np.arange(new_args[1].start, new_args[1].stop),
                axis=tuple(i for i in range(self.field_view.ndim) if i != 1)
            ),
        )
        return new_args

    def __getitem__(self, key):
        return self.field_view.__getitem__(self.shim_key(key))

    def __setitem__(self, key, value):
        return self.field_view.__setitem__(self.shim_key(key), value)


def run(*, qvapor, qliquid, qrain, qice, qsnow, qgraupel, temperature, density, delp, cloud_ice_nuclei, te, cvm, lcpk, icpk, tcpk, tcp3, dep, sub, qsi, dqdt, pidep0, pidep, qi_crt, sink1, sink2, tmp, dq, _domain_, _origin_):

    # --- begin domain boundary shortcuts ---
    _di_, _dj_, _dk_ = 0, 0, 0
    _dI_, _dJ_, _dK_ = _domain_
    # --- end domain padding ---

    qvapor = Field(qvapor, _origin_['qvapor'], (True, True, True))
    qliquid = Field(qliquid, _origin_['qliquid'], (True, True, True))
    qrain = Field(qrain, _origin_['qrain'], (True, True, True))
    qice = Field(qice, _origin_['qice'], (True, True, True))
    qsnow = Field(qsnow, _origin_['qsnow'], (True, True, True))
    qgraupel = Field(qgraupel, _origin_['qgraupel'], (True, True, True))
    temperature = Field(temperature, _origin_['temperature'], (True, True, True))
    density = Field(density, _origin_['density'], (True, True, True))
    delp = Field(delp, _origin_['delp'], (True, True, True))
    cloud_ice_nuclei = Field(cloud_ice_nuclei, _origin_['cloud_ice_nuclei'], (True, True, True))
    te = Field(te, _origin_['te'], (True, True, True))
    cvm = Field(cvm, _origin_['cvm'], (True, True, True))
    lcpk = Field(lcpk, _origin_['lcpk'], (True, True, True))
    icpk = Field(icpk, _origin_['icpk'], (True, True, True))
    tcpk = Field(tcpk, _origin_['tcpk'], (True, True, True))
    tcp3 = Field(tcp3, _origin_['tcp3'], (True, True, True))
    dep = Field(dep, _origin_['dep'], (True, True, False))
    sub = Field(sub, _origin_['sub'], (True, True, False))
    qsi = Field(qsi, _origin_['qsi'], (True, True, True))
    dqdt = Field(dqdt, _origin_['dqdt'], (True, True, True))
    pidep0 = Field(pidep0, _origin_['pidep0'], (True, True, True))
    pidep = Field(pidep, _origin_['pidep'], (True, True, True))
    qi_crt = Field(qi_crt, _origin_['qi_crt'], (True, True, True))
    sink1 = Field(sink1, _origin_['sink1'], (True, True, True))
    sink2 = Field(sink2, _origin_['sink2'], (True, True, True))
    tmp = Field(tmp, _origin_['tmp'], (True, True, True))
    dq = Field(dq, _origin_['dq'], (True, True, True))
    
    sz_pidep__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    q_solid__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    pidep__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.int64, (0, 0, 0))
    mask_140070139854032_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), bool, (0, 0, 0))
    qi_crt__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sink__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sz_sink1__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    pcb__444_89_35__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    delta_graupel__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    a__1fb_49_13__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    dqdt__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    tk__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    pca__444_89_35__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    qsnow__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    cloud_ice_nuclei__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    icpk__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    mask_140070134222512_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), bool, (0, 0, 0))
    tcp3__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    q__8f5_62_20__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    mask_140070139851968_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), bool, (0, 0, 0))
    mask_140070136950448_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), bool, (0, 0, 0))
    RETURN_VALUE__1fb_49_13__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    qi_gen__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    qliquid__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    return_val__066_9_8__8f5_62_20__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sz_qsi__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sz_tmp__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    qice__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    icpk__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    diff__1fb_125_35__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    b__1fb_125_35__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    delta_rain__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    delta_vapor__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    cvm__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    lcpk__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    delta_liquid__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    tcpk__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    qsi__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    tmp__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    dq__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    diff__1fb_49_13__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sz_qi_crt__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    tcp3__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    RETURN_VALUE__1fb_125_35__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sz_dqdt__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    RETURN_VALUE__066_9_8__8f5_62_20__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    qvapor__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sz_dq__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sz_pidep0__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    cvm__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    mask_140070138870128_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), bool, (0, 0, 0))
    temp__8f5_62_20__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    mu__444_89_35__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    sz_sink2__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    qrain__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    delta_snow__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    q_l__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    lcpk__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    tcpk__f2f_144_12__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    qgraupel__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    tc__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    dqdt__8f5_62_20__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    temperature__e56_244_16_gen_0 = Field.empty((_dI_ + 0, _dJ_ + 0, _dK_), np.float64, (0, 0, 0))
    

    with np.errstate(divide='ignore', over='ignore', under='ignore', invalid='ignore'):

    
    # --- begin vertical block ---
        k, K = _dk_, _dK_
        for k_ in range(k, K):

            # --- begin horizontal block --
            i, I = _di_ - 0, _dI_ + 0
            j, J = _dj_ - 0, _dJ_ + 0

            qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = qvapor[i:I, j:J, k_:k_ + 1]
            qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = qliquid[i:I, j:J, k_:k_ + 1]
            qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = qrain[i:I, j:J, k_:k_ + 1]
            qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = qice[i:I, j:J, k_:k_ + 1]
            qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = qsnow[i:I, j:J, k_:k_ + 1]
            qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = qgraupel[i:I, j:J, k_:k_ + 1]
            cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = cloud_ice_nuclei[i:I, j:J, k_:k_ + 1]
            temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = temperature[i:I, j:J, k_:k_ + 1]
            cvm__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = cvm[i:I, j:J, k_:k_ + 1]
            lcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = lcpk[i:I, j:J, k_:k_ + 1]
            icpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = icpk[i:I, j:J, k_:k_ + 1]
            tcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = tcpk[i:I, j:J, k_:k_ + 1]
            tcp3__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = tcp3[i:I, j:J, k_:k_ + 1]
            sz_qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = qsi[i:I, j:J, k_:k_ + 1]
            sz_dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = dqdt[i:I, j:J, k_:k_ + 1]
            sz_pidep0__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = pidep0[i:I, j:J, k_:k_ + 1]
            sz_pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = pidep[i:I, j:J, k_:k_ + 1]
            sz_qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = qi_crt[i:I, j:J, k_:k_ + 1]
            sz_sink1__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = sink1[i:I, j:J, k_:k_ + 1]
            sz_sink2__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = sink2[i:I, j:J, k_:k_ + 1]
            sz_tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = tmp[i:I, j:J, k_:k_ + 1]
            sz_dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = dq[i:I, j:J, k_:k_ + 1]
            sz_qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            sz_dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            sz_pidep0__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            sz_pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            sz_qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            sz_sink1__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            sz_sink2__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            sz_tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            sz_dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.float64(0.0)
            mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1] = (temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] < np.float64(273.15))
            pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], np.int64(0), pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], np.maximum((np.float64(273.15) - np.float64(160.0)), np.minimum(temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], (np.float64(273.15) + np.float64(102.0)))), temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            mask_140070134222512_gen_0[i:I, j:J, k_:k_ + 1] = (temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] < np.float64(273.15))
            return_val__066_9_8__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070134222512_gen_0[i:I, j:J, k_:k_ + 1]), (np.float64(611.21) * np.exp((((np.float64(-260.0) * np.log((temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / np.float64(273.15)))) + ((np.float64(2904599.0) * (temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] - np.float64(273.15))) / (temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * np.float64(273.15)))) / np.float64(461.5)))), return_val__066_9_8__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            return_val__066_9_8__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070134222512_gen_0[i:I, j:J, k_:k_ + 1]))), (np.float64(611.21) * np.exp((((np.float64(-2372.0) * np.log((temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / np.float64(273.15)))) + ((np.float64(3147911.8) * (temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] - np.float64(273.15))) / (temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * np.float64(273.15)))) / np.float64(461.5)))), return_val__066_9_8__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            RETURN_VALUE__066_9_8__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], return_val__066_9_8__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], RETURN_VALUE__066_9_8__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            q__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (RETURN_VALUE__066_9_8__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / ((np.float64(461.5) * temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]) * density[i:I, j:J, k_:k_ + 1])), q__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            mask_140070139851968_gen_0[i:I, j:J, k_:k_ + 1] = (temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] < np.float64(273.15))
            dqdt__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139851968_gen_0[i:I, j:J, k_:k_ + 1]), ((q__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * (np.float64(-260.0) + (np.float64(2904599.0) / temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))) / (np.float64(461.5) * temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])), dqdt__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            dqdt__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070139851968_gen_0[i:I, j:J, k_:k_ + 1]))), ((q__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * (np.float64(-2372.0) + (np.float64(3147911.8) / temp__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))) / (np.float64(461.5) * temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])), dqdt__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], q__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], dqdt__8f5_62_20__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], sz_qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], sz_dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] - qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], sz_dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / (np.float64(1.0) + (tcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))), tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], sz_tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1] = (qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] > np.float64(1e-15))
            cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), (np.bitwise_not(bool(bool(False))))), (np.int64(1) == np.int64(1))), (np.float64(53800000.0) * np.exp((np.float64(0.75) * np.log((qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * density[i:I, j:J, k_:k_ + 1]))))), cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), (np.bitwise_not(bool(bool(False))))), (np.bitwise_not((np.int64(1) == np.int64(1))))), (np.int64(1) == np.int64(2))), (np.exp(((-(np.float64(2.8))) + (np.float64(0.262) * (np.float64(273.15) - temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])))) * np.float64(1000.0)), cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), (np.bitwise_not(bool(bool(False))))), (np.bitwise_not((np.int64(1) == np.int64(1))))), (np.bitwise_not((np.int64(1) == np.int64(2))))), (np.int64(1) == np.int64(3))), (np.exp(((-(np.float64(0.639))) + (np.float64(12.96) * ((qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]) - np.float64(1.0))))) * np.float64(1000.0)), cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), (np.bitwise_not(bool(bool(False))))), (np.bitwise_not((np.int64(1) == np.int64(1))))), (np.bitwise_not((np.int64(1) == np.int64(2))))), (np.bitwise_not((np.int64(1) == np.int64(3))))), (np.int64(1) == np.int64(4))), ((np.float64(0.005) * np.exp((np.float64(0.304) * (np.float64(273.15) - temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])))) * np.float64(1000.0)), cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), (np.bitwise_not(bool(bool(False))))), (np.bitwise_not((np.int64(1) == np.int64(1))))), (np.bitwise_not((np.int64(1) == np.int64(2))))), (np.bitwise_not((np.int64(1) == np.int64(3))))), (np.bitwise_not((np.int64(1) == np.int64(4))))), ((np.float64(1e-05) * np.exp((np.float64(0.5) * (np.float64(273.15) - temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])))) * np.float64(1000.0)), cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            pca__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), bool(bool(False))), np.float64(1021225517.0653788), pca__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            pcb__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), bool(bool(False))), np.float64(1151.5921289318228), pcb__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            mu__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), bool(bool(False))), np.float64(3.35), mu__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), bool(bool(False))), ((pca__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / pcb__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]) * np.exp(((mu__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / (mu__444_89_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + np.float64(np.int64(3)))) * np.log(((np.float64(np.int64(6)) * density[i:I, j:J, k_:k_ + 1]) * qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))))), cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), bool(bool(False))), (cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / density[i:I, j:J, k_:k_ + 1]), cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), (((((np.float64(150.0) * dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]) * np.float64(4.0)) * np.float64(11.9)) * np.exp((np.float64(0.5) * np.log(((qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * density[i:I, j:J, k_:k_ + 1]) * cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))))) / ((((qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * density[i:I, j:J, k_:k_ + 1]) * np.power((tcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * cvm__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), np.float64(np.int64(2)))) / ((np.float64(0.024) * np.float64(461.5)) * np.power(temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], np.float64(np.int64(2))))) + (np.float64(1.0) / np.float64(2.25e-05)))).astype(np.int64), pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_pidep0__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070139854032_gen_0[i:I, j:J, k_:k_ + 1]), pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1].astype(np.float64), sz_pidep0__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1] = (dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] > np.float64(np.int64(0)))
            tc__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), (np.float64(273.15) - temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), tc__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qi_gen__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), (np.float64(4.92e-11) * np.exp((np.float64(1.33) * np.log((np.float64(1000.0) * np.exp((np.float64(0.1) * tc__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))))))), qi_gen__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), (np.int64(3) == np.int64(1))), (qi_gen__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / density[i:I, j:J, k_:k_ + 1]), qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), (np.bitwise_not((np.int64(3) == np.int64(1))))), (np.int64(3) == np.int64(2))), ((qi_gen__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * np.minimum(np.float64(1.0), (np.float64(0.1) * tc__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))) / density[i:I, j:J, k_:k_ + 1]), qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), (np.bitwise_not((np.int64(3) == np.int64(1))))), (np.bitwise_not((np.int64(3) == np.int64(2))))), (np.int64(3) == np.int64(3))), ((np.float64(1.82e-06) * np.minimum(np.float64(1.0), (np.float64(0.1) * tc__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))) / density[i:I, j:J, k_:k_ + 1]), qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(np.bitwise_and(np.bitwise_and(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), (np.bitwise_not((np.int64(3) == np.int64(1))))), (np.bitwise_not((np.int64(3) == np.int64(2))))), (np.bitwise_not((np.int64(3) == np.int64(3))))), ((np.maximum(qi_gen__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], np.float64(1.82e-06)) * np.minimum(np.float64(1.0), (np.float64(0.1) * tc__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))) / density[i:I, j:J, k_:k_ + 1]), qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], sz_qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), np.minimum(tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], np.minimum(np.maximum((qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] - qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1].astype(np.float64)), (tc__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / tcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))), sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_sink1__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], sz_sink1__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            dep[i:I, j:J] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]), (dep[i:I, j:J] + (sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * delp[i:I, j:J, k_:k_ + 1])), dep[i:I, j:J])
            b__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]))), np.float64(184.0), b__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            diff__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]))), np.where(((temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] - b__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]) > np.float64(np.int64(0))), (temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] - b__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), np.float64(np.int64(0))), diff__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            RETURN_VALUE__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]))), diff__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], RETURN_VALUE__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]))), (pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1].astype(np.float64) * np.minimum(np.float64(np.int64(1)), (RETURN_VALUE__1fb_125_35__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * np.float64(0.2)))).astype(np.int64), pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]))), pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1].astype(np.float64), sz_pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]))), np.maximum(pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1].astype(np.float64), np.maximum(tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], (-(qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])))), sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sz_sink2__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]))), sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], sz_sink2__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            sub[i:I, j:J] = np.where(np.bitwise_and(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (np.bitwise_not(mask_140070136950448_gen_0[i:I, j:J, k_:k_ + 1]))), (sub[i:I, j:J] - (sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * delp[i:I, j:J, k_:k_ + 1])), sub[i:I, j:J])
            delta_vapor__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (-(sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])), delta_vapor__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            delta_liquid__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], np.float64(0.0), delta_liquid__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            delta_rain__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], np.float64(0.0), delta_rain__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            delta_snow__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], np.float64(0.0), delta_snow__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            delta_graupel__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], np.float64(0.0), delta_graupel__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + delta_vapor__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + delta_liquid__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + delta_rain__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + sink__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + delta_snow__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + delta_graupel__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            q_l__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), q_l__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            q_solid__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], ((qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]) + qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), q_solid__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            cvm__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (((np.float64(1.0) + (qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * np.float64(1.9294822660441782))) + (q_l__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * np.float64(5.87833600445962))) + (q_solid__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * np.float64(2.934987108912271))), cvm__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            tk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (((te[i:I, j:J, k_:k_ + 1] - (np.float64(4562.7071632638845) * qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])) + (np.float64(-339.0882865305553) * ((qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]) + qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]))) / cvm__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), tk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            lcpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], ((np.float64(4562.7071632638845) + (np.float64(-3.948853738415442) * tk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])) / cvm__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), lcpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            icpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], ((np.float64(-339.0882865305553) + (np.float64(2.943348895547349) * tk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])) / cvm__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), icpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            tcpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], ((np.float64(4223.618876733329) + ((np.float64(-3.948853738415442) + np.float64(2.943348895547349)) * tk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])) / cvm__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), tcpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            a__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], np.float64(273.15), a__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            diff__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], np.where(((a__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] - tk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]) > np.float64(np.int64(0))), (a__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] - tk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]), np.float64(np.int64(0))), diff__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            RETURN_VALUE__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], diff__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], RETURN_VALUE__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            tcp3__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], (lcpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] + (icpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] * np.minimum(np.float64(1.0), (RETURN_VALUE__1fb_49_13__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] / (np.float64(273.15) - np.float64(233.14999999999998)))))), tcp3__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            cvm__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], cvm__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], cvm__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], tk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            lcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], lcpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], lcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            icpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], icpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], icpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            tcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], tcpk__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], tcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            tcp3__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1] = np.where(mask_140070138870128_gen_0[i:I, j:J, k_:k_ + 1], tcp3__f2f_144_12__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1], tcp3__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1])
            qvapor[i:I, j:J, k_:k_ + 1] = qvapor__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            qliquid[i:I, j:J, k_:k_ + 1] = qliquid__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            qrain[i:I, j:J, k_:k_ + 1] = qrain__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            qice[i:I, j:J, k_:k_ + 1] = qice__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            qsnow[i:I, j:J, k_:k_ + 1] = qsnow__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            qgraupel[i:I, j:J, k_:k_ + 1] = qgraupel__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            cloud_ice_nuclei[i:I, j:J, k_:k_ + 1] = cloud_ice_nuclei__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            temperature[i:I, j:J, k_:k_ + 1] = temperature__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            cvm[i:I, j:J, k_:k_ + 1] = cvm__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            lcpk[i:I, j:J, k_:k_ + 1] = lcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            icpk[i:I, j:J, k_:k_ + 1] = icpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            tcpk[i:I, j:J, k_:k_ + 1] = tcpk__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            tcp3[i:I, j:J, k_:k_ + 1] = tcp3__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            dep[i:I, j:J] = dep[i:I, j:J]
            sub[i:I, j:J] = sub[i:I, j:J]
            qsi[i:I, j:J, k_:k_ + 1] = sz_qsi__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            dqdt[i:I, j:J, k_:k_ + 1] = sz_dqdt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            pidep0[i:I, j:J, k_:k_ + 1] = sz_pidep0__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            pidep[i:I, j:J, k_:k_ + 1] = sz_pidep__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            qi_crt[i:I, j:J, k_:k_ + 1] = sz_qi_crt__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            sink1[i:I, j:J, k_:k_ + 1] = sz_sink1__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            sink2[i:I, j:J, k_:k_ + 1] = sz_sink2__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            tmp[i:I, j:J, k_:k_ + 1] = sz_tmp__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            dq[i:I, j:J, k_:k_ + 1] = sz_dq__e56_244_16_gen_0[i:I, j:J, k_:k_ + 1]
            # --- end horizontal block --

        # --- end vertical block ---
    