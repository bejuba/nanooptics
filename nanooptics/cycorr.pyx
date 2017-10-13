cimport cython
import numpy as _np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as _np


def corr(channel, timestamp, cutofftime=1e-6, resolution=4e-12, chan0=0, chan1=1, normalize=True):
    if channel.dtype != 'uint8':
        channel = _np.uint8(channel)
    timestamp = _np.uint64(timestamp / resolution)
    cutofftime = _np.int64(cutofftime / resolution)
    t = (_np.arange(0, 2 * cutofftime) - cutofftime + 1) * resolution
    g2 = optcorr(channel, timestamp, cutofftime, chan0, chan1)
    g2_error = _np.sqrt(g2)
    if normalize:
        measurement_time = timestamp[-1]
        counts0 = _np.sum([channel == 0])
        counts1 = _np.sum([channel == 1])
        norm_factor = (
            ( measurement_time - cutofftime )
            / ( counts0 * counts1 )
        )
        g2 = norm_factor * g2
        g2_error = norm_factor * g2_error
    return t, g2, g2_error

@cython.cdivision
@cython.boundscheck(False)
@cython.wraparound(False)
cdef optcorr(_np.ndarray[_np.uint8_t, ndim=1] channel,
             _np.ndarray[_np.uint64_t, ndim=1] timestamp,
             _np.uint64_t cutofftime,
             _np.int_t chan0=0,
             _np.int_t chan1=1):
    cdef _np.uint64_t last_t2_index = len(timestamp) - 1
    cdef _np.uint_t i, j, tau
    cdef g2_unnormalized = _np.zeros(2 * cutofftime, dtype=_np.int32)

    for i in range(last_t2_index):
        if channel[i] == chan0:
            for j in range(i+1, last_t2_index):
                tau = (timestamp[j] - timestamp[i])
                if tau >= cutofftime:
                    break
                if channel[j] == chan1:
                    g2_unnormalized[tau+cutofftime] += 1
        elif channel[i] == chan1:
            for j in range(i+1, last_t2_index):
                tau = (timestamp[j] - timestamp[i])
                if tau > cutofftime:
                    break
                if channel[j] == chan0:
                    g2_unnormalized[cutofftime-tau] += 1
    return g2_unnormalized

def corr3(channel,
          timestamp,
          resolution=4e-12,
          chan0=0,
          chan1=1, chan1_min=100e-9, chan1_max=200e-9,
          chan2=2, chan2_min=100e-9, chan2_max=200e-9):
    if channel.dtype != 'uint8':
        channel = _np.uint8(channel)
    timestamp = _np.uint64(timestamp / resolution)
    chan1_min = _np.int64(chan1_min / resolution)
    chan1_max = _np.int64(chan1_max / resolution)
    chan2_min = _np.int64(chan2_min / resolution)
    chan2_max = _np.int64(chan2_max / resolution)
    g2 = optcorr3(channel,
                  timestamp,
                  chan0,
                  chan1,
                  chan1_min,
                  chan1_max,
                  chan2,
                  chan2_min,
                  chan2_max)
    g2_error = _np.sqrt(g2)
    return g2, g2_error


@cython.cdivision
@cython.boundscheck(False)
@cython.wraparound(False)
cdef optcorr3(_np.ndarray[_np.uint8_t, ndim=1] channel,
             _np.ndarray[_np.uint64_t, ndim=1] timestamp,
             _np.int_t chan0,
             _np.int_t chan1,
             _np.int_t chan1_min,
             _np.int_t chan1_max,
             _np.int_t chan2,
             _np.int_t chan2_min,
             _np.int_t chan2_max):
    cdef _np.uint64_t last_t2_index = len(timestamp) - 1
    cdef _np.uint_t i, j, k, tau1, tau2
    cdef g2 = _np.zeros((chan1_max - chan1_min, chan2_max - chan2_min), dtype=_np.int32)

    for i in range(last_t2_index):
        if channel[i] == chan0:
            for j in range(i+1, last_t2_index):
                tau1 = timestamp[j] - timestamp[i]
                if ( (tau1 < chan1_min) or (tau1 >= chan1_max) ):
                    break
                elif channel[j] == chan1:
                    for k in range(j+1, last_t2_index):
                        tau2 = timestamp[k] - timestamp[i]
                        if ( (tau2 < chan2_min) or (tau2 >= chan2_max) ):
                            break
                        elif channel[k] == chan2:
                            g2[tau1-chan1_min][tau2-chan2_min] += 1
    return g2


## his has to be redone for second resolution
# cpdef timetrace(_np.ndarray[_np.uint_t, ndim=2] arrivaltimes, _np.uint_t binwidth):
#     cdef DTYPE_t lengtharrivaltimes = len(arrivaltimes) - 1
#     trace_double = _np.zeros((int(arrivaltimes[-1, 1] / binwidth), 2))
#     cdef _np.ndarray[_np.uint_t, ndim=2] trace = trace_double.astype(_np.uint)
#     cdef _np.uint_t i
#
#     for i in range(lengtharrivaltimes):
#         trace[arrivaltimes[i, 1] / binwidth, arrivaltimes[i, 0]] += 1
#
#     t = _np.arange(0, arrivaltimes[-1, 1] / binwidth, binwidth)
#     return t, trace
