from cython.parallel import prange
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
    cutofftime = _np.uint64(cutofftime / resolution)
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

def pcorr(channel, timestamp, cutofftime=1e-6, resolution=4e-12, chan0=0, chan1=1, normalize=True):
    if channel.dtype != 'uint8':
        channel = _np.uint8(channel)
    timestamp = _np.uint64(timestamp / resolution)
    cutofftime = _np.uint64(cutofftime / resolution)
    t = (_np.arange(0, 2 * cutofftime) - cutofftime + 1) * resolution
    g2 = poptcorr(channel, timestamp, cutofftime, chan0, chan1)
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
    cdef int i, j
    cdef _np.uint64_t tau
    cdef _np.uint64_t[:] g2_unnormalized = _np.zeros(2 * cutofftime, dtype=_np.int64)

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
                tau = timestamp[j] - timestamp[i]
                if tau > cutofftime:
                    break
                if channel[j] == chan0:
                    g2_unnormalized[cutofftime-tau] += 1
    return g2_unnormalized

@cython.cdivision
@cython.boundscheck(False)
@cython.wraparound(False)
cdef poptcorr(_np.ndarray[_np.uint8_t, ndim=1] channel,
             _np.ndarray[_np.uint64_t, ndim=1] timestamp,
             _np.uint64_t cutofftime,
             _np.int_t chan0=0,
             _np.int_t chan1=1):
    cdef int last_t2_index = len(timestamp) - 1
    cdef int i, j
    cdef _np.uint64_t tau
    cdef _np.uint64_t [:] g2_unnormalized = _np.zeros(2 * cutofftime, dtype=_np.int32)

    for i in prange(last_t2_index, nogil=True):
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
    cdef g2 = _np.zeros((chan1_max - chan1_min, chan2_max - chan2_min), dtype=_np.int)

    for i in range(last_t2_index):
        if channel[i] == chan0:
            for j in range(i+1, last_t2_index):
                tau1 = timestamp[j] - timestamp[i]
                if (tau1 < chan1_min) or (tau1 >= chan1_max):
                    break
                elif channel[j] == chan1:
                    for k in range(j+1, last_t2_index):
                        tau2 = timestamp[k] - timestamp[i]
                        if (tau2 < chan2_min) or (tau2 >= chan2_max):
                            break
                        elif channel[k] == chan2:
                            g2[tau1-chan1_min][tau2-chan2_min] += 1
    return g2


#
# def timefilter(channel, timestamp, cutoff, resolution, syncchan,
#                chan1,
#                chan1_min,
#                chan1_max,
#                chan2,
#                chan2_min,
#                chan2_max):
#     if channel.dtype != 'uint8':
#         channel = _np.uint8(channel)
#     cutoff = _np.uint64(cutoff / resolution)
#     timestamp = _np.uint64(timestamp / resolution)
#     chan, tt = cytimefilter(channel,
#                           timestamp,
#                             cutoff,
#                           syncchan,
#                           chan1,
#                           chan1_min,
#                           chan1_max,
#                           chan2,
#                           chan2_min,
#                           chan2_max)
#     return chan, tt * resolution


# @cython.cdivision
# @cython.boundscheck(False)
# @cython.wraparound(False)
# cdef cytimefilter(_np.ndarray[_np.uint8_t, ndim=1] channel,
#              _np.ndarray[_np.uint64_t, ndim=1] timestamp,
#              _np.int_t cutoff,
#              _np.int_t sync,
#              _np.int_t chan1,
#              _np.int_t chan1_min,
#              _np.int_t chan1_max,
#              _np.int_t chan2,
#              _np.int_t chan2_min,
#              _np.int_t chan2_max):
#     cdef _np.uint64_t last_t2_index = len(timestamp) - 1
#     cdef _np.uint_t i, j, tau
#
#     for i in range(last_t2_index):
#         if channel[i] == sync:
#             for j in range(i+1, last_t2_index):
#                 if channel[j] == chan1:
#                     tau = timestamp[j] - timestamp[i]
#                     if tau > cutoff:
#                          break
#                     if (tau > chan1_min) or (tau < chan1_max):
#                         timestamp[j] = 0
#                 elif channel[j] == chan2:
#                     tau = timestamp[j] - timestamp[i]
#                     if tau > cutoff:
#                          break
#                     if (tau > chan2_min) or (tau < chan2_max):
#                         timestamp[j] = 0
#
#         elif channel[i] == chan1:
#             for j in range(i+1, last_t2_index):
#                 if channel[j] == sync:
#                     tau = timestamp[i] - timestamp[j]
#                     if tau < cutoff:
#                          break
#                     if (tau > chan1_min) or (tau < chan1_max):
#                         timestamp[i] = 0
#
#         elif channel[i] == chan2:
#             for j in range(i+1, last_t2_index):
#                 if channel[j] == sync:
#                     tau = timestamp[i] - timestamp[j]
#                     if tau < cutoff:
#                          break
#                     if (tau > chan2_min) or (tau < chan2_max):
#                         timestamp[i] = 0
#
#
#     mask = timestamp != 0
#     return channel[mask], timestamp[mask]


def syncdiff(channel, timestamp, syncchan, reverse=False):
    if channel.dtype != _np.uint8:
        channel = _np.uint8(channel)
        print('internally converted channel to uint8 for cython')
    if timestamp.dtype != _np.double:
        timestamp = _np.double(timestamp)
        print('internally converted timestamp to double for cython')

    if reverse:
        channel = channel[::-1]
        timestamp = timestamp[::-1]

    syncchan = _np.uint8(syncchan)
    print('internally converted syncchan to uint8 for cython')

    if reverse:
        return _np.array(cysyncdiff(channel, timestamp, syncchan))[::-1]
    else:
        return _np.array(cysyncdiff(channel, timestamp, syncchan))


@cython.cdivision
@cython.boundscheck(False)
@cython.wraparound(False)
cdef cysyncdiff(_np.ndarray[_np.uint8_t, ndim=1] channel,
                _np.ndarray[_np.double_t, ndim=1] timestamp,
                _np.uint8_t syncchan):
    cdef _np.uint64_t last_t2_index = len(timestamp)
    cdef _np.int_t i
    cdef _np.double_t lastsync = timestamp[0]
    cdef _np.double_t[:] timediff = timestamp.copy()

    for i in range(last_t2_index):
        if channel[i] == syncchan:
            lastsync = timestamp[i]
        else:
            timediff[i] -= lastsync

    return timediff