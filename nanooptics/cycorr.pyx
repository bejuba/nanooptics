import numpy as np
cimport numpy as np
from cython.parallel import prange
cimport cython


def corr(channel, timestamp, cutofftime=1e-6, resolution=4e-12, chan0=0, chan1=1, normalize=True, parallelize = True):
    if channel.dtype != 'uint8':
        channel = np.uint8(channel)
    timestamp = np.uint64(timestamp / resolution)
    cutofftime = np.uint64(cutofftime / resolution)
    t = (np.arange(0, 2 * cutofftime) - cutofftime + 1) * resolution
    if parallelize:
        g2 = poptcorr(channel, timestamp, cutofftime, chan0, chan1)
    else:
        g2 = optcorr(channel, timestamp, cutofftime, chan0, chan1)
    g2_error = np.sqrt(g2)
    if normalize:
        measurement_time = timestamp[-1] - timestamp[0]
        counts0 = np.sum([channel == chan0])
        counts1 = np.sum([channel == chan1])
        # by default the result of _np.sum is an int32 this can lead to overflows when multiplying large count numbers
        # so we cast them to double to prevent this problem.
        norm_factor = (
            ( measurement_time - cutofftime )
            / ( np.double(counts0) * np.double(counts1) )
        )
        g2 = norm_factor * g2
        g2_error = norm_factor * g2_error
    return np.rec.array([t, g2, g2_error], dtype=[('tau', float), ('g2', float), ('yerr', float)])

@cython.cdivision
@cython.boundscheck(False)
@cython.wraparound(False)
cdef optcorr(np.ndarray[np.uint8_t, ndim=1] channel,
             np.ndarray[np.uint64_t, ndim=1] timestamp,
             np.uint64_t cutofftime,
             np.int_t chan0=0,
             np.int_t chan1=1):
    cdef int timestamp_len = len(timestamp)
    cdef int i, j
    cdef np.uint64_t tau
    cdef np.double_t[:] g2_unnormalized = np.zeros(2 * cutofftime, dtype=np.double)

    for i in range(timestamp_len):
        if channel[i] == chan0:
            for j in range(i+1, timestamp_len):
                tau = timestamp[j] - timestamp[i]
                if tau >= cutofftime:
                    break
                if channel[j] == chan1:
                    g2_unnormalized[tau+cutofftime] += 1
        elif channel[i] == chan1:
            for j in range(i+1, timestamp_len):
                tau = timestamp[j] - timestamp[i]
                if tau > cutofftime:
                    break
                if channel[j] == chan0:
                    g2_unnormalized[cutofftime-tau] += 1
    return g2_unnormalized

@cython.cdivision
@cython.boundscheck(False)
@cython.wraparound(False)
cdef poptcorr(np.ndarray[np.uint8_t, ndim=1] channel,
             np.ndarray[np.uint64_t, ndim=1] timestamp,
             np.uint64_t cutofftime,
             np.int_t chan0=0,
             np.int_t chan1=1):
    cdef int timestamp_len = len(timestamp)
    cdef int i, j
    cdef np.uint64_t tau
    cdef np.double_t[:] g2_unnormalized = np.zeros(2 * cutofftime, dtype=np.double)

    for i in prange(timestamp_len, nogil=True):
        if channel[i] == chan0:
            for j in range(i+1, timestamp_len):
                tau = timestamp[j] - timestamp[i]
                if tau >= cutofftime:
                    break
                if channel[j] == chan1:
                    g2_unnormalized[tau+cutofftime] += 1
        elif channel[i] == chan1:
            for j in range(i+1, timestamp_len):
                tau = timestamp[j] - timestamp[i]
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
        channel = np.uint8(channel)
    timestamp = np.uint64(timestamp / resolution)
    chan1_min = np.uint64(chan1_min / resolution)
    chan1_max = np.uint64(chan1_max / resolution)
    chan2_min = np.uint64(chan2_min / resolution)
    chan2_max = np.uint64(chan2_max / resolution)
    g2 = optcorr3(channel,
                  timestamp,
                  chan0,
                  chan1,
                  chan1_min,
                  chan1_max,
                  chan2,
                  chan2_min,
                  chan2_max)
    g2_error = np.sqrt(g2)
    return [g2, g2_error]


@cython.cdivision
@cython.boundscheck(False)
@cython.wraparound(False)
cdef optcorr3(np.ndarray[np.int_t, ndim=1] channel,
             np.ndarray[np.uint64_t, ndim=1] timestamp,
             np.int_t chan0,
             np.int_t chan1,
             np.uint64_t chan1_min,
             np.uint64_t chan1_max,
             np.int_t chan2,
             np.uint64_t chan2_min,
             np.uint64_t chan2_max):
    cdef int timestamp_len = len(timestamp)
    cdef int i, j, k
    cdef np.uint64_t tau1, tau2
    cdef np.ndarray[np.uint64_t, ndim=2] g2 = np.zeros((chan1_max - chan1_min, chan2_max - chan2_min), dtype=np.uint64)

    for i in range(timestamp_len):
        if channel[i] == chan0:
            for j in range(i+1, timestamp_len):
                tau1 = timestamp[j] - timestamp[i]
                if (tau1 < chan1_min) or (tau1 >= chan1_max):
                    break
                elif channel[j] == chan1:
                    for k in range(j+1, timestamp_len):
                        tau2 = timestamp[k] - timestamp[i]
                        if (tau2 < chan2_min) or (tau2 >= chan2_max):
                            break
                        elif channel[k] == chan2:
                            g2[tau1-chan1_min][tau2-chan2_min] += 1
    return g2


def syncdiff(channel, timestamp, syncchan, reverse=False):
    if channel.dtype != np.int32:
        channel = np.int32(channel)
        print('internally converted channel to int32 for cython')
    if timestamp.dtype != np.uint64:
        timestamp = np.uint64(timestamp)
        print('internally converted timestamp to uint64 for cython')

    if reverse:
        channel = channel[::-1]
        timestamp = timestamp[::-1]

    #syncchan = np.uint8(syncchan)
    #print('internally converted syncchan to uint8 for cython')

    if reverse:
        return np.array(cysyncdiff(channel, timestamp, syncchan, reverse))[::-1]
    else:
        return np.array(cysyncdiff(channel, timestamp, syncchan, reverse))

    
@cython.cdivision
@cython.boundscheck(False)
@cython.wraparound(False)
cdef cysyncdiff(np.ndarray[np.int32_t, ndim=1] channel,
                np.ndarray[np.uint64_t, ndim=1] timestamp,
                np.int32_t syncchan, reverse):
    cdef np.uint64_t last_t2_index = len(timestamp)
    cdef np.uint64_t i
    cdef np.uint64_t lastsync = timestamp[0]
    cdef np.uint64_t[:] timediff = timestamp.copy()

    if reverse:
        for i in range(last_t2_index):
            if channel[i] == syncchan:
                lastsync = timestamp[i]
            else:
                timediff[i] = lastsync-timediff[i]
    else:
        for i in range(last_t2_index):
            if channel[i] == syncchan:
                lastsync = timestamp[i]
            else:
                timediff[i] -= lastsync
    return timediff





