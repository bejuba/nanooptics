# Data analysis library for auto correlation experiments

import lmfit as _lmfit
import numpy as _np

# import helper modules
from . import read_tcspc as _read_tcspc
from . import cycorr


def g2(tau, t_antibunch=10, a_bunch=None, t_bunch=None):
    """"
    g2 function of an n-level system
    for n > 2, n-1 bunching times and amplitudes have to be provided
    
    
    Parameters
    ----------
    tau: float array like
        g2_time delays in seconds
    t_antibunch : float
        antibunching g2_time constant
    a_bunch : float
        n-2 bunching amplitudes
    t_bunch : float
        n-2 bunching g2_time constants
    
    Returns
    -------
    float
        value of the g2 function at tau
    """
    if a_bunch is None:
        a_bunch = []
    if t_bunch is None:
        t_bunch = []

    antibunching = (1 + _np.sum(a_bunch)) * _np.exp(- tau / t_antibunch)
    bunching = sum([a * _np.exp(-tau / t) for a, t in zip(a_bunch, t_bunch)])
    return 1 - antibunching + bunching


def g2_bg(tau, i, b, t_antibunch, a_bunch, t_bunch):
    return (i ** 2 * g2(tau, t_antibunch, a_bunch, t_bunch) + b ** 2 + 2 * i * b) / (i + b) ** 2


def g2_bgfactor(tau, tau0, bgfactor, t_antibunch, a_bunch=None, t_bunch=None):
    tau = abs(tau - tau0)
    return (g2(tau, t_antibunch, a_bunch, t_bunch) + bgfactor ** 2 + 2 * bgfactor) / (1 + bgfactor) ** 2


def gauss(x, sigma):
    return _np.exp(-x ** 2 / sigma ** 2)


def g2_bgfactor_gauss(tau, tau0, bgfactor, sigma, t_antibunch, a_bunch=None, t_bunch=None):
    # we do not use the convoluted data at the boarders of the data set
    left = tau < (tau[0] + 3 * sigma)
    right = tau > (tau[-1] - 3 * sigma)
    orig_left = g2_bgfactor(tau[left], tau0, bgfactor, t_antibunch, a_bunch, t_bunch)
    convolved = _np.convolve(g2_bgfactor(tau,
                                         tau0, bgfactor, t_antibunch, a_bunch, t_bunch),
                             gauss(tau,
                                   sigma)
                             / sum(gauss(tau, sigma)),
                             'same')
    orig_right = g2_bgfactor(tau[right], tau0, bgfactor, t_antibunch, a_bunch, t_bunch)
    convolved[left] = orig_left
    convolved[right] = orig_right
    return convolved


def g2_specdiff(x, t0, bgfactor, a, t1, t2, c, t_d):
    return g2_bgfactor(x, t0, bgfactor, a, t1, t2) * (1 - _np.exp(-_np.abs(x - t0) / t_d) * c ** 2 / 2)


def g2_specdiff_gauss(x, t0, bgfactor, a, t1, t2, c, t_d, sigma):
    return _np.convolve(g2_specdiff(x, t0, bgfactor, a, t1, t2, c, t_d),
                        gauss(x, sigma) / sum(gauss(x, sigma)),
                        'same')


def fit_g2(tau, g2, g2_error, bgfactor=0.1, tau0=0, t_antibunch=1e-9, a_bunch=[], t_bunch=[], sigma=None, tau0_fixed=False):
    # create fit parameters
    params = _lmfit.Parameters()
    params.add('tau0', value=tau0)
    if tau0_fixed:
        params['tau0'].set(vary=False)
    params.add('bgfactor', value=bgfactor, min=0)
    params.add('t_antibunch', value=t_antibunch)
    a_params = [
        _lmfit.Parameter('A_' + str(i), value=a, min=0)
        for i, a in enumerate(a_bunch)
    ]
    t_params = [
        _lmfit.Parameter('t_' + str(i), value=t)
        for i, t in enumerate(t_bunch)
    ]
    for p in (a_params + t_params):
        params.add(p)
    params.add('g2_0', expr='( bgfactor**2 + 2*bgfactor ) / ( 1 + bgfactor )**2')

    # wrap g2_normalized in order to make it compatible to lmfit.Model
    def g2_wrapped(tau, **kwargs):
        tau0 = kwargs['tau0']
        bgfactor = kwargs['bgfactor']
        t_antibunch = kwargs['t_antibunch']
        a_bunch = [kwargs[a.name] for a in a_params]
        t_bunch = [kwargs[t.name] for t in t_params]
        if sigma is None:
            return g2_bgfactor(tau, tau0, bgfactor, t_antibunch, a_bunch, t_bunch)
        else:
            return g2_bgfactor_gauss(tau, tau0, bgfactor, sigma, t_antibunch, a_bunch, t_bunch)

    model = _lmfit.Model(g2_wrapped, independent_vars=['tau'])
    result = model.fit(tau=tau,
                       data=g2,
                       params=params,
                       weights=1/g2_error,
                       )
    return result


def timetrace(timestamp, integration_time=100e-3, counts_per_second=True):
    bins = _np.floor_divide(timestamp[-1],integration_time)
    counts, t = _np.histogram(timestamp[timestamp < integration_time*bins], bins=bins)
    t = t[:-1]
    if counts_per_second:
        counts *= integration_time
    return t, counts


class TimeTagMeasurement:
    def __init__(self, pt2_filepath, g2_cutofftime=1e-6, g2_resolution=4e-12, chan0=0, chan1=1):
        self.header, self.channels, self.timestamps = _read_tcspc.read_pt2(pt2_filepath)
        self.g2_cutofftime = g2_cutofftime
        self.g2_resolution = g2_resolution
        self.g2_measurementTime = self.timestamps[-1]
        self.g2_time, self.g2_unnormalized, g2_error = cycorr.corr(self.channels,
                                                                   self.timestamps,
                                                                   g2_cutofftime,
                                                                   g2_resolution,
                                                                   chan0,
                                                                   chan1,
                                                                   normalize=False)
        self.counts0 = _np.sum([self.channels == chan0])
        self.counts1 = _np.sum([self.channels == chan1])
        self.g2_resolution_raw = self.g2_resolution
        self.g2_time_raw = self.g2_time
        self.g2_cutofftime_raw = self.g2_cutofftime
        self.g2_unnormalized_raw = self.g2_unnormalized
        self.g2_normalized = self.normalize()
        self.g2_unnormalized_error, self.g2_normalized_error = self.calc_errors()
        self.fitresult = None

    def reset(self):
        self.g2_resolution = self.g2_resolution_raw
        self.g2_time = self.g2_time_raw
        self.g2_cutofftime = self.g2_cutofftime_raw
        self.g2_unnormalized = self.g2_unnormalized_raw
        self.g2_normalized = self.normalize()
        self.g2_unnormalized_error, self.g2_normalized_error = self.calc_errors()

    def normalize(self):
        g2_normalized = (self.g2_unnormalized
                         / _np.hstack((self.g2_time[1] - self.g2_time[0],
                                       self.g2_time[1:] - self.g2_time[:-1])
                                      )
                         * (self.g2_measurementTime - self.g2_cutofftime)
                         / (self.counts0 * self.counts1)
                         )
        return g2_normalized

    def calc_errors(self):
        g2_unnormalized_error = _np.sqrt(self.g2_unnormalized)
        g2_normalized_error = (g2_unnormalized_error
                               / _np.hstack((self.g2_time[1] - self.g2_time[0],
                                             self.g2_time[1:] - self.g2_time[:-1])
                                            )
                               * (self.g2_measurementTime - self.g2_cutofftime)
                               / (self.counts0 * self.counts1)
                               )
        return g2_unnormalized_error, g2_normalized_error

    def bin_linear(self, binfactor=100):
        self.g2_resolution = binfactor * self.g2_resolution
        bins = int(_np.floor(self.g2_cutofftime / self.g2_resolution))
        g2_unnormalized_shaped = _np.reshape(self.g2_unnormalized_raw[0:(binfactor * bins)],
                                             (bins, binfactor))
        self.g2_unnormalized = g2_unnormalized_shaped.sum(axis=1)
        self.g2_time = _np.arange(self.g2_unnormalized.size) * self.g2_resolution
        self.g2_cutofftime = self.g2_unnormalized.size * self.g2_resolution
        self.normalize()

    def bin_log(self, bins=250, offset=0):
        self.g2_unnormalized = self.g2_unnormalized[int(offset / self.g2_resolution):]

        bin_time = _np.logspace(_np.log10(self.g2_resolution),
                                _np.log10(self.g2_cutofftime),
                                bins + 1
                                )
        idx = _np.array(bin_time / self.g2_resolution, dtype=_np.int)
        idx = _np.unique(idx)
        bins = idx.size - 1
        g2 = _np.zeros(bins)
        for i in range(bins):
            g2[i] = _np.sum(self.g2_unnormalized[idx[i]: idx[i + 1]])
        self.g2_unnormalized = g2
        self.g2_time = idx[:-1] * self.g2_resolution
        self.normalize()

    def fit(self, bgfactor=0.1, tau0=0, t_antibunch=1e-9, a_bunch=[], t_bunch=[], sigma=None, tau0_fixed=False):
        self.fitresult = fit_g2(self.g2_time,
                                self.g2_normalized,
                                self.g2_normalized_error,
                                bgfactor,
                                tau0,
                                t_antibunch,
                                a_bunch,
                                t_bunch,
                                sigma,
                                tau0_fixed
                                )
        return self.fitresult

    def exclude_flashback(self, time_min, time_max):
        mask = (self.g2_time < time_min) | (self.g2_time > time_max)
        self.g2_time = self.g2_time[mask]
        self.g2_normalized = self.g2_normalized[mask]
        self.g2_unnormalized = self.g2_unnormalized[mask]
        self.g2_normalized_error = self.g2_normalized_error[mask]
        self.g2_unnormalized_error = self.g2_unnormalized_error[mask]


class StartStopMeasurement:
    # this still needs work!!
    def __init__(self, path):
        self.coincidences = _np.loadtxt(path, skiprows=10)
        self.resolution = _np.genfromtxt(path, skip_header=7, max_rows=1)
        self.time = _np.linspace(
            0,
            self.coincidences.size * self.resolution,
            num=self.coincidences.size
        )

    def bin_linear(self, binfactor=100):
        self.resolution = binfactor * self.resolution
        bins = int(_np.floor(self.time[-1] / self.resolution))
        coincidences_shaped = _np.reshape(self.coincidences[0:(binfactor * bins)],
                                          (bins, binfactor))
        self.coincidences = coincidences_shaped.sum(axis=1)
        self.time = _np.arange(self.coincidences.size) * self.resolution

    def fit(self, fitrange=None, bgfactor=0.1, tau0=0, t_antibunch=10e-12, a_bunch=[], t_bunch=[],
            sigma=None, tau0_fixed=False, plot=False, verbose=False):

        if fitrange is None:
            rangemin = 0
            rangemax = self.time.max()
        else:
            rangemin = fitrange[0]
            rangemax = fitrange[1]

        fitrange = (self.time >= rangemin) & (self.time <= rangemax)
        x = self.time[fitrange]
        y = self.coincidences[fitrange]

        # create fit parameters
        params_g2 = _lmfit.Parameters()
        params_g2.add('tau0', value=tau0)
        params_g2.add('bgfactor', value=bgfactor, min=0)
        params_g2.add('t_antibunch', value=t_antibunch, min=0)
        self.a_params = [
            _lmfit.Parameter('a_' + str(i), value=a, min=0)
            for i, a in enumerate(a_bunch)
        ]
        self.t_params = [
            _lmfit.Parameter('t_' + str(i), value=t, min=0)
            for i, t in enumerate(t_bunch)
        ]
        if sigma is not None:
            params_g2.add('sigma', value=sigma, min=0)
        if tau0_fixed:
            params_g2['tau0'].set(vary=False)

        for p in (self.a_params + self.t_params):
            params_g2.add(p)

        def g2_wrapped(x, **kwargs):
            tau0 = kwargs['tau0']
            bgfactor = kwargs['bgfactor']
            t_antibunch = kwargs['t_antibunch']
            a_bunch = [kwargs[a.name] for a in self.a_params]
            t_bunch = [kwargs[t.name] for t in self.t_params]
            if sigma is None:
                return g2_bgfactor(x, tau0, bgfactor, t_antibunch, a_bunch, t_bunch)
            else:
                return g2_bgfactor_gauss(x, tau0, bgfactor, sigma, t_antibunch, a_bunch, t_bunch)

        self.model = _lmfit.models.ExponentialModel() * _lmfit.Model(g2_wrapped)
        params_exp = self.model.left.guess(y, x=x)
        self.result_exp = self.model.left.fit(y, x=x, params=params_exp)

        g2_guess_range_max = 2 * _np.argmax(self.result_exp.residual)

        self.result_g2 = self.model.right.fit(
            _np.divide(y[0:g2_guess_range_max],
                       self.model.left.eval(x=x[0:g2_guess_range_max],
                                           params=self.result_exp.params)
                       ),
            x=x[0:g2_guess_range_max],
            params=params_g2
        )

        self.result = self.model.fit(
            y,
            x=x,
            params=(self.result_exp.params
                    + self.result_g2.params),
            weights=1 / _np.sqrt(y),
            fit_kws={'nan_policy': 'omit'}
        )

        if plot:
            self.result_exp.plot()
            self.result_g2.plot()
            self.result.plot()
        if verbose:
            print(self.result_exp.fit_report())
            print(self.result_g2.fit_report())
            print(self.result.fit_report())
        return self.result
