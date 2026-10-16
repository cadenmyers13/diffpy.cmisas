#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Characteristic function calculator backed by sasmodels."""

__all__ = ["SASCF"]

import numpy
from numpy import ceil, log2, pi
from numpy.fft import fftfreq, ifft

from diffpy.srfit.fitbase.calculator import Calculator


class SASCF(Calculator):
    """Calculator class for characteristic functions from sas-models.

    This class wraps a sas.models.BaseModel to calculate I(Q) related to
    nanoparticle shape. This I(Q) is inverted to f(r) according to:
    f(r) = 1 / (4 pi r) * SINFT(I(Q)),
    where "SINFT" represents the sine Fourier transform.

    Attributes
    ----------
    _model
        BaseModel object this adapts.
    Managed Parameters
        These depend on the parameters of the BaseModel object held by _model.
        They are created from the 'params' attribute of the BaseModel. If a
        dispersion is set for the BaseModel, the dispersion "width" will be
        accessible under "<parname>_width", where <parname> is the name a
        parameter adjusted by dispersion.
    """

    def __init__(self, name, model):
        """Initialize the generator.

        Parameters
        ----------
        name : str
            The name for the SASCF.
        model : BaseModel
            The sas.models.BaseModel object this adapts.
        """
        Calculator.__init__(self, name)

        self._model = model

        from diffpy.cmisas.sasparameter import SASParameter

        # Wrap normal parameters
        for parname in model.params:
            par = SASParameter(parname, model)
            self.add_parameter(par)

        # Wrap dispersion parameters
        for parname in model.dispersion:
            name = parname + "_width"
            parname += ".width"
            par = SASParameter(name, model, parname)
            self.add_parameter(par)

        return

    def __call__(self, r):
        """Calculate the characteristic function from the transform of
        the BaseModel."""
        # Determine q-values.
        # We want very fine r-spacing so we can properly normalize f(r). This
        # equates to having a large qmax so that the Fourier transform is
        # finely spaced. We also want the calculation to be fast, so we pick
        # qmax such that the number of q-points is a power of 2. This allows us
        # to use the fft.
        #
        # The initial dr is somewhat arbitrary, but using dr = 0.01 allows for
        # the f(r) calculated from a particle of diameter 50, over r =
        # arange(1, 60, 0.1) to agree with the spherical_particle with Rw <
        # 1e-4%.
        #
        # We also have to make a q-spacing small enough to compute out to at
        # least the size of the signal.
        raise NotImplementedError(
            "As of release 3.2.0, SAS characteristic functions are not working"
            + " but we hope to have them working again in a future release."
        )

        dr = min(0.01, r[1] - r[0])
        ed = 2 * self._model.calculate_ER()

        # Check for nans. If we find any, then return zeros.
        if numpy.isnan(ed).any():
            y = numpy.zeros_like(r)
            return y

        rmax = max(ed, 2 * r[-1])
        dq = pi / rmax
        qmax = pi / dr
        numpoints = int(2 ** (ceil(log2(qmax / dq))))
        qmax = dq * numpoints

        # Calculate F(q) = q * I(q) from model
        q = fftfreq(int(qmax / dq)) * qmax
        fq = q * self._model.evalDistribution(q)

        # Calculate g(r) and the effective r-points
        rp = fftfreq(numpoints) * 2 * pi / dq
        # Note sine transform = imaginary part of ifft
        gr = ifft(fq).imag

        # Calculate full-fr for normalization
        assert rp[0] == 0.0
        frp = numpy.zeros_like(gr)
        frp[1:] = gr[1:] / rp[1:]

        # Inerpolate onto requested grid, do not use data after jump in rp
        assert numpoints % 2 == 0
        nhalf = numpoints / 2
        fr = numpy.interp(r, rp[:nhalf], gr[:nhalf])
        vmask = r != 0
        fr[vmask] /= r[vmask]

        # Normalize. We approximate fr[0] by using the fact that f(r) is linear
        # at low r. By definition, fr[0] should equal 1.
        fr0 = 2 * frp[2] - frp[1]
        fr /= fr0

        # Fix potential divide-by-zero issue, fr is 1 at r == 0
        fr[~vmask] = 1

        return fr


# End of file
