#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2014-2016  Simone Donadello
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# pylint: disable=E1101

from scipy.optimize import curve_fit
import numpy as np



def gaussian((X, Y), mx, my, sx, sy):
    """
    Returns the result of a Gaussian.

    Args:
        (X, Y) (tuple of np.array): matrices of the coordinates to be combined,
        usually results of np.meshgrid
        mx (float): horizontal center
        my (float): vertical center
        sx (float): horizontal sigma
        sy (float): vertical sigma
    """
    return np.exp(- (X-mx)**2 / (2*sx**2) - (Y-my)**2 / (2*sy**2))


def thomasfermi((X, Y), mx, my, rx, ry):
    """
    Returns the result of a Thomas-Fermi function (inverted parabola).

    Args:
        (X, Y) (tuple of np.array): matrices of the coordinates to be combined,
        usually results of np.meshgrid
        mx (float): horizontal center
        my (float): vertical center
        rx (float): horizontal TF radius
        ry (float): vertical TF radius
    """
    b = (1 - ((X-mx)/rx)**2 - ((Y-my)/ry)**2)
    b = np.maximum(b, 0)
    b = np.sqrt(b)
    return b**3


class Fitting(object):
    """
    Base class for fitting routines. It has some common methods, other must be
    overridden for the specific fitting types with child classes.
    """

    def __init__(self, img0, par0=None):
        """
        Initialize the fitting routine with a given image.

        Args:
            img0 (np.array): image for the fit
            par0 (list): initial guess for the fit (optional, not yet used)
        """
        self.img0 = img0
        self.par0 = par0 #TODO: not used

        #calculates the matrices with the x and y coordinates
        ny, nx = self.img0.shape
        x = np.arange(nx)
        y = np.arange(ny)
        self.X, self.Y = np.meshgrid(x, y)

        #the fitted image result is initialized to None
        self.fitted = None

        #list of the parameter string names, must be implemented
        self.par_names = tuple([])

    def guess_gauss_par0(self, slc_main, slc_max, slc_bkg):
        """
        Guess and returns the initial gaussian parameters from the slices

        Args:
            slc_main (tuple): tuple of 2 slices for the coordinates of the main
            area for the gaussian center guess
            slc_max (tuple): tuple of 2 slices for the coordinates of the
            maximal area for the gaussian amplitude guess
            slc_bkg (tuple): tuple of 2 slices for the coordinates of the
            background area for the gaussian offset guess

        Returns:
            Tuple with the Gaussian guessed parameters
        """

        height, width = self.img0[slc_main].shape

        #center
        xm = np.mean(slc_main[1].indices(self.img0.shape[1])[0:-1])
        ym = np.mean(slc_main[0].indices(self.img0.shape[0])[0:-1])

        offs = np.mean(self.img0[slc_bkg])
        amp1 = np.mean(self.img0[slc_max])

        return (offs, amp1, xm, ym, width/4.0, height/4.0)

    def guess_par0(self, *args, **kwargs):
        """
        Parameters guess for the specific function (must be overridden).
        """
        pass

    def function(self, *args, **kwargs):
        """
        Specific function (must be overridden).
        """
        pass

    def fit(self):
        """
        Performs the fitting operations and returns a dictionary with the
        fitted parameters.
        """
        frame = self.img0
        #TODO handle when there is a wrong fit and set fit options
        #TODO consider parameters uncertanties
        try:
            results = curve_fit(self.function, (self.X, self.Y), frame.ravel(),
                                p0=self.par0)
        except RuntimeError:
            print "Error while fitting"
            results = [self.par0, None]

        return dict(zip(self.par_names, results[0]))


class Gauss2d(Fitting):
    """
    Gaussian 2D fit.
    """
    def __init__(self, img0, par0=None):
        super(Gauss2d, self).__init__(img0, par0)
        self.par_names = ["offs", "amp1", "mx", "my", "sx", "sy"]

    def function(self, (X, Y), offs, amp1, mx, my, sx, sy):
        """
        Implements the gaussian fitting function.
        (see gaussian() and thomasfermi())
        """
        self.fitted = amp1*gaussian((X, Y), mx, my, sx, sy) + offs
        return self.fitted.ravel()

    def guess_par0(self, slc_main, slc_max, slc_bkg):
        """
        Implements the gaussian parameter guess from slices.
        (see Fitting.guess_gauss_par0())
        """
        offs, amp1, mx, my, sx, sy = self.guess_gauss_par0(slc_main,
                                                           slc_max,
                                                           slc_bkg)
        par0 = (offs, amp1, mx, my, sx, sy)

        self.par0 = par0
        return par0


class ThomasFermi2d(Fitting):
    """
    Thomas-Fermi 2D fit (inverted parabola).
    """
    def __init__(self, img0, par0=None):
        super(ThomasFermi2d, self).__init__(img0, par0)
        self.par_names = ["offs", "amp1", "mx", "my", "rx", "ry"]

    def function(self, (X, Y), offs, amp2, mx, my, rx, ry):
        """
        Implements the Thomas-Fermi fitting function.
        (see gaussian() and thomasfermi())
        """
        self.fitted = amp2*thomasfermi((X, Y), mx, my, rx, ry) + offs
        return self.fitted.ravel()

    def guess_par0(self, slc_main, slc_max, slc_bkg):
        """
        Implements the Thomas-Fermi parameter guess.
        (see Fitting.guess_gauss_par0())
        """
        offs, amp1, mx, my, sx, sy = self.guess_gauss_par0(slc_main,
                                                           slc_max,
                                                           slc_bkg)
        par0 = (offs, amp1, mx, my, sx*2.0, sy*2.0)

        self.par0 = par0
        return par0


class Bimodal2d(Gauss2d, ThomasFermi2d):
    """
    Gaussian+Thomas Fermi bimodal 2D fit.
    """
    def __init__(self, img0, par0=None):
        super(Bimodal2d, self).__init__(img0, par0)
        self.par_names = ["offs", "amp1", "mx", "my", "sx", "sy",
                          "amp2", "rx", "ry"]

    def function(self, (X, Y), offs, amp1, mx, my, sx, sy, amp2, rx, ry):
        """
        Implements the bimodal fitting function.
        (see gaussian() and thomasfermi())
        """
        self.fitted = amp1*gaussian((X, Y), mx, my, sx, sy) +\
                      amp2*thomasfermi((X, Y), mx, my, rx, ry) + offs
        return self.fitted.ravel()

    def guess_par0(self, slc_main, slc_max, slc_bkg):
        """
        Implements the bimodal parameter guess.
        (see Fitting.guess_gauss_par0())
        """
        offs, amp1, mx, my, sx, sy = self.guess_gauss_par0(slc_main,
                                                           slc_max,
                                                           slc_bkg)
        par0 = (offs, amp1/2.0, mx, my, sx, sy, amp1/2.0, sx*2.0, sy*2.0)

        self.par0 = par0
        return par0
