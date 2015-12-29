#!/usr/bin/python
# -*- coding: utf-8 -*-

# Multi Image FFT Filter
# <https://github.com/simondona/img-fft-filter-bec-tn>
# Copyright (C) 2014-2015  Simone Donadello
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

import pickle



class DefaultProgSettings(object):
    """
    The settings of the program are stored into this class in order to be
    saved and recalled.
    """

    def __init__(self, fname=".program.conf"):
        """
        Instance a new settings class with defaults values for a given
        filename.

        Args:
            fname (str): name of the file where these settings are stored
        """
        #filename where to save these configurations
        self.config_fname = str(fname)

        #frames dimensions
        self.frame_width = 0
        self.frame_height = 0
        self.frame_width_old = -1
        self.frame_height_old = -1
        self.frame_number = 1
        self.frame_colums = 1
        self.frame_rows = 1

        #number of the image in the sis file to be considered (0 or 1)
        self.sis_img_number = 1

        #skip n frames from the beginning and n frames from the end
        self.skip_frames = (0, 0)

        #image filename
        self.img_fname = "default.sis"

        #frames analysis toggles
        self.img_compensate_max = True
        self.img_compensate_bkg = True

        #mask filename
        self.mask_fname = ""

        #colormaps
        self.fft_cmap = "jet"
        self.frame_cmap = "gist_stern"

        #colors limits in the reprentation
        self.vlim_img_corr = (0.0, 0.0)
        self.vlim_fft_corr = (0.0, 0.0)

        #show toggles
        self.show_mask = True
        self.show_regions = True
        self.show_contour = True
        self.show_grid = True
        self.show_residuals = False

        #gaussian filter
        self.gauss_filter = (False, 0.0)

        #frame number used for the fft display: see the FFTCanvas definition
        #(-1 is the average, -2 is no frame)
        self.fft_frame_n = -1

        #sliced for the areas of the maximum, the main part, and the background
        self.slices_max = (slice(None, None), slice(None, None))
        self.slices_main = (slice(None, None), slice(None, None))
        self.slices_bkg = (slice(None, None), slice(None, None))

        #intervals of time between the frames
        self.time_delta = tuple()

        #movie settings
        self.movie_fps = 1.0 #frame per seconds
        self.movie_frames = True #frames images included in the movie
        self.movie_integrate_x = False #integral along x included in the movie
        self.movie_integrate_y = False #integral along y included in the movie
        self.movie_fading = False #smooth transition between frames
        self.movie_subframes = 10 #number of subframes for each transition
                                  #(3 will be used for the smoothing,
                                  #the remaining are fixed)
        self.movie_dpi = 100.0 #resolution of the movie
        self.movie_writer = "libav" #rendering writer (see create_movie
                                    #definition)

        #bitmap output settings
        #for the reorder encoding see ImageLab.save_frames_bitmap
        self.bitmap_reorder = "original"
        self.bitmap_hor_crop = (0, 0)
        self.bitmap_ver_crop = (0, 0)

        #sis save settings
        self.sis_save_fname = self.img_fname
        self.sis_save_suffix = "_filtered"
        #for the format encoding see ImageLab.sis_writeimg
        self.sis_save_format = "original"


    def save_settings(self):
        """
        Saves current settings into the config file.
        """
        try:
            with open(self.config_fname, "wb") as fid:
                pickle.dump(vars(self), fid)
        except IOError:
            print "Error while writing settings into file '" + \
                                                     self.config_fname + "'"

    def load_settings(self):
        """
        Loads settings into current class from the config file.
        """
        try:
            with open(self.config_fname, "rb") as fid:
                var = pickle.load(fid)
            for key, val in var.items():
                setattr(self, key, val)
        except IOError:
            print "Warning: no config file '" + self.config_fname + "' found."
