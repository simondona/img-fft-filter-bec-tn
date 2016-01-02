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

import numpy as np
import os

#use Agg for fast matplotlib rendering:
#must be called for first or it wont work
import matplotlib as mpl
mpl.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy import ndimage

from PyQt4 import QtCore

from libraries.results import ResultList



class ImageLab(QtCore.QObject, object):
    """
    This is the main class that stores all the numerical data, parameters and
    informations of a image of the program, and allows to elaborate it with
    fft filtering, normalization, curve fitting and other operations.
    """

    #signal emitted when a fit is completed returng the float value of the
    #progress fraction (from 0 to 1) over the total number of fit to do
    fit_done = QtCore.pyqtSignal(float)

    def __init__(self, filename,
                 frame_h=0,
                 frame_w=0,
                 frame_h_old=-1,
                 frame_w_old=-1,
                 frame_rows=1,
                 frame_cols=1,
                 frame_num=1,
                 bkg_slice=(slice(None, None), slice(None, None)),
                 max_slice=(slice(None, None), slice(None, None)),
                 main_slice=(slice(None, None), slice(None, None)),
                 cmap_img="jet",
                 cmap_fft="gist_stern",
                 sis_img_n=0,
                 vlim_img_corr=(0.0, 0.0),
                 vlim_fft_corr=(0.0, 0.0),
                 skip_frames=(0, 0),
                 fft_frame_num=-1,
                 compensate_max=True,
                 compensate_bkg=True,
                 show_residuals=False,
                 gauss_filter=(False, 0.0)):
        """
        Initialize the image class. In order to initialize the image data
        once everyting is loaded correctly, the function init() must be called
        to perform the computations.

        Args:
            filename (string): name of the sis file to open
            frame_h (int): single frame height
            frame_w (int): single frame width
            frame_h_old (int): single frame height of the previous image
            frame_w_old (int): single frame width of the previous image
            frame_rows (int): number of frame rows
            frame_cols (int): number of frame columns
            frame_num (int): total number of frames
            bkg_slice (tuple): tuple of two slice objects (one for vertical one
            for horizontal) that defines the image area of the background
            max_slice (tuple): slices of the maximum image area
            main_slice (tuple): slices of the area of the image intresting part
            cmap_img (string): name of a matplotlib colormap for the frames panel
            cmap_fft (string): name of a matplotlib colormap for the fft panel
            sis_img_n (int): 0 or 1 part number to consider in the sis file
            (each file has 2 parts)
            vlim_img_corr (tuple): tuple of floats with the lower/upper
            colorscale corrections for frames (in the scale of an original
            scale extension)
            vlim_fft_corr (tuple): tuple of floats with colorscale fft corrections
            skip_frames (tuple): tuple with the number of frames to skip in the
            grid at the beginning and at the end
            fft_frame_num (int): frame number to be used in for the ffr panel
            compensate_max (bool): set maximum normalization of frames on/off
            compensate_bkg (bool): set background compensation of frames on/off
            show_residuals (bool): set if the fit residuals must be showed (if
            a valid frame is present)
            gauss_filter (tuple): 2-element tuple with the enable toggle for
            the gaussian filter and the float sigma radius of the filter
        """
        super(ImageLab, self).__init__()

        self.filename = str(filename)

        self.images = np.array([[]]) #original images in the sis file
        self.frames = np.array([[]]) #frames
        self.frames_fit = np.array([[]]) #fitted frames
        self.frames_orig = np.array([[]]) #unaltered frames
        self.average_frame_orig = np.array([[]]) #average frame of originals

        #dimensions of a single frame
        #the stored old values are needed in order to know if there has been
        #a change in the frame dimensions after a new image is opened
        self.frame_h = int(frame_h)
        self.frame_w = int(frame_w)
        self.frame_h_old = int(frame_h_old)
        self.frame_w_old = int(frame_w_old)

        self.frame_rows = int(frame_rows)
        self.frame_cols = int(frame_cols)
        self.frame_num = int(frame_num)

        #single image that collect all frames (modifield or original)
        self.all_frames_img = np.array([[]])
        self.all_frames_img_orig = np.array([[]])

        #images of the fft and the filter mask
        self.fft_mask = np.array([[]])
        self.fft_frame = np.array([[]])

        #list of filtered points
        self.filter_points_list = []

        #slices of selected areas
        self.background_slices = bkg_slice
        self.maximum_slices = max_slice
        self.main_slices = main_slice

        #controls of image normalization
        self.compensate_bkg = compensate_bkg
        self.compensate_max = compensate_max
        self.show_residuals = show_residuals

        #gaussian filter: first element is the enable, the second is the
        #guassian sigma
        self.gauss_filter = gauss_filter

        #colorscales limits for plots
        self.vlim_img = (None, None)
        self.vlim_fft = (None, None)
        self.vlim_img_orig = (None, None) #unaltered image limits

        #corrections to the colorscale limits
        self.vlim_img_corr = tuple(float(v) for v in vlim_img_corr)
        self.vlim_fft_corr = tuple(float(v) for v in vlim_fft_corr)

        #colormaps for plots
        self.cmap_img = str(cmap_img)
        self.cmap_fft = str(cmap_fft)

        #image and frames parameters
        self.sis_img_n = int(sis_img_n)
        self.skip_frames = tuple(int(s) for s in skip_frames)
        #see _set_fft_frame() for number encoding
        self.fft_frame_num = int(fft_frame_num)

        #list of indices of the effectively valid frames in the total
        #rows x cols grid
        self.valid_frames_id = tuple()
        self.effective_frame_num = 0 #number of valid frames

        #collection of fit and coordinates results on the frames
        self.results = None


    def init(self, time_delta=tuple()):
        """
        Initialize the image data. Must be called after contructor.
        All computationally important initializations goes here.

        Args:
            time_delta (tuple): list of time intervals between the frames
            to be passed to the ResultList constructor
        """

        #check the validity of parameters for a limited maximum number of times
        #in order to avoid infinite loops if _check_parameters() fails to
        #correct something
        check_iters = 0
        while not self._check_parameters():
            check_iters += 1
            self._check_parameters()
            if check_iters > 10:
                print "Warning: probable wrong image paramers left"
                break

        #initializes the valid frames
        self._calc_frames_id()

        #initialize the results list
        self.results = ResultList(time_delta, self.valid_frames_id)
        self.update_time_delta() #must be called at every time_delta update

        #effectively load images and calculate frames
        self.sis_loadimg()
        self._set_frames()

        #set the frame on the fft panel
        self._set_fft_frame()
        #initialize the fft filter mask with ones
        self.fft_mask = np.ones_like(self.fft_frame)

        #if the size of the image changed since last image opened, guess the
        #new areas of the image (the old one might be out of range
        #or no more valid)
        if self.frame_h != self.frame_h_old \
                or self.frame_w != self.frame_w_old:
            self.guess_areas()

        #update just loaded image to calculate remaining fields
        self.update()


    def update(self):
        """
        This must be called each time a changes in the image parameters
        changes: the new frames are calculated.
        """
        #if the fft frame has been invalidated, reset it
        if self.fft_frame is None:
            self._set_fft_frame()

        #update frames and store the new total image
        self._create_mask()
        self._filter_frames()
        self._gaussian_filter()
        self._compensate_frames()
        if self.show_residuals:
            #if the residuals must be showed calculate them and then take the
            #difference between frames and fits
            self._set_fit_residuals()
            frames = self.frames - self.frames_fit
        else:
            frames = self.frames
        self.all_frames_img = self._get_all_frames_img(frames)

        #if they have been invalidated, recalculate the colorbars scales
        #from the first valid frame as reference
        fr_n = self.valid_frames_id[0]
        if None in self.vlim_img:
            self.vlim_img = self._calc_frame_vlim(self.frames[fr_n])

        if None in self.vlim_img_orig:
            self.vlim_img_orig = self._calc_frame_vlim(self.frames_orig[fr_n])

        if None in self.vlim_fft:
            self.vlim_fft = self._calc_fft_vlim()

    def update_time_delta(self, time_delta=None):
        """
        Updates the time intervals between frames in the results list.

        Args:
            time_delta (tuple): optional list of times intervals between
            frames. If None a list of ones will be considered.
        """
        #the time intervals in the results must be in the same number of
        #valid frames: if not initialize with ones
        if len(self.results.time_delta) != self.effective_frame_num:
            self.results.time_delta = tuple(1.0 for n in self.valid_frames_id)

        #if intervals are given update them
        if time_delta is not None:
            if len(time_delta) == self.effective_frame_num:
                self.results.time_delta = time_delta

        #propagate changes
        self.results.update_time()


    def _check_parameters(self):
        """
        This function checks the image parameters and eventually corrects them
        before the image initiazation. It is intended to avoid runtime errors.
        Returns True if all the controls were positive, False otherwise.
        """
        message = []
        control = True

        if self.frame_rows < 1:
            self.frame_rows = 1
            message.append(["Image error: wrong frame rows, corrected to ",
                            self.frame_rows])
            control = False

        if self.frame_cols < 1:
            self.frame_cols = 1
            message.append(["Image error: wrong frame columns, corrected to ",
                            self.frame_cols])
            control = False

        if self.frame_num < 1:
            self.frame_num = 1
            message.append(["Image error: wrong frame number, corrected to ",
                            self.frame_num])
            control = False

        if self.skip_frames[0] < 0:
            self.skip_frames[0] = 0
            message.append(["Image error: wrong skipped frames, corrected to ",
                            self.skip_frames])
            control = False

        if self.skip_frames[1] < 0:
            self.skip_frames[1] = 0
            message.append(["Image error: wrong skipped frames, corrected to ",
                            self.skip_frames])
            control = False

        if self.frame_num > self.frame_rows*self.frame_cols:
            self.frame_num = self.frame_rows*self.frame_cols
            message.append(["Image error: wrong frame number, corrected to ",
                            self.frame_num])
            control = False

        if self.frame_num - self.skip_frames[0] - self.skip_frames[1] < 1:
            self.skip_frames = (0, 0)
            message.append(["Image error: wrong skipped frames, corrected to ",
                            self.skip_frames])
            control = False

        height, width = self._get_sis_dimensions(self.filename)
        height = int(height/2) #the sis contains two images
        if self.frame_h > height or self.frame_h < 1:
            self.frame_h = height
            message.append(["Image error: wrong frame height, corrected to ",
                            self.frame_h])
            control = False

        if self.frame_w > width or self.frame_w < 1:
            self.frame_w = width
            message.append(["Image error: wrong frame width, corrected to ",
                            self.frame_w])
            control = False

        if self.frame_w*self.frame_cols > width:
            self.frame_cols = int(width/self.frame_w)
            message.append(["Image error: wrong frame columns, corrected to ",
                            self.frame_cols])
            control = False

        if self.frame_h*self.frame_rows > height:
            self.frame_rows = int(height/self.frame_h)
            message.append(["Image error: wrong frame rows, corrected to ",
                            self.frame_rows])
            control = False

        for msg in message:
            print str(msg[0]) + str(msg[1])

        return control


    def _get_sis_dimensions(self, filename):
        """
        Reads from the given sis file the stored image dimension and
        returns height, width.
        The dimension is the total one, considering the two sis parts as a
        unique one. The height must therefore divided by two as int(height/2).

        Args:
            filename (string): sis filename
        """
        with open(str(filename), 'rb') as fid:
            #skip unused data
            fid.read(10)

            #get dimensions of image
            height = int(np.fromfile(fid, dtype=np.uint16, count=1))
            width = int(np.fromfile(fid, dtype=np.uint16, count=1))

        return height, width


    def sis_loadimg(self):
        """
        Loads the sis file into the image.
        """
        #total image
        img = self._sis_read(self.filename)

        img_h, _ = img.shape
        img_h = int(img_h/2) #the sis contains two images

        #split the total image into the two parts
        img0 = img[:img_h].astype(np.float)
        img1 = img[img_h:].astype(np.float)

        #TODO: verify this normalization
        img0 = np.ma.masked_where(img0 == 0, (img0-1)/6553.6)
        img1 = np.ma.masked_where(img1 == 0, (img1-1)/6553.6)

        self.images = np.array([img0, img1])


    def _sis_read(self, filename):
        """
        Low-level interaction with the sis file for reading it.
        Returns the whole image, with the two parts in the sis as a single one.

        Args:
            filename (string): sis filename
        """
        try:
            with open(str(filename), 'rb') as fid:

                #get dimensions of image
                height, width = self._get_sis_dimensions(filename)
                fid.read(4)

                #skip unused data
                fid.read(10)
                fid.read(182)

                #xoff = int(np.fromfile(fid, dtype=np.uint16, count=1))
                #yoff = int(np.fromfile(fid, dtype=np.uint16, count=1))
                fid.read(4)

                img = np.fromfile(fid, dtype=np.uint16, count=width*height)

                img.shape = (height, width)

                return img.astype(np.uint16)

        except IOError:
            print "Error: wrong filename '" + filename + "'"
            return np.zeros((10, 10))


    def sis_writeimg(self, sis_fname, sis_format="original", sis_suffix="_"):
        """
        Saves the current image to a new sis file.

        Args:
            sis_fname (string): name of the file to save
            sis_format (string): geometry of the new sis file, "original" for a
            image with the same dimensions of the original one, "unique" for a
            single cropped image with just the frames in it, "single" for
            independet sis files for each frame
            sis_suffix (string): filename to add at the end of the filename
        """
        #create a list of images with filenames to save
        #remove extension from filename and add suffix
        sis_fname = os.path.splitext(sis_fname)[0] + sis_suffix
        img_to_write = []
        fname_to_write = []
        if sis_format == "unique":
            fname_to_write.append(sis_fname)
            img_to_write.append(self.all_frames_img)
        elif sis_format == "single":
            for frm_n, frm in enumerate(self.frames):
                img_to_write.append(frm)
                fname_to_write.append(sis_fname + "_%03d"%frm_n)
        elif sis_format == "original":
            #calculate the number of pixels to pad with zeros the frame panels
            #up to the dimension of the original total image
            _shapes = zip(self.images[self.sis_img_n].shape,
                              self.all_frames_img.shape)
            diff_shape = tuple((0, n_i-n_f) for n_i, n_f in _shapes)
            img = np.pad(self.all_frames_img, diff_shape,
                         'constant', constant_values=0)

            fname_to_write.append(sis_fname)
            img_to_write.append(img)

        for img, fname in zip(img_to_write, fname_to_write):
            #normalization and write
            self._sis_write((img+1)*6553.6, fname + ".sis")


    def _sis_write(self, image, filename):
        """
        Low-level interaction with the sis file for writing it.
        Writes the whole image, with the unused part filled with zeros.

        Args:
            image (np.array): the 2d-array that must be writed after conversion
            to 16-bit unsigned-integers (must be already normalized)
            filename (string): sis filename
        """
        #keep the double-image convention for sis files, filling the unused
        #with zeros
        if self.sis_img_n == 0:
            image = np.concatenate((image, np.zeros_like(image)))
        elif self.sis_img_n == 1:
            image = np.concatenate((np.zeros_like(image), image))

        with open(str(filename), 'wb') as fid:
            fid.write('0'*10) #skip unused

            height, width = image.shape
            size = np.array([height, width], dtype=np.uint16)
            size.tofile(fid)

            fid.write('0'*182)
            fid.write('0'*4) #skip offset

            image.astype(np.uint16).tofile(fid)


    def _calc_frames_id(self):
        """
        Initializes the valid_frames_id tuple with the id of the frames in the
        total cols x rows grid that must be considered taking care of the skip
        parameters.
        """
        self.valid_frames_id = []
        for n_frame in range(self.frame_num):
            if n_frame >= self.skip_frames[0] and \
               n_frame < self.frame_num - self.skip_frames[1]:
                self.valid_frames_id.append(n_frame)

        self.valid_frames_id = tuple(self.valid_frames_id)
        self.effective_frame_num = self.frame_num \
                                   - self.skip_frames[0] - self.skip_frames[1]


    def _set_frames(self):
        """
        Splits the total image into frames according to the frames grid.
        Valid frames that are not skypped are stored into, into a 3D np.array,
        otherwise the frame is filled with zeros.
        """
        #initialize frames array and the average of the frames
        self.frames = []
        self.average_frame_orig = np.zeros((self.frame_h, self.frame_w))

        for n_frame in range(self.frame_num):
            row = int(n_frame / self.frame_cols)
            col = int(n_frame % self.frame_cols)
            fr_h = self.frame_h
            fr_w = self.frame_w

            #each nominal frame is stored, the ones that aren't valid (skipped)
            #are filled with zeros
            if n_frame in self.valid_frames_id:
                frm = self.images[self.sis_img_n][row*fr_h:(row+1)*fr_h,
                                                  col*fr_w:(col+1)*fr_w]
                self.average_frame_orig = self.average_frame_orig + frm
            else:
                frm = np.zeros((self.frame_h, self.frame_w))

            self.frames.append(frm)

        #the frames are stored into a 3d array
        #the lenght of the array is the original frame number (with skipped)
        self.frames = np.array(self.frames)
        self.frames_orig = np.copy(self.frames) #store the unaltered frames

        self.all_frames_img_orig = self._get_all_frames_img(self.frames_orig)

        #calculate the average frame from the sum
        num = float(self.effective_frame_num)
        self.average_frame_orig = self.average_frame_orig / num


    def _set_fft_frame(self):
        """
        Sets the frame to be plotted in the fft frame from self.fft_frame_num.

        If self.fft_frame_num the frame is a valid frame id, that single frame
        is selected.
        If it is -1 the average of all frames is selected.
        If it is -2 no frame is selected and a frame of ones is considered.
        """
        #if the selection is invalid set the average
        if self.fft_frame_num not in [-2, -1]+list(self.valid_frames_id):
            self.fft_frame_num = -1

        if self.fft_frame_num == -1:
            self.fft_frame = np.fft.fft2(self.average_frame_orig)
        elif self.fft_frame_num == -2:
            self.fft_frame = np.ones_like(self.frames[self.valid_frames_id[0]])
        else:
            self.fft_frame = np.fft.fft2(self.frames_orig[self.fft_frame_num])


    def _create_mask(self):
        """
        Creates the fft filter mask from self.filter_points_list and stores it
        in self.fft_mask.
        """
        #create a list of single masks for each point
        masks = []
        for point in self.filter_points_list:
            gauss = point.get_gaussian_mask(self.fft_frame.shape)
            masks = masks + list(gauss)

        #multiply all the single masks with numpy for efficiency
        self.fft_mask = np.prod(1.0 - np.array(masks), axis=0)
        #shifts coordinates between zero-centered fft and zero-at-edges
        self.fft_mask = np.fft.fftshift(self.fft_mask)


    def _get_all_frames_img(self, frames):
        """
        Creates and return a unique np.array image from the given frames.
        Non valid frames are filled with zeros.

        Args:
            frames (np.array): a list of frames (2d np.array) or a 3d np.array
            with the frames to be merged
        """
        all_frames = np.zeros((self.frame_h*self.frame_rows,
                               self.frame_w*self.frame_cols))

        for n_frame in self.valid_frames_id:
            row = int(n_frame / self.frame_cols)
            col = int(n_frame % self.frame_cols)

            all_frames[row*self.frame_h : (row+1)*self.frame_h,
                       col*self.frame_w : (col+1)*self.frame_w] = frames[n_frame]

        return all_frames


    def _set_fit_residuals(self):
        """
        Calculates the frame images from the valid fitted results.
        """
        self.frames_fit = []
        for fr_n, frame in enumerate(self.frames):
            #if the result contain a valid fitted image consider it, or take
            #an empty frame
            if fr_n in self.valid_frames_id and \
                                fr_n in self.results.frames_list():
                fit_frame = self.results.get_result(fr_n).fit_frame
                if fit_frame is None:
                    fit_frame = np.zeros_like(frame)
            else:
                fit_frame = np.zeros_like(frame)

            self.frames_fit.append(fit_frame)
        self.frames_fit = np.array(self.frames_fit)


    def _compensate_frames(self):
        """
        If the options are activated, compensate all frames in order to have
        the same background and amplitude.
        """
        new_frames = []

        for frm in self.frames:
            if self.compensate_bkg:
                #offset compensation
                ref_bkg = np.mean(frm[self.background_slices])
                frm = frm - ref_bkg
            if self.compensate_max:
                #normalization
                ref_max = np.mean(frm[self.maximum_slices])
                frm = frm / ref_max

            new_frames.append(frm)

        self.frames = np.array(new_frames)

    def _gaussian_filter(self):
        """
        Apply a gaussian filter to the frames.
        """
        if self.gauss_filter[0]:
            frames_new = []
            for frm in self.frames:
                filtered = ndimage.gaussian_filter(frm, self.gauss_filter[1])
                frames_new.append(filtered)
            self.frames = np.array(frames_new)

    def _filter_frames(self):
        """
        Filter the frames with the fft mask.
        """
        #start from original frames
        frm = np.array(self.frames_orig)
        frm = np.fft.fft2(frm)
        frm = np.multiply(self.fft_mask, frm)
        frm = np.real(np.fft.ifft2(frm)) #take just the real part
        self.frames = frm


    def _calc_frame_vlim(self, frame):
        """
        Calculates and returns the non-corrected frames colorscale limits.

        Args:
            frame (np.array): 2d array of the frame to be used as reference
            for the colorscales
        """
        #use slices for calculate the limits
        vlim = (np.mean(frame[self.background_slices]),
                np.mean(frame[self.maximum_slices]))

        #if inverted for some reason correct the limits
        if vlim[1] < vlim[0]:
            vlim = (np.min(frame), np.max(frame))
            print "Error in frames vlim: check background and maximum regions"

        return vlim


    def _calc_fft_vlim(self):
        """
        Calculates and returns the non-corrected fft colorscale limits.
        """
        tmp_img = np.fft.fftshift(np.abs(self.fft_frame))
        shape = tmp_img.shape

        #truncate the part in the middle (at zero frequency in the shifted fft)
        #because usually it is huge
        zero_range = 20
        tmp_img[int(shape[0]/2-zero_range):int(shape[0]/2+zero_range),
                int(shape[1]/2-zero_range):int(shape[1]/2+zero_range)] = 0

        vmin = 0
        vmax = np.max(tmp_img)/2.0

        #if inverted for some reason correct the limits
        if vmax < vmin:
            vmax, vmin = np.min(np.abs(self.fft_frame)), \
                         np.max(np.abs(self.fft_frame))
            print "Error while setting fft vlim: check it"

        return (vmin, vmax)


    def get_vlim_img(self, corr=None, vlim=None):
        """
        Returns the frames colorscale limits corrected with the correction
        factors.

        Args:
            corr (tuple): optional floats with the correction to the limits in
            units of the colorscale extension; if None the ones stored in the
            class will be used
            vlim (tuple): optional input float colorscale limits to be
            corrected; if None the ones stored in the class will be used
        """
        if corr is None:
            corr = self.vlim_img_corr

        if vlim is None:
            vlim = self.vlim_img

        corr = tuple(float(c) for c in corr)
        vlim = tuple(float(v) for v in vlim)

        delta_v = vlim[1] - vlim[0]

        return (vlim[0] + delta_v*corr[0], vlim[1] + delta_v*corr[1])


    def get_vlim_fft(self, corr=None):
        """
        Returns the fft colorscale limits corrected with the correction
        factors.

        Args:
            corr (tuple): optional, see self.get_vlim_img()
        """
        if corr is None:
            corr = self.vlim_fft_corr

        return self.get_vlim_img(corr=corr, vlim=self.vlim_fft)


    def get_vlim_img_orig(self, corr=None):
        """
        Returns the original frames colorscale limits corrected with the
        correction factors.

        Args:
            corr (tuple): optional, see self.get_vlim_img()
        """
        if corr is None:
            corr = self.vlim_img_corr

        return self.get_vlim_img(vlim=self.vlim_img_orig, corr=corr)


    def _get_slices_centered(self, center, radius):
        """
        Given the center of an area, return a tuple of slices with the
        selection around the center with a certain radius. If the area is out
        of the frame dimensions the slices are cut to remain in the frame.

        Args:
            center (tuple): tuple of int with the center coordinates
            radius (int): radius around the center to select
        """
        sl0 = slice(max(0, center[0] - radius),
                    min(self.frame_h, center[0] + radius))
        sl1 = slice(max(0, center[1] - radius),
                    min(self.frame_w, center[1] + radius))
        return (sl0, sl1)


    def guess_areas(self):
        """
        Try to guess reasonable max/main/min areas of the image and stores
        the slices in the class. The first valid frame is taken as reference.
        """
        #min size of the guessed areas
        radius = 10

        #add a gaussian filter to the first valid frame to remove peaks
        frame = ndimage.gaussian_filter(self.frames[self.valid_frames_id[0]],
                                        radius)
        frame = np.ma.masked_array(frame, mask=False)

        #mask the frame at the borders in order to be shure that the selection
        #inside the frame will be of the desidered size of 2*radius
        hei, wid = frame.shape
        if hei > radius*2:
            #TODO: is really centered the slicing?
            frame[0:radius, :] = np.ma.masked
            frame[-radius:, :] = np.ma.masked
        if wid > radius*2:
            frame[:, 0:radius] = np.ma.masked
            frame[:, -radius:] = np.ma.masked

        #find minimum for background
        arg_min = np.unravel_index(np.ma.argmin(frame), frame.shape)
        arg_min = tuple(int(ar) for ar in arg_min)
        self.background_slices = self._get_slices_centered(arg_min, radius)

        #find maximum for normalization
        arg_max = np.unravel_index(np.ma.argmax(frame), frame.shape)
        arg_max = tuple(int(ar) for ar in arg_max)
        self.maximum_slices = self._get_slices_centered(arg_max, radius)

        #the main area is initializated to the whole frame on default
        self.main_slices = (slice(None, None), slice(None, None))

        #invalidate the colorbars
        self.vlim_img = (None, None)
        self.vlim_img_orig = (None, None)


    def do_fit(self, fitting):
        """
        Given a fitting class perfom the fit and save the results.

        Args:
            fitting (class): an in implementation of the Fitting class
        """
        #emit fit start
        self.fit_done.emit(0.0)
        for fit_count, n_frame in enumerate(self.valid_frames_id):
            #each frame fit creates a partial fit event
            self.fit_done.emit((fit_count + 0.5) / self.effective_frame_num)

            row = int(n_frame / self.frame_cols)
            col = int(n_frame % self.frame_cols)

            #initialize the fitting class with the current frame, do the fit
            #and save results
            fit = fitting(self.frames[n_frame])
            fit.guess_par0(slc_main=self.main_slices,
                           slc_max=self.maximum_slices,
                           slc_bkg=self.background_slices)
            results = fit.fit()
            fit_frame = fit.fitted
            #the coordinates of the fitted region are coorected with the
            #position of the frame in the grid
            xy = (fit.X+col*self.frame_w, fit.Y+row*self.frame_h)

            self.results.set_fit_result(n_frame, fit.par_names,
                                        results, xy, fit_frame)

        #emit fit finished
        self.fit_done.emit(1.0)


    def save_frames_bitmap(self, filename, reorder="original",
                           hor_crop=(0,0), ver_crop=(0,0)):
        """
        Save the frames panel image to a file.

        Args:
            filename (string): name of the image file to save
            reorder (string): specify the options for the ordering of frames
            in the output bitmap (vertical, horizontal, original)
            hor_crop (tuple): tuple of ints with the number of pixel that must
            be cutted from left/right of each frame
            ver_crop (tuple): tuple of ints with the number of pixel that must
            be cutted from top/bottom of each frame
        """
        vlim = self.get_vlim_img()

        #decide how to put the frames in order
        if reorder == "vertical":
            cols = 1
            rows = self.frame_rows*self.frame_cols
        elif reorder == "horizontal":
            cols = self.frame_cols*self.frame_rows
            rows = 1
        else:
            cols = self.frame_cols
            rows = self.frame_rows

        #take care of new cropped dimensions
        #non-valid frames are filled with zeros
        frame_h_cr = self.frame_h-ver_crop[0]-ver_crop[1]
        frame_w_cr = self.frame_w-hor_crop[0]-hor_crop[1]
        all_frames = np.zeros((rows*frame_h_cr, cols*frame_w_cr))

        for n_frame in range(self.frame_num):
            row = int(n_frame / self.frame_cols)
            col = int(n_frame % self.frame_cols)

            slice_h = slice(row*self.frame_h+ver_crop[0],
                            (row+1)*self.frame_h-ver_crop[1])
            slice_w = slice(col*self.frame_w+hor_crop[0],
                            (col+1)*self.frame_w-hor_crop[1])
            frm = self.all_frames_img[slice_h, slice_w]

            row_n = int(n_frame / cols)
            col_n = int(n_frame % cols)

            all_frames[row_n*frame_h_cr:(row_n+1)*frame_h_cr,
                       col_n*frame_w_cr:(col_n+1)*frame_w_cr] = frm

        plt.imsave(str(filename), all_frames,
                   vmin=vlim[0], vmax=vlim[1],
                   cmap=self.cmap_img)


    def save_fft_bitmap(self, filename, show_mask):
        """
        Save the fft panel image to a file.

        Args:
            filename (string): name of the image file to save
        """
        vlim = self.get_vlim_fft()

        img_fft = np.array(self.fft_frame)
        if show_mask:
            img_fft *= self.fft_mask
        img_fft = np.fft.fftshift(np.abs(img_fft))

        plt.imsave(str(filename), img_fft,
                   vmin=vlim[0], vmax=vlim[1],
                   cmap=self.cmap_fft)


    def create_movie(self, filename,
                     use_frames=True,
                     integral_x=False, integral_y=False,
                     fps=1.0, fading=False, subframes=10,
                     dpi=100, writer_name="libav"):
        """
        Creates and saves a movie from the frames.

        Args:
            filename (string): name of the movie
            use_frames (bool): put the frame image in the movie
            integral_x (bool): put the frame horizontal integral in the movie
            integral_y (bool): put the frame vertical integral in the movie
            fps (float): frame per second speed
            fading (bool): set the fading on/off between frames
            subframes (int): if the fading is active, each frame is divided
            into this number of subframes, and the transition
            will use 3 frames (min 3)
            dpi (float): resolution of the movie (max 200)
            writer_name (string): name of a valid movie writer (see the code
            forlibav, ffmpeg, imagemagick)
        """
        #check if plot can be done
        if not (use_frames or integral_x or integral_y):
            print "Error: nothing to plot in the video"
            return
        if not use_frames and integral_x and integral_y:
            print "Error: wrong movie elements selections"
            return

        #frame figure width in inches
        inch_size = 5.0

        #number of subrames (that must be > 3 for the fading)
        subframes = int(max(subframes, 4))

        #avoid too big resolution that would crash the cpu
        dpi = float(min(dpi, 200))

        #get an initial frame, and the dimension from the main slices
        frame = self.frames[self.valid_frames_id[0]]
        frame_h, frame_w = frame[self.main_slices].shape
        frame_ratio = float(frame_h) / float(frame_w)

        #calculate the movie size according to what must be plotted
        if use_frames:
            fig_size = [inch_size, inch_size*frame_ratio]
            if integral_x:
                fig_size[1] = fig_size[1] + inch_size/2.0
            if integral_y:
                fig_size[0] = fig_size[0] + inch_size/2.0
        elif integral_x:
            fig_size = [inch_size, inch_size/2.0]
        elif integral_y:
            fig_size = [inch_size/2.0, inch_size]

        plt.rcParams['figure.figsize'] = fig_size

        metadata = dict(title='Atom Movie', artist='Matplotlib',
                        comment='Video of ultracold gases')

        if fading:
            #if fading is active the fps must be multiplied for the subframes
            fps_real = float(fps)*subframes
        else:
            subframes = 1
            fps_real = float(fps)

        #initialize the writer
        if writer_name == "ffmpeg":
            writer = animation.FFMpegWriter(fps=float(fps_real),
                                            metadata=metadata)
        elif writer_name == "libav":
            writer = animation.AVConvWriter(fps=float(fps_real),
                                            metadata=metadata)
        elif writer_name == "imagemagick":
            writer = animation.ImageMagickWriter(fps=float(fps_real),
                                                 metadata=metadata)
        else:
            print "Error, movie writer not supported."
            print "Supported writers: libav, ffmpeg, imagemagick"
            return
        #TODO: check GIF first frame colors

        #initialize figure and axes
        fig = plt.figure()
        ax1 = None #frame axe
        axx = None #horizontal integral axe
        axy = None #vertical integral axe
        if integral_x and not integral_y and not use_frames:
            axx = plt.subplot2grid((1, 1), (0, 0))
        elif not integral_x and integral_y and not use_frames:
            axy = plt.subplot2grid((1, 1), (0, 0))
        elif integral_x and not integral_y and use_frames:
            ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
            axx = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
        elif not integral_x and integral_y and use_frames:
            ax1 = plt.subplot2grid((1, 3), (0, 1), colspan=2)
            axy = plt.subplot2grid((1, 3), (0, 0), sharey=ax1)
        elif integral_x and integral_y and use_frames:
            ax1 = plt.subplot2grid((3, 3), (0, 1), rowspan=2, colspan=2)
            axx = plt.subplot2grid((3, 3), (2, 1), colspan=2, sharex=ax1)
            axy = plt.subplot2grid((3, 3), (0, 0), rowspan=2, sharey=ax1)
        elif not integral_x and not integral_y and use_frames:
            ax1 = plt.subplot2grid((1, 1), (0, 0))

        #TODO control better space between subplots
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1,
                            wspace=None, hspace=None)
        fig.tight_layout(pad=0.1)


        def get_frame(frames, n_frame, n_subfr):
            """
            Calculates and returns the frame given the current frame id.
            """
            if n_subfr < subframes - 3 or not fading:
                #if the subframes are at the beginning or the fading is off
                frame = frames[n_frame]
            else:
                #id of the next frame (ciclical)
                if n_frame == self.valid_frames_id[-1]:
                    n_frame_next = self.valid_frames_id[0]
                else:
                    n_frame_next = n_frame + 1

                #basic fading algoritm between frames
                if n_subfr == subframes - 3:
                    frame = 0.75*frames[n_frame] + 0.25*frames[n_frame_next]
                elif n_subfr == subframes-2:
                    frame = 0.5*frames[n_frame] + 0.5*frames[n_frame_next]
                elif n_subfr == subframes-1:
                    frame = 0.25*frames[n_frame] + 0.75*frames[n_frame_next]
            return frame


        def calc_x_int(frame):
            """
            Calculate the horizontal integral of a frame. Returns the
            coordinates and the normalized integral.

            Args:
                frame (np.array): frame to integrate
            """
            x_int = np.sum(frame[self.main_slices], axis=0)

            #indices of the main and max selection
            ind_main = self.main_slices[1].indices(frame.shape[1])
            ind_max = self.maximum_slices[1].indices(frame.shape[1])

            #slices of the frame that will be used to calculate the
            #normalization of the integral
            norm_slice = [self.main_slices[0],
                          slice(max(ind_main[0], ind_max[0]),
                                min(ind_main[1], ind_max[1]))]

            norm = np.mean(np.sum(frame[norm_slice], axis=0))

            return np.arange(frame_w), x_int/norm


        def calc_y_int(frame):
            """
            Calculate the vertical integral of a frame. Returns the
            normalized integral and the coordinates.

            Args:
                frame (np.array): frame to integrate
            """
            y_int = np.sum(frame[self.main_slices], axis=1)

            ind_main = self.main_slices[0].indices(frame.shape[0])
            ind_max = self.maximum_slices[0].indices(frame.shape[0])

            norm_slice = [slice(max(ind_main[0], ind_max[0]),
                                min(ind_main[1], ind_max[1])),
                          self.main_slices[1]]

            norm = np.mean(np.sum(frame[norm_slice], axis=1))

            return y_int/norm, np.arange(frame_h)


        #create initial plots
        if use_frames:
            ax1_plot = ax1.imshow(frame[self.main_slices],
                                  vmin=self.get_vlim_img()[0],
                                  vmax=self.get_vlim_img()[1],
                                  aspect='auto')
            ax1.axis("off")
            ax1_plot.set_cmap(self.cmap_img)
        if integral_x:
            axx_plot, = axx.plot(*calc_x_int(frame), lw=4)
            axx.axis("off")
            axx.set_xlim(0, frame_w)
            axx.set_ylim(-0.05, 1.1)
        if integral_y:
            axy_plot, = axy.plot(*calc_y_int(frame), lw=4)
            axy.axis("off")
            axy.set_ylim(0, frame_h)
            axy.set_xlim(-0.05, 1.1)

        #reconstruct frames from the total image (this is done in order
        #to include in the viedeo all effect applied to the final image)
        frames = []
        for n_frame in range(self.frame_num):
            row = int(n_frame / self.frame_cols)
            col = int(n_frame % self.frame_cols)

            frm = self.all_frames_img[row*self.frame_h : (row+1)*self.frame_h,
                                      col*self.frame_w : (col+1)*self.frame_w]
            frames.append(frm)
        frames = np.array(frames)

        #write the movie from the figure
        with writer.saving(fig, str(filename), dpi):
            #iterate over frames and subframes
            for n_frame in self.valid_frames_id:
                for n_subfr in range(subframes):
                    frame = get_frame(frames, n_frame, n_subfr)

                    #update plot data and grab instead of pltting again
                    #for efficiency
                    if use_frames:
                        ax1_plot.set_data(frame[self.main_slices])
                    if integral_x:
                        axx_plot.set_data(*calc_x_int(frame))
                    if integral_y:
                        axy_plot.set_data(*calc_y_int(frame))

                    writer.grab_frame()


#TODO: generalize this class to non-gaussian masks
class FilterPoint(object):
    """
    Abstract point for the fft filter mask.
    """
    def __init__(self, mx=0, my=0, sx=0, sy=0, A=1, e=2):
        """
        Initialize one filter point.

        Args:
            mx (float): horizontal center
            my (float): vertical center
            sx (float): horizontal sigma
            sy (float): vertical sigma
            A (float): amplitude (dept) of the filter point
            e (float): exponent of the "gaussian"-like mask
        """
        self.mx = float(mx)
        self.my = float(my)
        self.sx = float(sx)
        self.sy = float(sy)
        self.A = float(A)
        self.e = float(e)

    def set_from_coords(self, x1, y1, x2, y2, A=1, e=2):
        """
        Given initial and final coordinates, set the filter point.

        Args:
            x1 (float): first x-coordinate
            y1 (float): first y-coordinate
            x2 (float): second x-coordinate
            y2 (float): second y-coordinate
            A (float): amplitude (dept) of the filter point
            e (float): exponent of the "gaussian"-like mask
        """
        self.mx = abs(float(x1+x2))/2
        self.my = abs(float(y1+y2))/2
        self.sx = abs(float(x1-x2))
        self.sy = abs(float(y1-y2))
        self.A = float(A)
        self.e = float(e)

    def get_gaussian_mask(self, shape):
        #TODO: check
        """
        Returns a tuple of symmetric np.arrays of the filter masks relative
        to the filter point.
        """
        y = np.arange(shape[0])
        x = np.arange(shape[1])

        X, Y = np.meshgrid(x, y)

        g1 = self.A*np.exp(-((X-self.mx)**self.e/(2*self.sx**self.e) \
                           + (Y-self.my)**self.e/(2*self.sy**self.e)))
        g2 = self.A*np.exp(-((shape[1]-X-self.mx)**self.e/(2*self.sx**self.e)
                           + (Y+self.my-shape[0])**self.e/(2*self.sy**self.e)))

        return (g1, g2)

    def get_point(self):
        """
        Returns a dictionary with the parameters of the filter point.
        """
        dic = {"mx": self.mx,
               "my": self.my,
               "sx": self.sx,
               "sy": self.sy,
               "A": self.A,
               "e": self.e}
        return dic
