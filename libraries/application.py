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

import os
import csv

from PyQt4 import QtGui

from libraries.imagelab import ImageLab, FilterPoint
from libraries.panels import FramesCanvas, FFTCanvas, VlimPanel
from libraries.dialogs import OpenSisPopup, SaveMoviePopup,\
                              TimeTablePopup, FitResultsPopup, SaveBitmapPopup,\
                              SaveSisPopup

from matplotlib.backends.backend_qt4agg\
     import NavigationToolbar2QT as NavigationToolbar

from libraries.settings import DefaultProgSettings
from libraries import fitting

from libraries.panels import FramesZoomPanel, FilterZoomPanel

import matplotlib.pyplot as plt

from functools import partial



class FrameStatus(object):
    """
    This class stores the status of the activities of the frames panel.
    """
    def __init__(self):
        self.pan = False
        self.min = False
        self.max = False
        self.pick = False
        self.pick_multi = False
        self.main = False

    def set_status(self, what):
        """
        Set the current status.

        Args:
            what (string): status to be set
                valid values: pan, min, max, main, pick
        """
        self.pan = what == "pan"
        self.min = what == "min"
        self.max = what == "max"
        self.main = what == "main"
        self.pick = what == "pick"
        self.pick_multi = what == "pick_multi"


class FilterStatus(object):
    """
    This class stores the status of the activities of the fft filter panel.
    """
    def __init__(self):
        self.pan = False
        self.filter = False

    def set_status(self, what):
        """
        Set the current status.

        Args:
            what (string): status to be set
                valid values: pan, filter
        """
        self.pan = what == "pan"
        self.filter = what == "filter"


#options to be passed to the open/save file dialogs
FILE_DIALOG_OPT = QtGui.QFileDialog.DontUseNativeDialog


class Application(object):
    """
    Main application class without GUI.
    """

    def __init__(self, args):
        """
        Initialize the main application.

        Args:
            args: the arguments given by the ArgumentParser
        """
        super(Application, self).__init__()

        self.is_init = False
        self.args = args

        #default settings of the program
        self.defaults = DefaultProgSettings()
        if not self.args.reset:
            self.defaults.load_settings()

        self.frame_status = FrameStatus()
        #self.frame_status.set_status("pan")

        self.filter_status = FilterStatus()
        #self.filter_status.set_status("pan")

        self._cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))

        self.image = ImageLab("")

        #open initial file
        if os.path.isfile(self.defaults.img_fname):
            self.open_sis(self.defaults.img_fname)
        else:
            self.open_sis()


        #create panels for the images
        self.frames_panel = FramesCanvas(self.image,
                                         width=5, height=4, dpi=100,
                                         show_areas=self.defaults.show_regions,
                                         show_contour=self.defaults.show_contour,
                                         show_grid=self.defaults.show_grid)
        self.filter_panel = FFTCanvas(self.image,
                                      width=5, height=4, dpi=100,
                                      show_mask=self.defaults.show_mask)


        self.vlim_frame_panel = \
                        VlimPanel(init_values=self.defaults.vlim_img_corr)
        self.vlim_fft_panel = \
                        VlimPanel(scale=100, slide_range=(-200, 200),
                                  init_values=self.defaults.vlim_fft_corr)

        self.results_window = None
        self.progress_bar = None

        self.filter_panel.set_image(self.image)
        self.frames_panel.set_image(self.image)

        if os.path.isfile(self.defaults.mask_fname):
            self.open_mask(self.defaults.mask_fname)

        self.fft_frame_menu = QtGui.QMenu()
        self.status_bar = QtGui.QStatusBar()

        self.toolbar_frame = NavigationToolbar(self.frames_panel, None)
        self.toolbar_filter = NavigationToolbar(self.filter_panel, None)

        self.nav_panel_filter = FilterZoomPanel(self._set_pan_filter,
                                                self._set_filter)
        self.nav_panel_frames = FramesZoomPanel(self._set_pan_frames,
                                                self._set_pick_coords)

        self._set_pan_frames()
        self._set_pan_filter()



    def open_sis(self, fname=""):
        """
        Open a sis image file.

        Args:
            fname (string): optional filename, if empty it will be asked
        """
        def_fname = self.image.filename
        if fname == "":
            #if no filename is give ask it
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Open SIS file',
                                                      def_fname.decode("utf-8"),
                                                      filter="SIS image (*.sis)",
                                                      options=FILE_DIALOG_OPT)
            #must convert QString to python string
            fname = str(fname)

            if fname != "":
                #if a filname was given ask details
                if OpenSisPopup.setParams(self.defaults,
                                          self.image,
                                          parent=self):
                    self._load_sis(fname)
                #otherwise if any image is currently loaded ask again
                elif self.image.filename == "":
                    self.open_sis()
            #otherwise if any image is currently loaded ask again
            elif self.image.filename == "":
                #confirm exit to avoid infinite dialogs if one doesn't want to open a file
                msg = "If you don't provide a valid file to open the program will be closed. Are you sure?"
                reply = QtGui.QMessageBox.question(self, 'Close program', msg,
                                                   QtGui.QMessageBox.Yes | QtGui.QMessageBox.No,
                                                   QtGui.QMessageBox.Yes)
                if reply == QtGui.QMessageBox.Yes:
                    exit()
                else:
                    self.open_sis()
        else:
            #if a filename was give try to open it
            self._load_sis(fname)

        self.defaults.img_fname = fname


    def _load_sis(self, fname):
        """
        Loads into the program a given sis image file.

        Args:
            fname (string): filename of the sis file
        """
        if fname != "":

            self.image = ImageLab(filename=fname,
                                  frame_w=self.defaults.frame_width,
                                  frame_h=self.defaults.frame_height,
                                  frame_w_old=self.defaults.frame_width_old,
                                  frame_h_old=self.defaults.frame_height_old,
                                  frame_cols=self.defaults.frame_colums,
                                  frame_rows=self.defaults.frame_rows,
                                  frame_num=self.defaults.frame_number,
                                  sis_img_n=self.defaults.sis_img_number,
                                  vlim_img_corr=self.defaults.vlim_img_corr,
                                  vlim_fft_corr=self.defaults.vlim_fft_corr,
                                  skip_frames=self.defaults.skip_frames,
                                  cmap_img=self.defaults.frame_cmap,
                                  cmap_fft=self.defaults.fft_cmap,
                                  fft_frame_num=self.defaults.fft_frame_n,
                                  compensate_bkg=self.defaults.img_compensate_bkg,
                                  compensate_max=self.defaults.img_compensate_max,
                                  bkg_slice=self.defaults.slices_bkg,
                                  main_slice=self.defaults.slices_main,
                                  max_slice=self.defaults.slices_max,
                                  show_residuals=self.defaults.show_residuals,
                                  gauss_filter=self.defaults.gauss_filter)

            self.image.init(time_delta=self.defaults.time_delta)

            #update some defaults values because they may be corrected during
            #the ImageLab initialization
            self.defaults.slices_bkg = self.image.background_slices
            self.defaults.slices_max = self.image.maximum_slices
            self.defaults.slices_main = self.image.main_slices
            self.defaults.fft_frame_n = self.image.fft_frame_num
            self.defaults.frame_width_old = self.defaults.frame_width
            self.defaults.frame_height_old = self.defaults.frame_height
            self.defaults.frame_height = self.image.frame_h
            self.defaults.frame_width = self.image.frame_w

            #load the filter mask file if it exists
            if os.path.isfile(self.defaults.mask_fname):
                self.open_mask(self.defaults.mask_fname)

            #only if the GUI has been already initialized uptdate frames
            #and menu
            if self.is_init:
                self.status_bar.showMessage("File: \"%s\"" % fname)
                self.filter_panel.set_image(self.image)
                self.frames_panel.set_image(self.image)
                self.update_plots()
                self._update_fft_frame_menu()


    def open_mask(self, fname=""):
        """
        Loads a mask file for the fft filter.

        Args:
            fname (string): optional filename, if empty it will be asked
        """
        if fname == "":
            def_fname = os.path.splitext(self.image.filename)[0] + ".mask"
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Open mask file',
                                                      def_fname.decode("utf-8"),
                                                      filter="Mask CSV (*.mask)",
                                                      options=FILE_DIALOG_OPT)
            fname = str(fname)
        if fname != "":
            with open(fname, "rb") as fid:
                reader = csv.DictReader(fid)
                self.image.filter_points_list = \
                                    [FilterPoint(**row) for row in reader]

            self.defaults.mask_fname = fname

            #only if the GUI has been already initialized update frames
            if self.is_init:
                self.update_plots()


    def save_sis(self):
        """
        Saves the frames to sis files.
        """
        def_fname = self.image.filename
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save sis files',
                                                  def_fname.decode("utf-8"),
                                                  filter="SIS (*.sis)",
                                                  options=FILE_DIALOG_OPT)
        fname = str(fname)
        if fname != "":
            self.defaults.sis_save_fname = fname

            if SaveSisPopup.setParams(self.defaults, self.image, parent=self):
                self.image.sis_writeimg(self.defaults.sis_save_fname,
                                        self.defaults.sis_save_format,
                                        self.defaults.sis_save_suffix)


    def save_mask(self):
        """
        Saves the mask of the current fft filter to a file.
        """
        def_fname = os.path.splitext(self.image.filename)[0] + ".mask"
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save mask file',
                                                  def_fname.decode("utf-8"),
                                                  filter="Mask CSV (*.mask)",
                                                  options=FILE_DIALOG_OPT)
        fname = str(fname)
        if fname != "":
            with open(fname, "wb") as fid:
                writer = csv.DictWriter(fid,
                                        ("mx", "my", "sx", "sy", "A", "e"))
                writer.writeheader()
                writer.writerows([p.get_point() for p in self.image.filter_points_list])

            self.defaults.mask_fname = fname


    def save_movie(self):
        """
        Save the movie of the frames to a file.
        """
        def_fname = os.path.splitext(self.image.filename)[0] + ".mp4"
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save movie',
                                                  def_fname.decode("utf-8"),
                                                  filter="MPEG (*.mp4);;GIF (*.gif)",
                                                  options=FILE_DIALOG_OPT)
        fname = str(fname)

        if fname != "":
            #check the class of the selected file extension
            test_mpeg = fname.lower().endswith((".mp4", ".avi", ".mkv"))
            test_img = fname.lower().endswith(".gif")

            if test_mpeg:
                #if the extension if of the mpeg class check in the defaults if
                #a valid encoder was already selected, otherwise set a default
                #one (this is to help systems where one of libav or ffmpeg is
                #not present)
                if self.defaults.movie_writer not in ["libav", "ffmpeg"]:
                    self.defaults.movie_writer = "libav"
            elif test_img:
                #same for the image class
                self.defaults.movie_writer = "imagemagick"
            else:
                print "Error, utput movie format not supported."
                print "Supported formats: MP4, AVI, MKV, GIF"

            #if the format was valid get the parameters and write the movie
            if test_mpeg or test_img:
                if SaveMoviePopup.setParams(self.defaults,
                                            self.image,
                                            parent=self):
                    self._write_movie(fname)


    def _write_movie(self, fname):
        """
        Writes the frames movie to a file.

        Args:
            fname (string): filename of the movie
        """
        if fname != "":
            self.image.create_movie(filename=fname,
                                    fps=self.defaults.movie_fps,
                                    integral_x=self.defaults.movie_integrate_x,
                                    integral_y=self.defaults.movie_integrate_y,
                                    fading=self.defaults.movie_fading,
                                    use_frames=self.defaults.movie_frames,
                                    subframes=self.defaults.movie_subframes,
                                    dpi=self.defaults.movie_dpi,
                                    writer_name=self.defaults.movie_writer)
            self.status_bar.showMessage("Video saved to " + fname)


    def save_bitmap(self):
        """
        Saves the current frames bitmap to a file.
        """
        def_fname = os.path.splitext(self.image.filename)[0] + ".png"
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save bitmap',
                                                  def_fname.decode("utf-8"),
                                                  filter="PNG (*.png);;JPEG (*.jpg)",
                                                  options=FILE_DIALOG_OPT)
        fname = str(fname)
        if fname != "":
            if SaveBitmapPopup.setParams(self.defaults,
                                        self.image,
                                        parent=self):
                self.image.save_frames_bitmap(filename=fname,
                                              reorder=self.defaults.bitmap_reorder,
                                              hor_crop=self.defaults.bitmap_hor_crop,
                                              ver_crop=self.defaults.bitmap_ver_crop)

    def save_fft_bitmap(self):
        """
        Saves the current fft bitmap to a file.
        """
        def_fname = os.path.splitext(self.image.filename)[0] + "-fft.png"
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save bitmap',
                                                  def_fname.decode("utf-8"),
                                                  filter="PNG (*.png);;JPEG (*.jpg)",
                                                  options=FILE_DIALOG_OPT)
        fname = str(fname)
        if fname != "":
            self.image.save_fft_bitmap(filename=fname, show_mask=self.filter_panel.show_mask)


    def update_plots(self):
        """
        Updates all plots (fft and frames).
        """
        self.filter_panel.update_figure()
        self.frames_panel.update_figure(update_img=False)
        #TODO check when the update img is necessary for both


    def do_fit(self, fit_type="gauss"):
        """
        Handles the call for a fit.

        Args:
            fit_type (string): identifies the type of the fit
                "gauss" = Gaussian fit
                "tf" = Thomas-Fermi fit
                "bimodal" = Gaussian + Thomas-Fermin fit
        """

        #set progress bar updates
        self.progress_bar.setVisible(True)
        self.image.fit_done.connect(self._on_fitupdate)

        #initialize and call the fitting class
        if fit_type == "gauss":
            fit_class = fitting.Gauss2d
        elif fit_type == "tf":
            fit_class = fitting.ThomasFermi2d
        elif fit_type == "bimodal":
            fit_class = fitting.Bimodal2d
        else:
            print "Warning: wrong fit type, using gaussian"
            fit_class = fitting.Gauss2d

        self.image.do_fit(fit_class)

        #at the end update the contours in the image and results
        if self.image.show_residuals:
            #if the fit residuals are selected the image must be updated
            #after a new fit
            self.frames_panel.update_figure(update_img=True)
        self.frames_panel.update_contour()
        self.show_results()

        #hide progress bar
        self.progress_bar.setVisible(False)


    def _on_fitupdate(self, fract):
        """
        Updates the progress bar with the fraction passed as argument.

        Args:
            fract (float): number between 0 and 1 to be set in the progress bar
        """
        self.progress_bar.setValue(100.0 * float(fract))


    def show_results(self, new_window=True):
        """
        Shows the results of the fit in a window.

        Args:
            new_window (bool): open a new window if not already open
        """
        #enters if a window is present or it can be opened
        if (new_window) or\
           (not new_window and self.results_window is not None):

            #if the window is not opened open a new one
            if self.results_window == None:
                self.results_window = FitResultsPopup(parent=self)
                self.results_window.\
                        button_save.clicked.connect(self.save_results)
                self.results_window.show()

            #build the string conatining the results table
            results = "\t".join(self.image.results.keys()) + "\r\n"
            for n_id in sorted(self.image.results.frames_list()):
                dat = self.image.results.get_data(n_id)

                #if None is present write nothing
                row = "\t".join([format(dat[k], ".4g") \
                                 if dat[k] is not None \
                                 else "" \
                                 for k in self.image.results.keys()])
                results += row + "\r\n"

            #update the text
            self.results_window.text.setText(results)

            def _delete_window():
                """
                Deletes the results window for the application internal
                proprieties.
                """
                self.results_window = None

            #must use the custom window destroy for consistency
            self.results_window.destroyed.connect(_delete_window)


    def show_time_table(self):
        """
        Shows the configuration popup for setting the time intervals between
        frames.
        """
        if TimeTablePopup.setParams(self.defaults,
                                    self.image,
                                    parent=self):
            #update results window if needed
            self.show_results(new_window=False)


    def save_results(self):
        """
        Save the results of fit and coordinates pick to a file.
        """
        #check if there is something to save
        if len(self.image.results.frames_list()) > 0:
            def_fname = os.path.splitext(self.image.filename)[0] + ".csv"

            fname = QtGui.QFileDialog.getSaveFileName(self, 'Save fit results',
                                                      def_fname.decode("utf-8"),
                                                      filter="CSV (*.csv)",
                                                      options=FILE_DIALOG_OPT)
            fname = str(fname)

            if fname != "":
                with open(fname, "wb") as fid:
                    rows = [self.image.results.get_data(n_id) \
                            for n_id in self.image.results.frames_list()]

                    writer = csv.DictWriter(fid, self.image.results.keys())
                    writer.writeheader()
                    writer.writerows(rows)
        else:
            print "Warning: nothing to save"


    def _show_contour(self):
        """
        Called when there the user changes the fit contours visibility.
        """
        self.defaults.show_contour = not self.defaults.show_contour
        self.frames_panel.set_contour_visible(self.defaults.show_contour)

    def _show_residuals(self):
        """
        Called when there the user changes the fit residuals visibility.
        """
        self.defaults.show_residuals = not self.defaults.show_residuals
        self.image.show_residuals = self.defaults.show_residuals
        self.frames_panel.update_figure(update_img=True)

    def _show_areas(self):
        """
        Called when there the user changes the image areas visibility.
        """
        self.defaults.show_regions = not self.defaults.show_regions
        self.frames_panel.set_areas_visible(self.defaults.show_regions)

    def _show_grid(self):
        """
        Called when there the user changes the image grid visibility.
        """
        self.defaults.show_grid = not self.defaults.show_grid
        self.frames_panel.set_grid_visible(self.defaults.show_grid)

    def _show_mask(self):
        """
        Called when there the user changes the fft mask visibility.
        """
        self.filter_panel.show_mask = not self.filter_panel.show_mask
        self.defaults.show_mask = self.filter_panel.show_mask
        self.filter_panel.update_figure(update_img=False)

    def _compensate_max(self):
        """
        Called when there the user changes the image maximum normalization.
        """
        self.image.compensate_max = not self.image.compensate_max
        self.defaults.img_compensate_max = self.image.compensate_max
        #vlim must be updated
        self.image.vlim_img = (None, None)
        self.update_plots()

    def _compensate_bkg(self):
        """
        Called when there the user changes the image background compensation.
        """
        self.image.compensate_bkg = not self.image.compensate_bkg
        self.defaults.img_compensate_bkg = self.image.compensate_bkg
        #vlim must be updated
        self.image.vlim_img = (None, None)
        self.update_plots()

    def _set_cmap_fft(self, cmap):
        """
        Called when there the user changes fft color map.

        Args:
            cmap (string): color map name
        """
        self.image.cmap_fft = cmap
        self.filter_panel.update_figure(update_img=False)
        self.defaults.fft_cmap = cmap

    def _set_cmap_frame(self, cmap):
        """
        Called when there the user changes frames color map.

        Args:
            cmap (string): color map name
        """
        self.image.cmap_img = cmap
        self.frames_panel.update_figure(update_img=False)
        self.defaults.frame_cmap = cmap

    def _guess_areas(self):
        """
        Called when there the user ask for guessing the image areas.
        """
        self.image.guess_areas()
        self.frames_panel.update_areas()
        self.frames_panel.update_figure()
        self.defaults.slices_bkg = self.image.background_slices
        self.defaults.slices_max = self.image.maximum_slices
        self.defaults.slices_main = self.image.main_slices

    def _get_coord_frame(self, coords):
        """
        From the coordinates get the number of the frame, column and row.

        Args:
            coords (tuple): tuple of float coordinates
        """
        col = int(coords[0] / self.image.frame_w)
        row = int(coords[1] / self.image.frame_h)
        frame_num = row*self.image.frame_cols + col

        #return a negative number if the coordinates are negative
        if coords[0] < 0 or coords[1] < 0:
            frame_num = col = row = -1

        return (frame_num, col, row)

    def _set_min_area(self):
        """Changes the frame panel mode to the background area setting."""
        if self.frame_status.pan:
            self.toolbar_frame.pan()
        self.frame_status.set_status("min")
        self.nav_panel_frames.status.setText("set minimum")

    def _set_max_area(self):
        """Changes the frame panel mode to the maximum area setting."""
        if self.frame_status.pan:
            self.toolbar_frame.pan()
        self.frame_status.set_status("max")
        self.nav_panel_frames.status.setText("set maximum")

    def _set_main_area(self):
        """Changes the frame panel mode to the main area setting."""
        if self.frame_status.pan:
            self.toolbar_frame.pan()
        self.frame_status.set_status("main")
        self.nav_panel_frames.status.setText("set main")

    def _set_pick_coords(self, multi=False):
        """
        Changes the frame panel mode to the pick coordinates setting.
        If called while pick is active it goes in multiple selection mode.

        Args:
            multi (bool): keep the pick coordinates mode on for multiple
            selections
        """
        if self.frame_status.pan:
            self.toolbar_frame.pan()
        if self.frame_status.pick or multi:
            self.frame_status.set_status("pick_multi")
            self.nav_panel_frames.status.setText("multi pick")
        else:
            self.frame_status.set_status("pick")
            self.nav_panel_frames.status.setText("pick coordinates")

    def _set_pan_frames(self):
        """Changes the frame panel mode to pan/zoom."""
        if not self.frame_status.pan:
            self.toolbar_frame.pan()
        self.frame_status.set_status("pan")
        self.nav_panel_frames.status.setText("Pan/Zoom")

    def _set_pan_filter(self):
        """Changes the filter panel mode to pan/zoom."""
        if not self.filter_status.pan:
            self.toolbar_filter.pan()
        self.filter_status.set_status("pan")
        self.nav_panel_filter.status.setText("Pan/Zoom")

    def _set_filter(self):
        """Changes the filter panel mode to filtering mode."""
        if self.filter_status.pan:
            self.toolbar_filter.pan()
        self.filter_status.set_status("filter")
        self.nav_panel_filter.status.setText("Filter")

    def _set_gauss_sigma(self, text_control):
        """
        Changes the gaussian filter sigma width.

        Args:
            text_control (QLineEdit): the control where to take the sigma value
        """
        gauss_filter = list(self.image.gauss_filter)
        gauss_filter[1] = float(text_control.text())
        self.image.gauss_filter = tuple(gauss_filter)
        self.defaults.gauss_filter = tuple(gauss_filter)
        self.frames_panel.update_figure()

    def _set_gauss_toggle(self, text_control, toggle):
        """
        Enables the gaussian filter.

        Args:
            text_control (QLineEdit): the control where to take the sigma value
            toggle (bool): filter on/off
        """
        gauss_filter = list(self.image.gauss_filter)
        gauss_filter[0] = bool(toggle)
        self.image.gauss_filter = tuple(gauss_filter)
        self.defaults.gauss_filter = tuple(gauss_filter)
        self.frames_panel.update_figure()
        text_control.setEnabled(toggle)


    def _on_frames_click(self, event):
        """
        Called by the mouse click event on the frames panel.
        """
        #store first coordinates if not in pan mode
        if not self.frame_status.pan:
            self.frames_panel.coords0 = (event.xdata, event.ydata)
            self.frames_panel.selections["current"].set_visible(True)
        else:
            self.frames_panel.coords0 = (None, None)

        self.frames_panel.coords1 = self.frames_panel.coords0

    def _on_frames_move(self, event):
        """
        Called when the mouse moves on the frames panel.
        """
        #first and current coordinates must be valid
        if None not in self.frames_panel.coords0 + (event.xdata, event.ydata)\
                                    and not self.frame_status.pan:
            #get the frame number of the first and current coordinates
            frame0, _, _ = self._get_coord_frame(self.frames_panel.coords0)
            frame1, _, _ = self._get_coord_frame((event.xdata, event.ydata))

            #if current coordinates are on the same frame update second coords
            if frame0 == frame1 and frame0 >= 0 and frame1 >= 0:
                self.frames_panel.coords1 = (event.xdata, event.ydata)
                self.frames_panel.update_current_sel()


    def _on_frames_release(self, event):
        """
        Called by the mouse release event on the frames panel.
        """
        #first check if all the coordinates are valid
        if None not in self.frames_panel.coords0 + self.frames_panel.coords1:
            frame_num, _, _ = self._get_coord_frame(self.frames_panel.coords1)

            #get coords relative to a single frame
            x01 = [int(self.frames_panel.coords0[0] % self.image.frame_w),
                   int(self.frames_panel.coords1[0] % self.image.frame_w)]
            y01 = [int(self.frames_panel.coords0[1] % self.image.frame_h),
                   int(self.frames_panel.coords1[1] % self.image.frame_h)]
            x01.sort()
            y01.sort()

            #check if the selection is null, and in that case add a minimum
            #1 pixel size
            if x01[0] == x01[1]:
                x01[1] += 1
                print "Warning: null horizontal selection correted to one-sized"
            if y01[0] == y01[1]:
                y01[1] += 1
                print "Warning: null vertical selection correted to one-sized"

            #invalidate stored coords
            self.frames_panel.coords0 = self.frames_panel.coords1 = (None, None)

            #delete (by setting it to null size) rectangle of current selection
            if not self.frame_status.pan:
                self.frames_panel.selections["current"].set_visible(False)

            #update respective slices
            if self.frame_status.min:
                self.image.background_slices = (slice(*y01), slice(*x01))
                self.defaults.slices_bkg = self.image.background_slices

            if self.frame_status.max:
                self.image.maximum_slices = (slice(*y01), slice(*x01))
                self.defaults.slices_max = self.image.maximum_slices

            if self.frame_status.main:
                self.image.main_slices = (slice(*y01), slice(*x01))
                self.defaults.slices_main = self.image.main_slices

            #renormalize
            if self.frame_status.max or self.frame_status.min:
                #vlim must be updated
                self.image.vlim_img = (None, None)
                self.image.vlim_img_orig = (None, None)
                self.frames_panel.update_figure()

            #pick coordinates
            if self.frame_status.pick or self.frame_status.pick_multi:

                m_x = (x01[0] + x01[1]) / 2.0
                m_y = (y01[0] + y01[1]) / 2.0
                wid = abs(x01[0] - x01[1])
                hei = abs(y01[0] - y01[1])

                pick_keys = ["pick_x", "pick_y", "pick_w", "pick_h"]
                pick_values = [m_x, m_y, wid, hei]
                self.image.results.set_pick(frame_num, pick_keys,
                                            dict(zip(pick_keys, pick_values)))

                message = 'Coordinates: (%.1f, %.1f)' % (m_x, m_y)
                message += ' - Width: %d, Height: %d' % (wid, hei)
                self.status_bar.showMessage(message)

                #update results window if needed
                self.show_results(new_window=False)

            #finally restore pan mode
            if not self.frame_status.pan:
                self.frames_panel.update_areas()
                if not self.frame_status.pick_multi:
                    self._set_pan_frames()


    def _on_filter_click(self, event):
        """
        Called by the mouse click event on the filter panel.
        """
        if self.filter_status.filter:
            self.filter_panel.coords0 = (event.xdata, event.ydata)
            self.filter_panel.selections["current"].set_visible(True)
        else:
            self.filter_panel.coords0 = (None, None)
        self.filter_panel.coords1 = self.filter_panel.coords0


    def _on_filter_move(self, event):
        """
        Called when the mouse moves on the filter panel.
        """
        #first and current coordinates must be valid
        if None not in self.filter_panel.coords0 + (event.xdata, event.ydata)\
                                    and not self.filter_status.pan:

            self.filter_panel.coords1 = (event.xdata, event.ydata)
            self.filter_panel.update_current_sel()


    def _on_filter_release(self, event):
        """
        Called by the mouse release event on the filter panel.
        """
        if None not in self.filter_panel.coords0 + self.filter_panel.coords1:

            x01 = [self.filter_panel.coords0[0], self.filter_panel.coords1[0]]
            y01 = [self.filter_panel.coords0[1], self.filter_panel.coords1[1]]
            x01.sort()
            y01.sort()

            if x01[0] == x01[1]:
                print "Warning: null horizontal selection"
            if y01[0] == y01[1]:
                print "Warning: null vertical selection"

            if self.filter_status.filter:
                #reset the selection circle
                self.filter_panel.selections["current"].set_visible(False)
                self.image.filter_points_list.append(FilterPoint())
                self.image.filter_points_list[-1].\
                                set_from_coords(x01[0], y01[0], x01[1], y01[1])
                self.update_plots()

            #invalidate stored coords
            self.filter_panel.coords0 = self.filter_panel.coords1 = (None, None)


    def _test_frames_vlim(self, vlim):
        """
        Returns the validity of frames colorscale with the given correction.
        If valid the colorscale is update in the image class.

        Args:
            vlim (tuple): floats of the corrections to the colorscale
        """
        test_min, test_max = self.image.get_vlim_img(corr=vlim)
        test_min_orig, test_max_orig = self.image.get_vlim_img_orig(corr=vlim)

        test_cond = test_max > test_min and test_max_orig > test_min_orig

        if test_cond:
            self.image.vlim_img_corr = vlim
            self.defaults.vlim_img_corr = vlim
            self.frames_panel.update_figure(update_img=False)

        return test_cond


    def _test_fft_vlim(self, vlim):
        """
        Returns the validity of fft colorscale with the given correction.
        If valid the colorscale is update in the image class.

        Args:
            vlim (tuple): floats of the corrections to the colorscale
        """
        test_min, test_max = self.image.get_vlim_fft(corr=vlim)

        test_cond = test_max > test_min

        if test_cond:
            self.image.vlim_fft_corr = vlim
            self.defaults.vlim_fft_corr = vlim
            self.filter_panel.update_figure(update_img=False)

        return test_cond


    def _update_fft_frame_menu(self):
        """
        Initializes the menu that selects the frame to show in the fft filter
        panel.
        """
        self.fft_frame_menu.clear()
        fft_frames = QtGui.QActionGroup(self)

        selections = [-2, -1]+list(self.image.valid_frames_id)

        #if the default frame is no more present set the average as preselected
        if self.defaults.fft_frame_n in selections:
            active = self.defaults.fft_frame_n
        else:
            active = -1

        for frame_n in selections:
            #"no frame" and "average frame" are codified in the ImageLab class,
            #other frames are codified by their id
            if frame_n == -2:
                name = "None"
            elif frame_n == -1:
                name = "Average"
            else:
                name = str(frame_n)
            act = QtGui.QAction(name, self)
            act.triggered.connect(partial(self._select_fft_frame, frame_n))
            act.setCheckable(True)
            if frame_n == active:
                act.setChecked(True)
            fft_frames.addAction(act)

        self.fft_frame_menu.addActions(fft_frames.actions())


    def _select_fft_frame(self, num):
        """
        Selects the frame number for the fft panel visualization.

        Args:
            num (int): the number of the frame selected
        """
        num = int(num)
        if num != self.image.fft_frame_num:
            self.image.fft_frame_num = num
            self.defaults.fft_frame_n = num

            #invalidate current fft frame and renormalize
            self.image.fft_frame = None
            self.image.vlim_fft = (None, None)
            self.filter_panel.update_figure()


    def _filter_reset(self):
        """
        Resets the  fft filter deleting all the mask.
        """
        self.image.filter_points_list = []
        self.update_plots()


    def _filter_cancel_last(self):
        """
        Deletes just the last mask point inserted in the fft filter.
        """
        self.image.filter_points_list = self.image.filter_points_list[0:-1]
        self.update_plots()
