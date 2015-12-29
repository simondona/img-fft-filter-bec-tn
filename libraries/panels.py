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

# pylint: disable=E1101

from PyQt4 import QtGui, QtCore
from functools import partial

from matplotlib.backends.backend_qt4agg \
                        import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import ColorConverter as color_conv
to_rgba = color_conv().to_rgba
import matplotlib as mpl
mpl.rcParams['font.size'] = 9

import numpy as np



class VlimPanel(QtGui.QVBoxLayout, object):
    """
    This layout contains the sliders controls for the colorscale
    min/max corrections.
    """
    def __init__(self,
                 scale=100,
                 slide_range=(-100, 100),
                 init_values=(0.0, 0.0),
                 parent=None):
        """
        Initializes the colorbar panel.

        Args:
            scale (int): scale factor of the sliders values, the sliders
            values will be normalized with this value
            slide_range (tuple): slider int ranges
            init_values (tuple): initial float correction values
            parent (QObject): parent of the layout
        """
        super(VlimPanel, self).__init__(parent)
        self.scale = scale

        #text controls of the min/max corrections
        self.text_max = QtGui.QLineEdit("%.2f"%init_values[1])
        self.text_max.setFixedWidth(45)
        self.text_min = QtGui.QLineEdit("%.2f"%init_values[0])
        self.text_min.setFixedWidth(45)

        #sliders of the min/max corrections
        self.sl_max = QtGui.QSlider(QtCore.Qt.Vertical)
        self.sl_max.setRange(*slide_range)
        self.sl_max.setValue(int(init_values[1]*self.scale))
        self.sl_max.setGeometry(30, 40, 150, 30)

        self.sl_min = QtGui.QSlider(QtCore.Qt.Vertical)
        self.sl_min.setRange(*slide_range)
        self.sl_min.setValue(int(init_values[0]*self.scale))
        self.sl_min.setGeometry(30, 40, 150, 30)

        #old correct values of the correction
        self.vmin_old, self.vmax_old = init_values

        self.addWidget(QtGui.QLabel("Vmax"))
        self.addWidget(self.text_max)
        self.addWidget(self.sl_max)

        self.addWidget(QtGui.QLabel("Vmin"))
        self.addWidget(self.text_min)
        self.addWidget(self.sl_min)


    def connect_events(self, test_handle):
        """
        Connects the events to the sliders and texts.

        Args:
            test_handle (function): handler to the function that test the
            validity of the passed corrections and returns a bool value
        """
        self.sl_max.valueChanged.connect(partial(self.on_vlim_slide_update,
                                                 test_handle=test_handle))
        self.sl_min.valueChanged.connect(partial(self.on_vlim_slide_update,
                                                 test_handle=test_handle))
        self.text_max.editingFinished.connect(partial(self.on_vlim_text_update,
                                                      test_handle=test_handle))
        self.text_min.editingFinished.connect(partial(self.on_vlim_text_update,
                                                      test_handle=test_handle))


    def get_vlim_text_corr(self):
        """
        Returns a tuple with the colorscale min/max corrections of the
        text controls.
        """
        return (float(self.text_min.text()),
                float(self.text_max.text()))

    def get_vlim_slider_corr(self):
        """
        Returns a tuple with the colorscale min/max corrections of the
        slider controls.
        """
        return (float(self.sl_min.value())/self.scale,
                float(self.sl_max.value())/self.scale)

    def restore_old(self):
        """
        Restores the old valid values of the colorscale corrections as the
        corrent ones.
        """
        self.sl_min.setValue(int(self.vmin_old*self.scale))
        self.sl_max.setValue(int(self.vmax_old*self.scale))
        self.text_min.setText("%.2f"%float(self.vmin_old))
        self.text_max.setText("%.2f"%float(self.vmax_old))

    def set_old(self):
        """
        Stores the current correction values as the old valid values.
        """
        self.vmin_old = float(self.text_min.text())
        self.vmax_old = float(self.text_max.text())

    def on_vlim_slide_update(self, test_handle):
        """
        Called when there is an update in the slider values: it tests the
        validity of the new values and eventually propagates changes.

        Args:
            test_handle (function): see connect_events(test_handle)
        """
        vlim = self.get_vlim_slider_corr()
        if test_handle(vlim):
            self.set_old()
            self.text_min.setText("%.2f"%float(vlim[0]))
            self.text_max.setText("%.2f"%float(vlim[1]))
        else:
            self.restore_old()

    def on_vlim_text_update(self, test_handle):
        """
        Called when there is an update in the text control values: it tests the
        validity of the new values and eventually propagates changes.

        Args:
            test_handle (function): see connect_events(test_handle)
        """
        vlim = self.get_vlim_text_corr()
        if test_handle(vlim):
            self.set_old()
            self.sl_min.setValue(int(vlim[0]*self.scale))
            self.sl_max.setValue(int(vlim[1]*self.scale))
        else:
            self.restore_old()


class StatusLabel(QtGui.QLabel, object):
    """
    Extension of QLabel with a custom setText method.
    """
    def __init__(self, text):
        super(StatusLabel, self).__init__(text)

    def setText(self, text):
        """
        This custom setText() add a string like "Status:" to the label text.
        """
        super(StatusLabel, self).setText("Status: "+text)


class FilterZoomPanel(QtGui.QHBoxLayout, object):
    """
    Layout of the panel with the zoom controls for the FFT filter panel.
    """
    def __init__(self, pan_handle, filter_handle, parent=None):
        """
        Initialize the panel.

        Args:
            pan_button (function): handle tot the pan/zoom function to be
            called by the button
            filter_handle (function): handle tot the set filter function
            to be called by the button
            parent (QObject): parent of the layout
        """
        super(FilterZoomPanel, self).__init__(parent)

        pan_button = QtGui.QPushButton("Pan/Zoom")
        pan_button.setMinimumWidth(120)
        filter_button = QtGui.QPushButton("Filter")
        filter_button.setMinimumWidth(120)
        self.addWidget(pan_button)
        self.addWidget(filter_button)

        self.status = StatusLabel("")
        self.addWidget(self.status)

        #fill of empty space at the end
        self.insertStretch(-1, 0)

        pan_button.clicked.connect(pan_handle)
        filter_button.clicked.connect(filter_handle)


class FramesZoomPanel(QtGui.QHBoxLayout, object):
    """
    Layout of the panel with the zoom controls for the frames panel.
    """
    def __init__(self, pan_handle, pick_handle, parent=None):
        """
        Initialize the panel.

        Args:
            pan_button (function): handle tot the pan/zoom function to be
            called by the button
            pick_handle (function): handle tot the pick coordinates function
            to be called by the button
            parent (QObject): parent of the layout
        """
        super(FramesZoomPanel, self).__init__(parent)

        pan_button = QtGui.QPushButton("Pan/Zoom")
        pan_button.setMinimumWidth(120)
        pick_button = QtGui.QPushButton("Pick Coordinates")
        pick_button.setMinimumWidth(120)
        self.addWidget(pan_button)
        self.addWidget(pick_button)

        self.status = StatusLabel("")
        self.addWidget(self.status)

        #fill of empty space at the end
        self.insertStretch(-1, 0)

        pan_button.clicked.connect(pan_handle)
        pick_button.clicked.connect(partial(pick_handle, multi=False))


class CrossLine2D(mpl.lines.Line2D):
    """
    This Line2D extension helps to build a 2D line from a coordinate selection
    with the possibility of specifying a "sign" that defines if the line
    crosses the selection along one diagonal or the other.
    """
    def __init__(self, xy, width, height, color, linewidth, turn):
        """
        Initializes the line.

        Args:
            xy (tuple): coordinates of the first corner of the selection
            width (int): width of the selection
            height (int): height of the selection
            color (string): color of the painted rectangle
            linewidth (float): linewidth of the line
            turn (int): the sign of this number defines which diagonal is
            this lise, one if positive, the other if negative
        """
        if turn >= 0:
            super(CrossLine2D, self).__init__((xy[1], xy[1]+width),
                                              (xy[0], xy[0]+height),
                                              color=to_rgba(color, alpha=0.5),
                                              linewidth=1)
        else:
            super(CrossLine2D, self).__init__((xy[1]+width, xy[1]),
                                              (xy[0], xy[0]+height),
                                              color=to_rgba(color, alpha=0.5),
                                              linewidth=1)
        self.turn = turn

    def set_xydata_turn(self, xy0, xy1):
        """
        Sets the line coordinates according to the selection coordinates
        taking care of the sign of the line.

        Args:
            xy0 (tuple): tuple of initial coordinates of the selection
            xy1 (tuple): tuple of final coordinates of the selection
        """
        if self.turn >= 0:
            self.set_xdata((xy0[0], xy1[0]))
            self.set_ydata((xy0[1], xy1[1]))
        else:
            self.set_xdata((xy1[0], xy0[0]))
            self.set_ydata((xy0[1], xy1[1]))


class SelectionFigure(object):
    """
    This generic class contains the figures relative to a selection on
    the panels and that must be painted. Child classes must implement
    constructor.
    """
    def __init__(self, xy=(0, 0), width=0, height=0, color="white", num_axes=1):
        """
        Initializes the selection figure. The constructor must be overridden
        in child classes in order to update correctly the figures lists.

        Args:
            xy (tuple): coordinates of the first corner of the selection
            width (int): width of the selection
            height (int): height of the selection
            color (string): color of the painted rectangle
            num_axes (int): number of axes to create multiple figures, one for
            each axes (otherwise adding the same patch to different axes
            wont' work)
        """
        #the lists of different types of figures are intended as lists of list:
        #each axes have his own list of figures
        self.rectangles = []
        self.circles = []
        self.lines = []
        for axe in range(num_axes):
            pass
            #in the implementation consider the number of axes correctly

    def set_visible(self, show):
        """
        Sets the figures of the selection visibility.

        Args:
            show (bool): True/False visibility
        """
        for rects in self.rectangles:
            for rec in rects:
                rec.set_visible(show)

        for lines in self.lines:
            for lin in lines:
                lin.set_visible(show)

        for circs in self.circles:
            for cir in circs:
                cir.set_visible(show)

    def update_from_coords(self, xy0, xy1):
        """
        Updates the figures of the selection according to the selection
        coordinates.

        Args:
            xy0 (tuple): first xy coordinates of the selection
            xy1 (tuple): second xy coordinates of the selection
        """
        for rects in self.rectangles:
            for rec in rects:
                rec.set_xy(xy0)
                rec.set_width(xy1[0] - xy0[0])
                rec.set_height(xy1[1] - xy0[1])

        for lines in self.lines:
            for lin in lines:
                lin.set_xydata_turn(xy0, xy1)

        for circs in self.circles:
            for cir in circs:
                cir.center = ((xy1[0] + xy0[0]) / 2.0, (xy1[1] + xy0[1]) / 2.0)
                cir.width = xy1[0] - xy0[0]
                cir.height = xy1[1] - xy0[1]

    def update_from_slices(self, slices, image, frame_n):
        """
        Updates the figures from slices defining the selection.

        Args:
            slices (tuple): tuple of slices for the vertical/horizontal
            selection
            image (ImageLab): the image instance of the program
            frame_n (int): the frame number relative to the slices
        """
        id0 = slices[0].indices(image.frame_h)
        id1 = slices[1].indices(image.frame_w)

        row = int(frame_n / image.frame_cols)
        col = int(frame_n % image.frame_cols)

        x0 = image.frame_w*col
        y0 = image.frame_h*row

        xy0 = (id1[0] + x0, id0[0] + y0)
        xy1 = (id1[1] + x0, id0[1] + y0)

        self.update_from_coords(xy0, xy1)


class SelectionRectangle(SelectionFigure):
    """
    This class contains the rectangles relative to a selection on the panels
    and that must be painted.
    """
    def __init__(self, xy=(0, 0), width=0, height=0,
                 color="white", num_axes=1):
        super(SelectionRectangle, self).__init__(xy, width, height,
                                                 color, num_axes)
        for axe in range(num_axes):
            #selection rectangle
            rect = mpl.patches.Rectangle(xy,
                                         width=width,
                                         height=height,
                                         facecolor=to_rgba(color, alpha=0.2),
                                         edgecolor=to_rgba(color, alpha=0.8),
                                         linewidth=1)
            #one rectangle for each frames panel (original and corrected)
            self.rectangles.append([rect])

            #crossing lines in the center of the selection rectangle
            line1 = CrossLine2D(xy, width, height, color, linewidth=1, turn=1)
            line2 = CrossLine2D(xy, width, height, color, linewidth=1, turn=-1)

            #two crossing lines for each frames panel (original and corrected)
            self.lines.append([line1, line2])


class SelectionCircle(SelectionFigure):
    """
    This class contains the circles relative to a selection on the panels
    and that must be painted.
    """
    def __init__(self, xy=(0, 0), width=0, height=0,
                 color="white", num_axes=1):
        super(SelectionCircle, self).__init__(xy, width, height,
                                              color, num_axes)
        for axe in range(num_axes):
            circle = mpl.patches.Ellipse(xy,
                                         width=width,
                                         height=height,
                                         facecolor=to_rgba(color, alpha=0.2),
                                         edgecolor=to_rgba(color, alpha=0.8),
                                         linewidth=1)
            self.circles.append([circle])


class FFTCanvas(FigureCanvas):
    """
    This extends the matplotlib figure canvas for building the fft-filter
    panel.
    """
    def __init__(self, image, parent=None,
                 width=5, height=4, dpi=100,
                 show_mask=True):
        """
        Initializes the panel.

        Args:
            image (ImageLab): the current image of the program
            parent (QtObject): parent of the panel
            width (float): matplotlib figure dimension in inches
            height (float): matplotlib figure dimension in inches
            dpi (float): matplotlib figure resolution
            show_mask (bool): if the mask should be showed or not in the plot
        """
        #matplotlib figure
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(FFTCanvas, self).__init__(self.fig)
        self.setParent(parent)

        self.setSizePolicy(QtGui.QSizePolicy.Expanding,
                           QtGui.QSizePolicy.Expanding)
        self.updateGeometry()

        self.plot = None
        self.axes = None

        #here the first and second coordinates of a selection can be stored
        self.coords0 = (None, None)
        self.coords1 = (None, None)

        self.show_mask = show_mask
        self.selections = {}

        self.image = None
        self.set_image(image)


    def set_image(self, image):
        """
        Sets and plots an image to the figure.

        Args:
            image (ImageLab): image instance of the program
        """
        self.image = image

        #prepare the axes and plots
        self.axes = self.fig.add_subplot(111)
        #We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.plot = self.axes.imshow(np.fft.fftshift(np.abs(self.image.fft_frame)),
                                     vmin=self.image.get_vlim_fft()[0],
                                     vmax=self.image.get_vlim_fft()[1],
                                     cmap=self.image.cmap_fft)

        self.axes.set_title("FFT + filter mask")


        #set the current selection figures
        self._create_selections()

        self.draw()


    def _create_selections(self):
        """
        Create the figures of the current selection.
        """
        circle = SelectionCircle(color="white", num_axes=1)
        self.selections["current"] = circle
        for circs in circle.circles:
            for cir in circs:
                self.axes.add_artist(cir)


    def update_current_sel(self):
        """
        Updates the current selection figure with the coordinates stored in
        the class.
        """
        self.selections["current"].update_from_coords(self.coords0,
                                                      self.coords1)
        self.draw()


    def update_figure(self, update_img=True):
        """
        Updates the plot (without creating a new one to be faster).

        Args:
            update_img (bool): update or not the ImageLab data (if data aren't
            changed should be False in order to be faster)
        """
        if update_img:
            self.image.update()

        if self.show_mask:
            draw = np.abs(self.image.fft_frame*self.image.fft_mask)
        else:
            draw = np.abs(self.image.fft_frame)

        self.plot.set_data(np.fft.fftshift(draw))
        self.plot.set_clim(vmin=self.image.get_vlim_fft()[0],
                           vmax=self.image.get_vlim_fft()[1])
        self.plot.set_cmap(self.image.cmap_fft)

        self.draw()



class FramesCanvas(FigureCanvas):
    """
    This extends the matplotlib figure canvas for building the frames panel.
    """
    def __init__(self, image, parent=None,
                 width=5, height=4, dpi=100,
                 show_areas=False, show_contour=False, show_grid=False):
        """
        Initializes the panel.

        Args:
            image (ImageLab): the current image of the program
            parent (QtObject): parent of the panel
            width (float): matplotlib figure dimension in inches
            height (float): matplotlib figure dimension in inches
            dpi (float): matplotlib figure resolution
            show_areas (bool): if the frames areas (main/max/background)
            should be showed or not in the plot
            show_contour (bool): if the fitted function contour lines
            should be showed or not in the plot
            show_grid (bool): if the grid separating the frames
            should be showed or not in the plot
        """
        #matplotlib figure
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(FramesCanvas, self).__init__(self.fig)
        self.setParent(parent)

        self.setSizePolicy(QtGui.QSizePolicy.Expanding,
                           QtGui.QSizePolicy.Expanding)
        self.updateGeometry()

        self.plot1 = None
        self.plot2 = None
        self.axes1 = None
        self.axes2 = None

        #here the first and second coordinates of a selection can be stored
        self.coords0 = (None, None)
        self.coords1 = (None, None)

        self.show_areas = show_areas
        self.show_contour = show_contour
        self.show_grid = show_grid

        self.selections = {}
        self.contours = []
        self.grid = []

        self.image = None
        self.set_image(image)


    def set_image(self, image):
        """
        Sets and plots an image to the figure.

        Args:
            image (ImageLab): image instance of the program
        """
        self.image = image

        #prepare the axes and plots
        self.axes1 = self.fig.add_subplot(211)
        self.axes2 = self.fig.add_subplot(212,
                                          sharex=self.axes1,
                                          sharey=self.axes1)
        #We want the axes cleared every time plot() is called
        self.axes1.hold(False)
        self.axes2.hold(False)

        self.plot1 = self.axes1.imshow(self.image.all_frames_img_orig,
                                       vmin=self.image.get_vlim_img_orig()[0],
                                       vmax=self.image.get_vlim_img_orig()[1],
                                       cmap=self.image.cmap_img)
        self.plot2 = self.axes2.imshow(self.image.all_frames_img,
                                       vmin=self.image.get_vlim_img()[0],
                                       vmax=self.image.get_vlim_img()[1],
                                       cmap=self.image.cmap_img)

        self.axes1.set_title("original frames")
        self.axes2.set_title("filtered and compensated frames")


        #the image selection areas must be created updated with actual
        #slices values
        self._create_selections()
        self.update_areas()

        #build new grid
        self._create_grid()

        #reset contours (maybe not needed with the current set_contour_visible
        #implementation)
        #TODO: check if it is possible to dinamically update contour lines as
        #others figures
        self._clear_contour()

        #set visibility
        self.set_areas_visible(self.show_areas)
        self.set_grid_visible(self.show_grid)
        self.set_contour_visible(self.show_contour)

        self.draw()


    def _clear_contour(self):
        """
        Initializes the contour lines, clearing the list of lines and plotted
        elements
        """
        for coll in self.contours[:]:
            if coll in self.axes2.collections:
                self.axes2.collections.remove(coll)
            self.contours.remove(coll)


    def _create_selections(self):
        """
        Creates a collection of the figures for the current/max/main/background
        selection areas.
        """
        #set the current selection figures
        rect = SelectionRectangle(color="white", num_axes=2)
        self.selections["current"] = rect
        for n_ax, axe in enumerate([self.axes1, self.axes2]):
            for rec in rect.rectangles[n_ax]:
                axe.add_artist(rec)
            for lin in rect.lines[n_ax]:
                axe.add_artist(lin)

        #set the max/main/backgroud areas
        for name, color in [["maximum", "red"],
                            ["main", "green"],
                            ["background", "blue"]]:
            #selections is a dictionary of dictionaries structured as:
            #[select type][frame num]
            self.selections[name] = {}

            #each frame has its max/main/background areas
            for n_ind in self.image.valid_frames_id:

                rect = SelectionRectangle(color=color, num_axes=2)
                self.selections[name][n_ind] = rect
                for n_ax, axe in enumerate([self.axes1, self.axes2]):
                    for rec in rect.rectangles[n_ax]:
                        axe.add_artist(rec)


    def _create_grid(self):
        """
        Creates a collection of the lines for the frames grid.
        """
        self.grid = []
        totsize = self.image.all_frames_img_orig.shape
        framesize = (self.image.frame_h, self.image.frame_w)

        for ax in [self.axes1, self.axes2]:
            y = [(n_id + 1) * framesize[0] \
                        for n_id in range(self.image.frame_rows-1)]
            line = ax.hlines(y, 0, totsize[1], colors="yellow")
            self.grid.append(line)

            x = [(n_id + 1) * framesize[1] \
                        for n_id in range(self.image.frame_cols-1)]
            line = ax.vlines(x, 0, totsize[0], colors="yellow")
            self.grid.append(line)


    def set_areas_visible(self, show):
        """
        Sets the main/max/background areas figures visibility.

        Args:
            show (bool): areas visibility
        """
        self.show_areas = show
        for n_ind in self.image.valid_frames_id:
            self.selections["main"][n_ind].set_visible(show)
            self.selections["background"][n_ind].set_visible(show)
            self.selections["maximum"][n_ind].set_visible(show)

        self.draw()


    def set_contour_visible(self, show):
        """
        Sets the fitted contours lines visibility.

        Args:
            show (bool): contour visibility
        """
        self.show_contour = show
        self.update_contour()


    def set_grid_visible(self, show):
        """
        Sets the frames grid lines visibility.

        Args:
            show (bool): grid visibility
        """
        self.show_grid = show
        for line in self.grid:
            line.set_visible(show)

        self.draw()


    def update_contour(self):
        """
        The contour lines of the fitted function are created and plotted.
        """
        if self.show_contour:
            self._clear_contour()

            #must keep hold on to plot over the figure
            self.axes2.hold(True)

            #create contours for each fit result
            for n_ind in self.image.results.frames_list():

                res = self.image.results.get_result(n_ind)
                if res.fit is not None:

                    #contour lines
                    number_lines = 4

                    #fit offset
                    if "offs" in res.keys():
                        offset = res.fit["offs"]
                    else:
                        offset = 0

                    #sum fit amplitudes (bimodal may have 2 amplitudes)
                    tot_amp = 0
                    for key in res.keys():
                        if key.startswith("amp"):
                            tot_amp += res.fit[key]

                    #calculate contour levels (and remove first-last)
                    levels = np.linspace(offset, tot_amp, number_lines+2)[1:-1]

                    contours = self.axes2.contour(res.fit_xy[0],
                                                  res.fit_xy[1],
                                                  res.fit_frame,
                                                  colors='white',
                                                  levels=levels)
                    self.contours = self.contours + contours.collections

            #restore hold off
            self.axes2.hold(False)
        else:
            self._clear_contour()

        self.draw()


    def update_areas(self):
        """
        Update max/main/background areas according to image slices.
        """
        rects = self.selections
        for n_ind in self.image.valid_frames_id:
            rects["main"][n_ind].update_from_slices(self.image.main_slices,
                                                    self.image,
                                                    n_ind)
            rects["background"][n_ind].update_from_slices(self.image.background_slices,
                                                          self.image,
                                                          n_ind)
            rects["maximum"][n_ind].update_from_slices(self.image.maximum_slices,
                                                       self.image,
                                                       n_ind)

        self.draw()


    def update_current_sel(self):
        """
        Updates the current selection figure with the coordinates stored in
        the class.
        """
        self.selections["current"].update_from_coords(self.coords0,
                                                      self.coords1)
        self.draw()

    def update_figure(self, update_img=True):
        """
        Updates the plot (without creating a new one to be faster).

        Args:
            update_img (bool): update or not the ImageLab data (if data aren't
            changed should be False in order to be faster)
        """
        if update_img:
            self.image.update()

        self.plot1.set_data(self.image.all_frames_img_orig)
        self.plot2.set_data(self.image.all_frames_img)
        self.plot1.set_clim(vmin=self.image.get_vlim_img_orig()[0],
                            vmax=self.image.get_vlim_img_orig()[1])
        self.plot2.set_clim(vmin=self.image.get_vlim_img()[0],
                            vmax=self.image.get_vlim_img()[1])
        self.plot1.set_cmap(self.image.cmap_img)
        self.plot2.set_cmap(self.image.cmap_img)

        self.draw()
