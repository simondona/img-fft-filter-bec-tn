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

from functools import partial
from PyQt4 import QtGui, QtCore

from libraries.settings import DefaultProgSettings
from libraries.application import Application



class ApplicationWindow(Application, QtGui.QMainWindow):
    """
    Main application class with GUI.
    """

    def __init__(self, args):
        """
        Initialize the main application.

        Args:
            args: the arguments given by the ArgumentParser
        """
        super(ApplicationWindow, self).__init__(args)

        #initialize the GUI
        self.init_ui()
        self.update_plots()


    def init_menu(self):
        """
        Initializes the menu bar.
        """
        #file menu
        file_menu = QtGui.QMenu('&File', self)
        file_menu.addAction('&Open SIS', self.open_sis,
                            QtCore.Qt.CTRL + QtCore.Qt.Key_O)
        file_menu.addAction('&Save SIS', self.save_sis,
                            QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        file_menu.addSeparator()
        file_menu.addAction('&Quit + Reset configuration',
                            self.quit_and_reset,
                            QtCore.Qt.CTRL +  QtCore.Qt.SHIFT + QtCore.Qt.Key_Q)
        file_menu.addAction('&Quit + Save configuration', self.quit_and_save,
                            QtCore.Qt.CTRL + QtCore.Qt.Key_Q)


        #help menu
        help_menu = QtGui.QMenu('&Help', self)
        help_menu.addAction('&About', self.about)


        #mask sub-menu
        mask_menu = QtGui.QMenu('&Mask', self)
        mask_menu.addAction('&Open mask', self.open_mask)
        mask_menu.addAction('&Save mask', self.save_mask)

        mask_menu.addSeparator()
        act = mask_menu.addAction('&Show mask', self._show_mask)
        act.setCheckable(True)
        act.setChecked(self.defaults.show_mask)

        #fft menu
        fft_menu = QtGui.QMenu('&FFT && Mask', self)
        self.fft_frame_menu = QtGui.QMenu('&Show FFT frame', self)
        fft_cmap_menu = QtGui.QMenu('&Colormap', self)
        fft_menu.addMenu(mask_menu)
        fft_menu.addMenu(self.fft_frame_menu)
        fft_menu.addSeparator()
        fft_menu.addMenu(fft_cmap_menu)
        fft_cmap_menu.setStyleSheet("* { menu-scrollable: 1 }")
        self._update_fft_frame_menu()
        fft_menu.addSeparator()
        fft_menu.addAction("Save fft bitmap", self.save_fft_bitmap)

        #TODO tooltips
        #mask_menu.setToolTip("Mask Menu")


        #images menu
        img_menu = QtGui.QMenu('&Image', self)

        act = img_menu.addAction('&Compensate background',
                                 self._compensate_bkg)
        act.setCheckable(True)
        act.setChecked(self.defaults.img_compensate_bkg)
        act = img_menu.addAction('&Normalize maximum',
                                 self._compensate_max)
        act.setCheckable(True)
        act.setChecked(self.defaults.img_compensate_max)

        img_menu.addSeparator()
        act = img_menu.addAction('&Set maximum area', self._set_max_area)
        act.setIcon(QtGui.QIcon('icons/icon_max.png'))
        act.setIconVisibleInMenu(True)

        act = img_menu.addAction('&Set background area', self._set_min_area)
        act.setIcon(QtGui.QIcon('icons/icon_bkg.png'))
        act.setIconVisibleInMenu(True)

        act = img_menu.addAction('&Set main area', self._set_main_area)
        act.setIcon(QtGui.QIcon('icons/icon_main.png'))
        act.setIconVisibleInMenu(True)

        img_menu.addAction('&Reset areas', self._guess_areas)

        img_menu.addSeparator()
        act = img_menu.addAction('&Show areas', self._show_areas)
        act.setCheckable(True)
        act.setChecked(self.defaults.show_regions)

        act = img_menu.addAction('&Show grid', self._show_grid)
        act.setCheckable(True)
        act.setChecked(self.defaults.show_grid)

        img_menu.addSeparator()
        frame_cmap_menu = QtGui.QMenu('&Colormap', self)
        #TODO: look better solutions for the long scrollable menu
        frame_cmap_menu.setStyleSheet("* { menu-scrollable: 1 }")
        img_menu.addMenu(frame_cmap_menu)


        #output menu
        output_menu = QtGui.QMenu('&Output', self)
        output_menu.addAction('&Save Movie',
                              self.save_movie)
        output_menu.addAction('&Save frames Bitmap',
                              self.save_bitmap)
        output_menu.addAction('&Save Results',
                              self.save_results)


        #fitting menu
        fit_menu = QtGui.QMenu('&Fit && Pick', self)
        fit_menu.addAction('&Gaussian fit',
                           partial(self.do_fit, fit_type="gauss"))
        fit_menu.addAction('&Thomas-Fermi fit',
                           partial(self.do_fit, fit_type="tf"))
        fit_menu.addAction('&Bimodal fit',
                           partial(self.do_fit, fit_type="bimodal"))
        fit_menu.addSeparator()

        act = fit_menu.addAction('&Show contour',
                                 self._show_contour)
        act.setCheckable(True)
        act.setChecked(self.defaults.show_contour)

        act = fit_menu.addAction('&Show fit residuals',
                                 self._show_residuals)
        act.setCheckable(True)
        act.setChecked(self.defaults.show_residuals)

        fit_menu.addSeparator()
        fit_menu.addAction('&Pick Coordinates (single)',
                           partial(self._set_pick_coords, multi=False))
        fit_menu.addAction('&Pick Coordinates (multi)',
                           partial(self._set_pick_coords, multi=True))
        fit_menu.addSeparator()
        fit_menu.addAction('&Show Results', self.show_results)
        fit_menu.addSeparator()
        fit_menu.addAction('&Set timing', self.show_time_table)

        fft_cmaps = QtGui.QActionGroup(self)
        frame_cmaps = QtGui.QActionGroup(self)
        for cmap in self._cmaps:
            act = QtGui.QAction(cmap, self)
            act.triggered.connect(partial(self._set_cmap_fft, cmap))
            act.setCheckable(True)
            if cmap == self.defaults.fft_cmap:
                act.setChecked(True)
            fft_cmaps.addAction(act)

            act = QtGui.QAction(cmap, self)
            act.triggered.connect(partial(self._set_cmap_frame, cmap))
            act.setCheckable(True)
            if cmap == self.defaults.frame_cmap:
                act.setChecked(True)
            frame_cmaps.addAction(act)

        fft_cmap_menu.addActions(fft_cmaps.actions())
        frame_cmap_menu.addActions(frame_cmaps.actions())


        #create the menu bar
        menubar = self.menuBar()
        menubar.addMenu(file_menu)
        menubar.addSeparator()
        menubar.addMenu(fft_menu)
        menubar.addMenu(img_menu)
        menubar.addMenu(fit_menu)
        menubar.addSeparator()
        menubar.addMenu(output_menu)
        menubar.addSeparator()
        menubar.addMenu(help_menu)


    def init_ui(self):
        """
        Initializes the GUI.
        """
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowIcon(QtGui.QIcon('icons/icon.jpg'))

        self.init_menu()


        #create panels
        main_widget = QtGui.QWidget(self)

        main_layout = QtGui.QHBoxLayout()
        main_widget.setLayout(main_layout)

        right_panel = QtGui.QVBoxLayout()
        right_big_panel = QtGui.QHBoxLayout()
        left_panel = QtGui.QVBoxLayout()
        left_big_panel = QtGui.QHBoxLayout()

        frame_tools_panel = QtGui.QHBoxLayout()

        self.frames_panel.setParent(main_widget)
        self.filter_panel.setParent(main_widget)


        #tool buttons
        frame_tools_panel.addWidget(QtGui.QLabel("Frames:"))
        create_movie_button = QtGui.QPushButton("Create Movie")
        create_movie_button.clicked.connect(self.save_movie)
        create_movie_button.setMinimumWidth(120)
        frame_tools_panel.addWidget(create_movie_button)

        save_bitmap_button = QtGui.QPushButton("Save Bitmap")
        save_bitmap_button.clicked.connect(self.save_bitmap)
        save_bitmap_button.setMinimumWidth(120)
        frame_tools_panel.addWidget(save_bitmap_button)
        frame_tools_panel.insertStretch(-1, 0)


        #buttons for the control of fft filter
        filter_tools_panel = QtGui.QHBoxLayout()
        reset_button = QtGui.QPushButton("RESET")
        reset_button.setMinimumWidth(120)
        cancel_last_button = QtGui.QPushButton("Cancel Last")
        cancel_last_button.setMinimumWidth(120)
        filter_tools_panel.addWidget(QtGui.QLabel("Mask points:"))
        filter_tools_panel.addWidget(cancel_last_button)
        filter_tools_panel.addWidget(reset_button)
        gauss_fil_sigma = QtGui.QLineEdit(str(self.image.gauss_filter[1]))
        gauss_fil_sigma.setFixedWidth(50)
        gauss_fil_sigma.setEnabled(self.image.gauss_filter[0])
        gauss_fil_box = QtGui.QCheckBox("Gauss filter")
        gauss_fil_box.setChecked(self.image.gauss_filter[0])
        filter_tools_panel.addWidget(gauss_fil_box)
        filter_tools_panel.addWidget(gauss_fil_sigma)
        filter_tools_panel.insertStretch(-1, 0)

        gauss_fil_box.toggled.connect(partial(self._set_gauss_toggle,
                                              gauss_fil_sigma))
        gauss_fil_sigma.editingFinished.connect(partial(self._set_gauss_sigma,
                                                        gauss_fil_sigma))
        reset_button.clicked.connect(self._filter_reset)
        cancel_last_button.clicked.connect(self._filter_cancel_last)


        #the images navigation toobars
        self.toolbar_filter.hide()
        self.toolbar_filter.setParent(main_widget)

        self.toolbar_frame.hide()
        self.toolbar_frame.setParent(main_widget)

        #activate the status bar
        self.status_bar = self.statusBar()


        #populate windows dictionary
        self.progress_bar = QtGui.QProgressBar(self.status_bar)
        self.status_bar.addPermanentWidget(self.progress_bar)
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(False)

        self.fft_frame_menu.setParent(self)

        #add the panels
        left_panel.addLayout(filter_tools_panel)
        left_panel.addWidget(self.filter_panel)
        left_panel.addLayout(self.nav_panel_filter)
        left_big_panel.addLayout(self.vlim_fft_panel)
        left_big_panel.addLayout(left_panel)

        right_panel.addLayout(frame_tools_panel)
        right_panel.addWidget(self.frames_panel)
        right_panel.addLayout(self.nav_panel_frames)
        right_big_panel.addLayout(right_panel)
        right_big_panel.addLayout(self.vlim_frame_panel)

        main_layout.addLayout(left_big_panel)
        main_layout.addLayout(right_big_panel)


        #connect vlim events
        self.vlim_frame_panel.connect_events(test_handle=\
                                                    self._test_frames_vlim)
        self.vlim_fft_panel.connect_events(test_handle=self._test_fft_vlim)


        #connect panels click events
        self.filter_panel.fig.canvas.mpl_connect('button_press_event',
                                                 self._on_filter_click)
        self.filter_panel.fig.canvas.mpl_connect('button_release_event',
                                                 self._on_filter_release)

        self.frames_panel.fig.canvas.mpl_connect('button_press_event',
                                                 self._on_frames_click)
        self.frames_panel.fig.canvas.mpl_connect('button_release_event',
                                                 self._on_frames_release)


        #connect panels motion events
        self.frames_panel.fig.canvas.mpl_connect('motion_notify_event',
                                                 self._on_frames_move)
        self.filter_panel.fig.canvas.mpl_connect('motion_notify_event',
                                                 self._on_filter_move)

        self.setCentralWidget(main_widget)

        #gui initialization finished
        self.is_init = True


    def quit_and_reset(self):
        """
        Resets and save settings to default before closing the application.
        """
        self.defaults = DefaultProgSettings()
        self.close()

    def quit_and_save(self):
        """
        Save current settings before closing the application.
        """
        self.close()

    def closeEvent(self, event):
        """
        Reimplements the event called when closing the application by saving
        the current settings before.
        """
        self.defaults.save_settings()
        event.accept()

    def about(self):
        """
        The program informations dialog box.
        """
        QtGui.QMessageBox.about(self, "About",
                                """\
MIFF - Multi Image FFT Filter

https://github.com/simondona/img-fft-filter-bec-tn

FFT filter, 2D fit and video creator for '.sis' images:
* remove interference fringes via FFT filtering masks
* fit images to 2D functions
* normalize images and filter noise
* create movies from multiple frames

Uses Python, Scipy, Matplotlib, Numpy, PyQt
Video support uses Libav (or FFmpeg) and ImageMagick

Copyright (C) 2014-2016  Simone Donadello
License: GNU GPL v3""")
