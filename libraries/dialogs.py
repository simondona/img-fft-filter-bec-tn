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

from PyQt4 import QtGui, QtCore

from functools import partial



class MyGroupBox(QtGui.QGroupBox, object):
    """
    This class creates a custom GroupBox that helps the serial addition of
    fields to the group.
    """
    def __init__(self, title, layout, parent=None):
        """
        Initializes the GroupBox.

        Args:
            title (string): title of the group
            layout (QLayout): layout of the group
            parent (QObject): optional parent
        """
        super(MyGroupBox, self).__init__(title, parent=parent)
        self.setLayout(layout)

    def setWidgets(self, widgs=tuple(), w_args=None):
        """
        Adds a list of widgets with their options to the layout of the group

        Args:
            widgs (list): list of QWidgets to be added to the group
            w_args (list): optional parameters for the widgets to be added;
            leave it to None if no arguments are needed
        """
        for num, wi in enumerate(widgs):
            if w_args != None:
                self.layout().addWidget(wi, *w_args[num])
            else:
                self.layout().addWidget(wi)


class ParamsDialog(QtGui.QDialog, object):
    """
    General abstract class for the settings dialogs to be called statically.
    """
    def __init__(self, defaults, image, parent=None):
        """
        Initializes the dialog.

        Args:
            defaults (DefaultProgSettings): instance of the program
            default settings
            image (ImageLab): instance of the program image
            parent (QObject): optional parent of the dialog
        """
        super(ParamsDialog, self).__init__(parent)
        self.image = image
        self.defaults = defaults

        #modal dialog blocks the program execution until an answer is given
        self.setModal(True)

        self.setFocus()
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

    def _accept(self):
        """
        Must be called by the OK button event.
        """
        self.updateParams()
        self.accept()

    def updateParams(self):
        """
        Implement the parameters update in the child classes.
        """
        pass

    @classmethod
    def setParams(cls, defaults, image, parent=None):
        """
        Static method that opens a new dialog and if accepted stores the
        parameters in the setting defaults.

        Args:
            defaults (DefaultProgSettings): instance of the program
            default settings
            image (ImageLab): instance of the program image
            parent (QObject): optional parent of the dialog
        """
        #the current class (eventually the child class) is passed as the first
        #argument cls if the methos is a classmethod
        dialog = cls(defaults, image, parent)
        result = dialog.exec_()
        return result == QtGui.QDialog.Accepted


class OpenSisPopup(ParamsDialog):
    """
    Popup window for setting the image parameters before opening it.
    """
    def __init__(self, defaults, image, parent=None):
        super(OpenSisPopup, self).__init__(defaults, image, parent)

        self.setWindowTitle("Configure SIS file to open")

        layout = QtGui.QVBoxLayout(self)
        self.resize(200, 500)

        #single frame width and height
        self.frames_width_text = QtGui.QLineEdit(str(defaults.frame_width))
        self.frames_height_text = QtGui.QLineEdit(str(defaults.frame_height))
        group_dim = MyGroupBox("Frame dimensions", QtGui.QGridLayout(), self)
        group_dim.setWidgets([QtGui.QLabel("width:"), self.frames_width_text,
                              QtGui.QLabel("height:"), self.frames_height_text],
                             [(0, 0), (0, 1), (1, 0), (1, 1)])
        layout.addWidget(group_dim)

        #total number of frames, columns and rows
        self.frames_number_text = QtGui.QLineEdit(str(defaults.frame_number))
        self.frames_cols_text = QtGui.QLineEdit(str(defaults.frame_colums))
        self.frames_rows_text = QtGui.QLineEdit(str(defaults.frame_rows))
        group_num = MyGroupBox("Frames gometry", QtGui.QGridLayout(), self)
        group_num.setWidgets([QtGui.QLabel("number:"), self.frames_number_text,
                              QtGui.QLabel("columns:"), self.frames_cols_text,
                              QtGui.QLabel("rows:"), self.frames_rows_text],
                             [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)])
        layout.addWidget(group_num)

        #part of the SIS file to consider
        self.img_num_select = QtGui.QComboBox()
        self.img_num_select.addItems(["0", "1"])
        self.img_num_select.setCurrentIndex(defaults.sis_img_number)
        group_sis = MyGroupBox("SIS file", QtGui.QHBoxLayout(), self)
        group_sis.setWidgets([QtGui.QLabel("part number:"),
                              self.img_num_select])
        layout.addWidget(group_sis)

        #skip frames at beginning and at end of the series
        self.skip_begin_text = QtGui.QLineEdit(str(defaults.skip_frames[0]))
        self.skip_end_text = QtGui.QLineEdit(str(defaults.skip_frames[1]))
        group_skip = MyGroupBox("Skip frames", QtGui.QGridLayout(), self)
        group_skip.setWidgets([QtGui.QLabel("from begin:"), self.skip_begin_text,
                               QtGui.QLabel("from end:"), self.skip_end_text],
                              [(0, 0), (0, 1), (1, 0), (1, 1)])
        layout.addWidget(group_skip)

        #buttons
        ok_button = QtGui.QPushButton("OK")
        layout.addWidget(ok_button)
        cancel_button = QtGui.QPushButton("cancel")
        layout.addWidget(cancel_button)

        ok_button.clicked.connect(self._accept)
        cancel_button.clicked.connect(self.reject)


    def updateParams(self):
        """
        Sets the fields values into the settings.
        """
        self.defaults.frame_width = int(self.frames_width_text.text())
        self.defaults.frame_height = int(self.frames_height_text.text())
        self.defaults.frame_colums = int(self.frames_cols_text.text())
        self.defaults.frame_rows = int(self.frames_rows_text.text())
        self.defaults.frame_number = int(self.frames_number_text.text())
        self.defaults.sis_img_number = int(self.img_num_select.currentText())

        self.defaults.skip_frames = (int(self.skip_begin_text.text()),
                                     int(self.skip_end_text.text()))


class SaveSisPopup(ParamsDialog):
    """
    Popup window for setting the parameters of the sis file that must be saved.
    """
    def __init__(self, defaults, image, parent=None):
        super(SaveSisPopup, self).__init__(defaults, image, parent)

        self.setWindowTitle("Configure SIS file to save")

        layout = QtGui.QVBoxLayout(self)
        self.resize(200, 200)

        #select the geometry of the sis files that must be saved
        group_sel = MyGroupBox("Output format", QtGui.QGridLayout(), self)
        self.format_orig = QtGui.QRadioButton("original format")
        self.format_unique = QtGui.QRadioButton("frames to single image")
        self.format_single = QtGui.QRadioButton("separated frames")
        group_sel.setWidgets([self.format_orig,
                              self.format_unique,
                              self.format_single],
                             [(0, 0), (1, 0), (2, 0)])
        layout.addWidget(group_sel)

        #the suffix that must be added to the filenames in the case of single
        #frames saving selected
        self.fname_suffix = QtGui.QLineEdit(str(defaults.sis_save_suffix))
        layout.addWidget(self.fname_suffix)

        #set from previous selection
        #for the format encoding see ImageLab.sis_writeimg
        if defaults.sis_save_format == "original":
            self.format_orig.setChecked(True)
        elif defaults.sis_save_format == "unique":
            self.format_unique.setChecked(True)
        elif defaults.sis_save_format == "single":
            self.format_single.setChecked(True)
        else:
            self.format_orig.setChecked(True)

        #buttons
        ok_button = QtGui.QPushButton("OK")
        layout.addWidget(ok_button)
        cancel_button = QtGui.QPushButton("cancel")
        layout.addWidget(cancel_button)

        ok_button.clicked.connect(self._accept)
        cancel_button.clicked.connect(self.reject)


    def updateParams(self):
        """
        Sets the fields values into the settings.
        """
        #for the format encoding see ImageLab.sis_writeimg
        if self.format_orig.isChecked():
            self.defaults.sis_save_format = "original"
        elif self.format_unique.isChecked():
            self.defaults.sis_save_format = "unique"
        elif self.format_single.isChecked():
            self.defaults.sis_save_format = "single"
        self.defaults.sis_save_suffix = str(self.fname_suffix.text())


class SaveMoviePopup(ParamsDialog):
    """
    Popup window for setting the movie parameters while saving it.
    """
    def __init__(self, defaults, image, parent=None):
        super(SaveMoviePopup, self).__init__(defaults, image, parent)

        self.setWindowTitle("Configure movie")
        self.resize(200, 400)

        layout = QtGui.QGridLayout(self)


        #fps part
        self.movie_fps_text = QtGui.QLineEdit(str(defaults.movie_fps))

        group_speed = MyGroupBox("Speed (FPS)", QtGui.QHBoxLayout(), self)
        group_speed.setWidgets([self.movie_fps_text])
        layout.addWidget(group_speed, 0, 0, 1, 2)


        #plotting selections part
        self.movie_frames_box = QtGui.QCheckBox("bitmap frames")
        self.movie_frames_box.setChecked(defaults.movie_frames)

        self.movie_integral_x_box = QtGui.QCheckBox("x-integral")
        self.movie_integral_x_box.setChecked(defaults.movie_integrate_x)

        self.movie_integral_y_box = QtGui.QCheckBox("y-integral")
        self.movie_integral_y_box.setChecked(defaults.movie_integrate_y)

        group_plot = MyGroupBox("Plot", QtGui.QVBoxLayout(), self)
        group_plot.setWidgets([self.movie_frames_box,
                               self.movie_integral_x_box,
                               self.movie_integral_y_box])
        layout.addWidget(group_plot, 2, 0, 3, 2)


        #fading part: each frames is multiplied in n subframes, 3 of them
        #are used for the smooth transition
        self.movie_subframes_text = QtGui.QLineEdit(str(defaults.movie_subframes))

        self.movie_fading_box = MyGroupBox("Fading", QtGui.QHBoxLayout(), self)
        self.movie_fading_box.setWidgets([QtGui.QLabel("subframes:"),
                                          self.movie_subframes_text])
        self.movie_fading_box.setCheckable(True)
        self.movie_fading_box.setChecked(defaults.movie_fading)
        layout.addWidget(self.movie_fading_box, 5, 0, 2, 2)


        #dpi resolution part
        self.movie_dpi_text = QtGui.QLineEdit(str(defaults.movie_dpi))
        group_dpi = MyGroupBox("Resolution (DPI)", QtGui.QVBoxLayout(), self)
        group_dpi.setWidgets([self.movie_dpi_text,
                              QtGui.QLabel("(max=200dpi, width=5in)")])
        layout.addWidget(group_dpi, 7, 0, 2, 0)


        #movie-writer (encoder) part
        self.movie_writer_choice = QtGui.QComboBox()
        self.movie_writer_choice.addItem("LibAV", "libav")
        self.movie_writer_choice.addItem("FFmpeg", "ffmpeg")
        self.movie_writer_choice.addItem("ImageMagick", "imagemagick")
        choice = self.movie_writer_choice.findData(defaults.movie_writer)
        self.movie_writer_choice.setCurrentIndex(choice)
        group_writer = MyGroupBox("Movie writer", QtGui.QVBoxLayout(), self)
        group_writer.setWidgets([self.movie_writer_choice])
        layout.addWidget(group_writer, 9, 0, 2, 0)


        #buttons
        ok_button = QtGui.QPushButton("OK")
        layout.addWidget(ok_button, 11, 0, 1, 2)
        cancel_button = QtGui.QPushButton("cancel")
        cancel_button.clicked.connect(self.close)
        layout.addWidget(cancel_button, 12, 0, 1, 2)

        ok_button.clicked.connect(self._accept)
        cancel_button.clicked.connect(self.reject)


    def updateParams(self):
        """
        Sets the fields values into the settings.
        """
        self.defaults.movie_fps = float(self.movie_fps_text.text())
        subframes = self.movie_subframes_text.text()
        self.defaults.movie_subframes = int(float(subframes))
        self.defaults.movie_dpi = float(self.movie_dpi_text.text())
        self.defaults.movie_integrate_x = self.movie_integral_x_box.isChecked()
        self.defaults.movie_integrate_y = self.movie_integral_y_box.isChecked()
        self.defaults.movie_fading = self.movie_fading_box.isChecked()
        self.defaults.movie_frames = self.movie_frames_box.isChecked()

        sel = self.movie_writer_choice.currentIndex()
        self.defaults.movie_writer = self.movie_writer_choice.itemData(sel)


class SaveBitmapPopup(ParamsDialog):
    """
    Popup window for setting the movie parameters while saving it.
    """
    def __init__(self, defaults, image, parent=None):
        super(SaveBitmapPopup, self).__init__(defaults, image, parent)

        self.setWindowTitle("Configure bitmap")
        self.resize(200, 200)

        layout = QtGui.QGridLayout(self)

        #single frame crop
        self.bitmap_h_crop_texts = [QtGui.QLineEdit(str(defaults.bitmap_hor_crop[0])),
                                    QtGui.QLineEdit(str(defaults.bitmap_hor_crop[1]))]
        self.bitmap_v_crop_texts = [QtGui.QLineEdit(str(defaults.bitmap_ver_crop[0])),
                                   QtGui.QLineEdit(str(defaults.bitmap_ver_crop[1]))]
        group_crop = MyGroupBox("Frame cropping", QtGui.QGridLayout(), self)
        group_crop.setWidgets([QtGui.QLabel("hozontal:"),
                               self.bitmap_h_crop_texts[0],
                               self.bitmap_h_crop_texts[1],
                               QtGui.QLabel("vertical:"),
                               self.bitmap_v_crop_texts[0],
                               self.bitmap_v_crop_texts[1]],
                               [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)])
        layout.addWidget(group_crop)

        #set reordering options of frames
        #for the reorder encoding see ImageLab.save_frames_bitmap
        items = ["original", "vertical", "horizontal"]
        self.bitmap_reorder = QtGui.QComboBox()
        self.bitmap_reorder.addItems(items)
        if defaults.bitmap_reorder in items:
            selected = items.index(defaults.bitmap_reorder)
        else:
            selected = 0
        self.bitmap_reorder.setCurrentIndex(selected)
        group_order = MyGroupBox("Order frames", QtGui.QHBoxLayout(), self)
        group_order.setWidgets([QtGui.QLabel("display:"),
                                self.bitmap_reorder])
        layout.addWidget(group_order)

        #buttons
        ok_button = QtGui.QPushButton("OK")
        layout.addWidget(ok_button, 11, 0, 1, 2)
        cancel_button = QtGui.QPushButton("cancel")
        cancel_button.clicked.connect(self.close)
        layout.addWidget(cancel_button, 12, 0, 1, 2)

        ok_button.clicked.connect(self._accept)
        cancel_button.clicked.connect(self.reject)


    def updateParams(self):
        """
        Sets the fields values into the settings.
        """
        self.defaults.bitmap_reorder = str(self.bitmap_reorder.currentText())
        self.defaults.bitmap_hor_crop = (int(self.bitmap_h_crop_texts[0].text()),
                                        int(self.bitmap_h_crop_texts[1].text()))
        self.defaults.bitmap_ver_crop = (int(self.bitmap_v_crop_texts[0].text()),
                                        int(self.bitmap_v_crop_texts[1].text()))


class TimeTablePopup(ParamsDialog):
    """
    Popup window for setting the time intervals between frames.
    """
    def __init__(self, defaults, image, parent=None):
        super(TimeTablePopup, self).__init__(defaults, image, parent)

        self.setWindowTitle("Time table")

        layout = QtGui.QGridLayout(self)
        self.resize(200, 300)

        #series of fields for the time intervals
        group = QtGui.QGroupBox("Frame time intervals", self)
        grid_layout = QtGui.QGridLayout()
        group.setLayout(grid_layout)

        self.frames_list = []
        for n_id, n_frame in enumerate(self.image.valid_frames_id):
            text = QtGui.QLineEdit(str(self.image.results.time_delta[n_id]))
            text.setMinimumWidth(10)
            self.frames_list.append(text)
            grid_layout.addWidget(QtGui.QLabel(str(n_frame)+":"),
                                  int(n_id/2), 2*int(n_id%2))
            grid_layout.addWidget(text, int(n_id/2), 1+2*int(n_id%2))

            text.textEdited.connect(partial(self.update_text, n_id))

        layout.addWidget(group)

        #controls the propagation of the last time entered to all the
        #following fields
        self.propagate_check = QtGui.QCheckBox("Propagate last edit to the end")
        self.propagate_check.setCheckable(True)
        self.propagate_check.setChecked(True)
        layout.addWidget(self.propagate_check)

        #buttons
        button_ok = QtGui.QPushButton("OK")
        button_cancel = QtGui.QPushButton("cancel")
        layout.addWidget(button_ok)
        layout.addWidget(button_cancel)

        button_ok.clicked.connect(self._accept)
        button_cancel.clicked.connect(self.reject)

    def update_text(self, n_id):
        """
        Called when a value is edited. Propagates field changes if required.
        """
        if self.propagate_check.isChecked():
            value = 0
            for n0, fr in enumerate(self.frames_list):
                if n0 == n_id:
                    value = float(fr.text())

            for n0, fr in enumerate(self.frames_list):
                if n0 > n_id:
                    fr.setText(str(value))

    def updateParams(self):
        """
        Sets the fields values into the settings.
        """
        lis = []
        for fr in self.frames_list:
            lis.append(float(fr.text()))
        lis = tuple(lis)
        self.defaults.time_delta = lis
        self.image.update_time_delta(lis)


class FitResultsPopup(QtGui.QDialog, object):
    """
    Popup window to report the results of the fits and the coordinates.
    """
    def __init__(self, parent=None):
        super(FitResultsPopup, self).__init__(parent)

        self.setFocus()
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        self.setWindowTitle("Fit results")

        layout = QtGui.QVBoxLayout(self)
        self.resize(600, 300)

        #TODO convert to a real table
        self.text = QtGui.QTextEdit(self)
        self.text.setWordWrapMode(QtGui.QTextOption.NoWrap)
        layout.addWidget(self.text)

        bottomLayout = QtGui.QHBoxLayout()
        layout.addLayout(bottomLayout)

        self.button_save = QtGui.QPushButton("save")
        button_close = QtGui.QPushButton("close")
        bottomLayout.addWidget(self.button_save)
        bottomLayout.addWidget(button_close)

        button_close.clicked.connect(self.reject)
