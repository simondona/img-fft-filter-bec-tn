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

#pylint: disable-msg=E1101

PROGNAME = "MIFF - Multi Image FFT Filter"
PROGVERSION = "1.1.5"

import os, sys
#this sequence is needed in order to get path unicode characters working on
#some terminals and systems
reload(sys)
sys.setdefaultencoding("utf-8")
#change path
os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))

import argparse
from PyQt4 import QtGui, QtCore

from libraries.applicationwindow import ApplicationWindow


def main():
    """
    Run the program.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reset", action="store_true",
                        help="reset configurations")
    args = parser.parse_args()

    print "Multi Image FFT Filter"
    print "https://github.com/simondona/img-fft-filter-bec-tn"
    print "author: Simone Donadello - license: GNU GPL v3"
    print

    app = QtGui.QApplication(sys.argv)
    logo = QtGui.QPixmap("icons/splash.jpg")
    splash = QtGui.QSplashScreen(logo, QtCore.Qt.WindowStaysOnTopHint)
    splash.show()
    splash.showMessage(" ")

    win = ApplicationWindow(args)

    win.setWindowTitle("%s (%s)" % (PROGNAME, PROGVERSION))
    win.showMaximized()
    win.setFocus(True)

    splash.finish(win)

    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
