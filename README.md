# Multi Image FFT Filter
This software reads, filters and renders the images in the `.sis` file format,
used in the atomic physics experiment of the ultracold gases laboratory at the University of Trento - Italy (BEC research group).

https://github.com/simondona/img-fft-filter-bec-tn


### Features:
* interactively **remove interference fringes** from the experimental `.sis` images using **FFT** filtering masks
* remove **image noise** with Gaussian filtering
* **normalize** images to their dynamic range
* change images **colormaps** and ranges
* **fit** images to 2D functions and show the residuals
* serial edit of **multiple frames** contained into a single image
* **save** edited images to files
* convert multiple frames into **videos**


## Install
```
git clone https://github.com/simondona/img-fft-filter-bec-tn
```

### Dependencies:
* Python 2.7
* PyQT 4.11
* Scipy 0.14
* Matplotlib 1.3
* Numpy 1.8
* ImageMagick 6.8
* Libav 6.11 (or FFmpeg 2.5)


### Linux (Ubuntu/Debian)
```
sudo apt-get install python-qt4 python-scipy python-numpy python-matplotlib libav-tools imagemagick
```

### Linux (Fedora)
```
sudo yum install python python-qt4 scipy numpy python-matplotlib ImageMagick ffmpeg
sudo pip install --upgrade matplotlib
```

### Windows
Download and install:
* https://www.python.org/download/releases/2.7/
* http://www.riverbankcomputing.com/software/pyqt/download
* http://www.scipy.org/install.html
* http://sourceforge.net/projects/numpy/files/
* http://matplotlib.org/downloads.html
* http://www.imagemagick.org/script/binary-releases.php#windows
* https://libav.org/download.html

Create a text file in `C:\Users\username\.matplotlib\matplotlibrc` (without extension)
where `username` is the Windows user, with the following text:
```
animation.convert_path: C:\Program Files\ImageMagick-6.9.0-Q16\convert.exe
animation.avconv_path: C:\libav\usr\bin\avconv.exe
```
(eventually update the correct paths to ImageMagick `convert.exe` and LibAV `avconv.exe`)


## Usage
Launch `./img-fft-filter.py` or `python img-fft-filter.py`.
Use `-r` argument to reset the current configuration if the startup has problems.


## Author
* [Simone Donadello](https://github.com/simondona/)


## License
This software is licensed under the [GNU General Public License version 3](http://www.gnu.org/licenses/gpl-3.0.html).

```
Copyright (C) 2014-2016  Simone Donadello

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
```
