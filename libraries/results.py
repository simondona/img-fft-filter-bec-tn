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



class ResultList(object):
    """
    This class contains the list of results objects and the methods
    to modify them.
    """

    def __init__(self, time_delta, frames_id):
        """
        Initialize with an empty results list and time intervals.

        Args:
            frames_id (tuple): list of the frame ids
        """
        self.results = []

        #list of time intervals: each time it is modified update_time() must
        #be called in order to propagate the time to the results
        self.time_delta = time_delta

        self.frames_id = frames_id

    def frames_list(self):
        """
        Returns the list of frames contained in the results list.
        """
        return tuple(res.n_frame for res in self.results)

    def set_fit_result(self, n_frame, res_keys, res_pars, res_xy, res_z):
        """
        Add or update a fit result in the results list.

        Args:
            n_frame (int): frame number of the result
            res_keys (list): ordered keys of the dictionary of the fit result
            res_pars (dict): dictionary of the fit results
            res_xy (tuple): tuple of np.array with the xy coordinates of the
            fitted region
            res_z (np.array): the values of the fitted function
        """
        #if result already exists get it, otherwise create an empty result
        #and append it
        if n_frame in self.frames_list():
            res = self.results[self._get_id(n_frame)]
        else:
            res = Result(n_frame)
            self.results.append(res)
            self.results.sort()

        #set the result
        res.fit_keys = res_keys
        res.fit = res_pars
        res.fit_xy = res_xy
        res.fit_frame = res_z

        #update keys since they can be non uniform now
        self.update_results()

    def set_pick(self, n_frame, pick_keys, pick_coords):
        """
        Add or update a coordinates pick in the results list.

        Args:
            n_frame (int): frame number of the result
            pick_keys (list): ordered keys of the dictionary of the pick
            pick_coords (dict): dictionary of the picked coordinates
        """
        #if result already exists get it, otherwise create an empty result
        #and append it
        if n_frame in self.frames_list():
            res = self.results[self._get_id(n_frame)]
        else:
            res = Result(n_frame)
            self.results.append(res)
            self.results.sort()

        #set the result
        res.pick_keys = pick_keys
        res.pick = pick_coords

        #update keys since they can be non uniform now
        self.update_results()

    def update_results(self):
        """
        Iterate over results in order to get uniform results
        with the same keys and updates the frame times.
        """
        for _r1 in self.results:
            for _r2 in self.results:
                _r2.fill_key_values(pick_keys=_r1.pick_keys,
                                    fit_keys=_r1.fit_keys)

        self.update_time()

    def update_time(self):
        """
        Updates the frame time in the results from the list of time intervals
        stored in the class.
        """
        time = 0
        if len(self.frames_id) == len(self.time_delta):
            for n_l, n_frame in enumerate(self.frames_id):
                res = self.get_result(n_frame)
                if res is not None:
                    #if the frame has a result set the time
                    res.time = time
                time += self.time_delta[n_l]
        else:
            print "Error: wrong time-list format"

    def _get_id(self, n_frame):
        """
        Returns the id in the results list relative to a frame number id.
        If the frame doesn not already have a result in the list None is
        returned.

        Args:
            n_frame (int): frame id
        """
        lst = self.frames_list()
        if n_frame in lst:
            return lst.index(n_frame)
        else:
            return None

    def get_data(self, n_frame):
        """
        Returns the dictionary with all the data contained in the result
        relative to the specified frame.
        Returns None if the frame doesn't have a result.

        Args:
            n_frame (int): frame number
        """
        r_id = self._get_id(n_frame)
        if r_id is not None:
            return self.results[r_id].get_data()
        else:
            return None

    def get_result(self, n_frame):
        """
        Returns the result instance relative to the specified frame.
        Returns None if the frame doesn't have a result.

        Args:
            n_frame (int): frame number
        """
        r_id = self._get_id(n_frame)
        if r_id is not None:
            return self.results[r_id]
        else:
            return None

    def keys(self):
        """
        Returns the list of the ordered keys of the results in the list,
        assuming that they are uniform (call update_results() if not).
        """
        if len(self.results) > 0:
            return self.results[0].keys()
        else:
            return []


class Result(object):
    """
    This class contains the results (e.g. of the fitting and the picked
    coordinates) for a specific frame.
    """

    def __init__(self, n_frame):
        """
        Initialize and empty result for a specific frame number.

        Args:
            n_frame (int): frame number of the result
        """
        self.n_frame = int(n_frame)

        #the results must be initialized with None if not present,
        #and the keys with an empty list
        self.fit = None
        self.fit_keys = []
        self.fit_xy = None
        self.fit_frame = None

        self.pick = None
        self.pick_keys = []

        self.time = 0

    def __cmp__(self, other):
        """
        Compares and order two results using frame number values.
        """
        return self.n_frame.__cmp__(other.n_frame)

    def fill_key_values(self, fit_keys, pick_keys):
        """
        If the provided keys are not present in the result, they are added
        and the result is filled with None for each field. This is needed when
        different results must have uniform keys.

        Args:
            fit_keys (list): the list of keys relative to the fit results
            to be (eventually) added
            pick_keys (list): the list of keys relative to the picked
            coordinates to be (eventually) added
        """
        #the keys are updated only if the result is empty and if the keys
        #are spcified
        if self.fit is None and len(fit_keys) > 0:
            #dictionary with None in every key
            self.fit = dict(zip(fit_keys, tuple(None for _k in fit_keys)))
            self.fit_keys = fit_keys

        if self.pick is None and len(pick_keys) > 0:
            #dictionary with None in every key
            self.pick = dict(zip(pick_keys, tuple(None for _k in pick_keys)))
            self.pick_keys = pick_keys

    def keys(self):
        """
        Returns a list of the ordered keys contained in the result.
        """
        return ["frame", "time"] + self.fit_keys + self.pick_keys

    def get_data(self):
        """
        Returns a dictionary with all the results.
        """
        values = [self.n_frame, self.time] + \
                 [self.fit[k] for k in self.fit_keys] + \
                 [self.pick[k] for k in self.pick_keys]
        return dict(zip(self.keys(), values))
