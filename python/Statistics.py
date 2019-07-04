#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2019 gr-orbitsense author.
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

import numpy
from gnuradio import gr

class Statistics(gr.sync_block):
    """
    docstring for block Statistics
    """
    def __init__(self, samp_rate, modulation):
        gr.sync_block.__init__(self,
            name="Statistics",
            in_sig=[numpy.complex64],
            out_sig=[numpy.complex64])
        self.samp_rate = samp_rate
        self.modulation = modulation
        self.noise_power = 0
        self.signal_noise_power = 0
        self.snr = 0
        self.flag = 1

    def work(self, input_items, output_items):
        in0 = input_items[0]
        out = output_items[0]
        # <+signal processing here+>
        if self.flag == 1:
          self.noise_power = self.power_estimation(in0, len(input_items[0]))
          print "noise power: ", self.noise_power
        self.signal_noise_power = self.power_estimation(in0, len(input_items[0]))
        self.snr = 10 * numpy.log10((self.signal_noise_power - self.noise_power) /
                                    self.noise_power)
        print self.snr
        out[:] = in0
        return len(output_items[0])

    def power_estimation(self, samples, input_items):
        self.flag = 0
        return numpy.sum(numpy.absolute(samples) ** 2) / input_items

