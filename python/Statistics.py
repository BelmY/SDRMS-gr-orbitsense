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
import pmt
from gnuradio import gr

class Statistics(gr.sync_block):
    """
    docstring for block Statistics
    """
    def __init__(self, samp_rate, signal_on, modulation, nf_estimation_time):
        gr.sync_block.__init__(self,
            name="Statistics",
            in_sig=[numpy.complex64],
            out_sig=[numpy.complex64])
        self.samp_rate = samp_rate
        self.signal_on = signal_on
        self.modulation = modulation
        self.nf_estimation_time = nf_estimation_time
        self.nf_estimation_samples = nf_estimation_time * samp_rate
        self.noise_samples = 0
        self.signal_samples = 0
        self.noise_power = 0
        self.noise_power_list = []
        self.signal_noise_power = 0
        self.signal_power_list = []
        self.snr = 0
        self.detection_ratio = 0
        self.false_alarm_ratio = 0.0
        self.detected = 0.0
        self.detected_signal = 0
        self.detected_no_signal = 0
        self.flag2 = 0
        self.message_port_register_in(pmt.intern("data_in"))
        self.set_msg_handler(pmt.intern("data_in"), self.handle_msg)

    def __del__(self):
        print "EKALESTHKA"

    def work(self, input_items, output_items):
        in0 = input_items[0]
        out = output_items[0]
        # <+signal processing here+>
        if not self.signal_on and self.noise_samples < self.nf_estimation_samples:
            self.noise_samples += len(input_items[0])
            self.noise_power_estimation(in0, len(input_items[0]))
        elif self.signal_on and self.signal_samples < self.nf_estimation_samples:
            self.signal_samples += len(input_items[0])
            self.power_estimation(in0, len(input_items[0]))
        
        '''
        if self.signal_on and self.detected > 0:
            print "Detection Ratio: ", self.detected_signal/self.detected
        elif self.signal_on and self.detected == 0:
            print "Detection Ratio: ", 0
        '''
        '''
        if not self.signal_on and self.detected > 0:
            print "False Alarm: ", self.detected_no_signal/self.detected
        elif not self.signal_on and self.detected == 0:
            print "False Alarm: ", 0
        '''
        if self.detected >= 1000:
            print "Detection Ratio: ", self.detected_signal/self.detected
            print "Algorithm responses: ", self.detected
            print "Detected: ", self.detected_signal
            return -1
        
        out[:] = in0
        return len(output_items[0])

    def noise_power_estimation(self, samples, input_items):
        estimated_power = numpy.sum(numpy.absolute(samples) ** 2) / input_items
        self.noise_power_list.append(estimated_power)
        if self.signal_on or self.noise_samples >= self.nf_estimation_samples:
            self.noise_power = sum(self.noise_power_list) / len(self.noise_power_list)
            print "noise power: ", self.noise_power

    def power_estimation(self, samples, input_items):
        estimated_power = numpy.sum(numpy.absolute(samples) ** 2) / input_items
        self.signal_power_list.append(estimated_power)
        if self.signal_samples >= self.nf_estimation_samples:
            self.signal_noise_power = sum(self.signal_power_list) / len(self.signal_power_list)
            print "signal - noise: ", (self.signal_noise_power - self.noise_power)
            self.snr = 10 * numpy.log10((self.signal_noise_power - self.noise_power) /
                                      self.noise_power)
            print "SNR: ", self.snr

        return estimated_power 

    def stats(self, msg):
        if pmt.to_python(msg) or self.signal_on:
            self.detected += 1
        if pmt.to_python(msg) and self.signal_on:
            self.detected_signal += 1
        if pmt.to_python(msg) and not self.signal_on:
            #print msg
            self.detected_no_signal += 1

    def handle_msg(self, msg):
        self.stats(msg)
        pass

    def set_signal(self, on):
        self.signal_on = on
        print self.signal_on
