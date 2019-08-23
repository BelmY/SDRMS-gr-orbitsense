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
import math
import os
import pandas
import matplotlib

matplotlib.use('TKAgg', warn = True, force = True)

import matplotlib.pyplot as plt

from gnuradio import gr
from gnuradio import analog

class Statistics(gr.sync_block):
    """
    docstring for block Statistics
    """
    def __init__(self, samp_rate, signal_on, signal_estimation_time, iterations,\
                 snr_list, export_csv, csv_file_path, export_plots, plot_file_path,\
                 modulation, smoothing_factor, signal_bandwidth, false_alarm, window_size):
        gr.sync_block.__init__(self,
            name="Statistics",
            in_sig=[numpy.complex64],
            out_sig=[numpy.complex64])
        # Initialize private members
        self.samp_rate = samp_rate
        self.signal_on = signal_on
        self.signal_estimation_time = signal_estimation_time
        self.iterations = iterations
        self.snr_list = snr_list
        self.export_csv = export_csv
        self.csv_file_path = csv_file_path
        self.export_plots = export_plots
        self.plot_file_path = plot_file_path
        self.modulation = modulation
        self.smoothing_factor = smoothing_factor
        self.signal_bandwidth = signal_bandwidth
        self.false_alarm = false_alarm
        self.window_size = window_size
        self.signal_estimation_samples = self.signal_estimation_time * samp_rate
        self.snr_index = 0
        self.signal_samples = 0.0
        self.signal_power = 0.0
        self.signal_power_list = []
        self.noise_power = 0.0
        self.signal_noise_power = 0.0
        self.signal_power_list = []
        #self.snr = 0.0
        self.detection_ratio = 0.0
        self.detection_ratio_list = []
        self.false_alarm_ratio = 0.0
        self.detected = 0.0
        self.detected_signal = 0.0
        self.detected_no_signal = 0.0
        # Initialize AWGN
        self.noise = analog.fastnoise_source_c(analog.GR_GAUSSIAN, 0.0, 0, 8192)
        # Initialize message port
        self.message_port_register_in(pmt.intern("data_in"))
        self.set_msg_handler(pmt.intern("data_in"), self.handle_msg)

    def work(self, input_items, output_items):
        in0 = input_items[0]
        out = output_items[0]
        # <+signal processing here+>
        out[:] = in0

        if self.signal_samples < self.signal_estimation_samples:
            self.signal_samples += len(in0)
            self.signal_power_estimation(in0, len(in0))
        else:
            if self.noise_power == 0.0:
                self.noise_power = self.calculate_noise_power(self.signal_power, self.snr_list[self.snr_index])
                self.noise.set_amplitude(self.noise_power)

            for i in range(0, len(out)):
                out[i] = in0[i] + self.noise.sample()    
                    
        if self.detected >= self.iterations:
            print "---------------------------------------------------------"
            print "SNR: ", self.snr_list[self.snr_index]
            print "Detection Ratio: ", self.detected_signal/self.detected
            print "Algorithm responses: ", self.detected
            print "Detected: ", self.detected_signal
            self.detection_ratio = self.detected_signal/self.detected
            self.detection_ratio_list.append(self.detection_ratio)
            self.noise_power = 0.0
            self.snr_index += 1
            self.detected = 0.0
            self.detected_signal = 0.0
            self.detected_no_signal = 0

        if len(self.snr_list) == self.snr_index:
            if (self.export_csv):
                self.export_csv_detection_ratio(self.snr_list, self.detection_ratio_list, self.csv_file_path)
            if (self.export_plots):
                print(matplotlib.rcParams['axes.titlesize'])
                self.plot(self.snr_list, self.detection_ratio_list, self.plot_file_path)
            return -1
        
        return len(output_items[0])

    def signal_power_estimation(self, samples, input_items):
        estimated_power = numpy.sum(numpy.absolute(samples) ** 2) / input_items
        self.signal_power_list.append(estimated_power)
        if self.signal_samples >= self.signal_estimation_samples:
            self.signal_power = sum(self.signal_power_list) / len(self.signal_power_list)

    def calculate_noise_power(self, signal_power, snr):
        noise_power = signal_power * (10.0**(-snr/10.0))
        return math.sqrt(noise_power)
    
    def stats(self, msg):
        if pmt.to_python(msg) or self.signal_on:
            self.detected += 1
        if pmt.to_python(msg) and self.signal_on:
            self.detected_signal += 1
        if pmt.to_python(msg) and not self.signal_on:
            self.detected_no_signal += 1

    def handle_msg(self, msg):
        if self.signal_samples >= self.signal_estimation_samples:
            self.stats(msg)

    def set_signal(self, on):
        self.signal_on = on

    def plot(self, snr_list, detection_ratio, path):
        matplotlib.rcParams['figure.facecolor'] = 'w'
        matplotlib.rcParams['legend.loc'] = 'upper left'
        
        fig = plt.figure(figsize=[20, 15])
        
        plt.plot(snr_list, detection_ratio, marker='o')
        
        plt.suptitle('Detection Ratio per SNR')
        title = 'sampling rate: ' + str(self.samp_rate)
        if self.smoothing_factor:
            title += ', smoothing factor: ' + str(self.smoothing_factor)
        if self.signal_bandwidth:
            title += ', signal bandwidth: ' + str(self.signal_bandwidth)
        if self.false_alarm:
            title += ', probability false alarm: ' + str(self.false_alarm)

        plt.title(title)
        
        plt.xlabel('SNR (dB)')
        plt.ylabel('Detection Ratio')
        plt.grid(b = True)
        
        plt.legend([str(self.window_size) + ' sps'])
        plt.yticks(numpy.arange(0.0, 1.1, 0.1))
        
        if (os.path.isdir(path)):
            path += '/' + 'Detection_ratio_' + str(self.modulation) + '_' + str(self.samp_rate) + 'sps_' +\
                    str(self.false_alarm) + 'prob_fa_' + str(self.window_size) + '_window' + '.jpg'
            plt.savefig(path)
        plt.show(block = True)

    def export_csv_detection_ratio(self, snr_list, detection_ratio_list, path):
        if (os.path.isdir(path)):
            path += '/' + 'Detection_ratio_' + str(self.modulation) + '_' + str(self.samp_rate) + 'sps_' +\
                    str(self.false_alarm) + 'prob_fa_' + str(self.window_size) + '_window' + '.csv'
        df1 = pandas.DataFrame(snr_list, columns = ['SNR'])
        df2 = pandas.DataFrame(detection_ratio_list, columns = ['Detection_Ratio'])

        data = pandas.concat([df1, df2], axis = 1)
        data.to_csv(path, index = False)
