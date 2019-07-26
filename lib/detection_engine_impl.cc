/* -*- c++ -*- */
/* 
 * Copyright 2019 gr-orbitsense author.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <boost/math/special_functions/erf.hpp>
#include "detection_engine_impl.h"
#include <thread>
#include <chrono>

namespace gr
{
  namespace orbitsense
  {

    detection_engine::sptr
    detection_engine::make (const size_t fft_size, uint8_t method,
                            float energy_thresh_dB, uint8_t nf_est,
                            float noise_floor_val, float noise_floor_time,
                            const double sampling_rate, uint8_t window,
                            const size_t num_samples, uint8_t smoothing_factor,
                            float false_alarm)
    {
      return gnuradio::get_initial_sptr (
          new detection_engine_impl (fft_size, method, energy_thresh_dB, nf_est,
                                     noise_floor_val, noise_floor_time,
                                     sampling_rate, window, num_samples,
                                     smoothing_factor, false_alarm));
    }

    /*
     * The private constructor
     */
    detection_engine_impl::detection_engine_impl (const size_t fft_size,
                                                  uint8_t method,
                                                  float energy_thresh_dB,
                                                  uint8_t nf_est,
                                                  float noise_floor_val,
                                                  float noise_floor_time,
                                                  const double sampling_rate,
                                                  uint8_t window,
                                                  const size_t num_samples,
                                                  uint8_t smoothing_factor,
                                                  float false_alarm) :
            gr::sync_block ("detection_engine",
                            gr::io_signature::make (1, 1, sizeof(gr_complex)),
                            gr::io_signature::make (0, 0, 0)),
            d_fft_size (fft_size),
            d_method (method),
            d_energy_thresh_dB (energy_thresh_dB),
            d_nf_est (nf_est),
            d_noise_floor_val (noise_floor_val),
            d_noise_floor_time (noise_floor_time),
            d_sampling_rate (sampling_rate),
            d_window (window),
            d_num_samples (num_samples),
            d_smoothing_factor (smoothing_factor),
            d_false_alarm_probability (false_alarm)
    {
      message_port_register_out (pmt::mp ("data_out"));
      /* Process in a per-FFT basis */
      set_output_multiple (d_fft_size);

      switch (d_method)
        {
        case ENERGY_DETECTION:
          d_energy_detection = new energy_detection (d_fft_size,
                                                     d_energy_thresh_dB,
                                                     d_nf_est,
                                                     d_noise_floor_val,
                                                     d_noise_floor_time,
                                                     d_sampling_rate, d_window);
          break;
        case COVARIANCE:
          set_output_multiple (d_num_samples);
          d_prev_samples = new gr_complex[d_smoothing_factor]();
          /* Initialize CAV instance */
          d_cav_engine = new cav::CAV(d_num_samples,
                                    d_smoothing_factor,
                                    d_false_alarm_probability);
          break;
      }
    }

    /*
     * Our virtual destructor.
     */
    detection_engine_impl::~detection_engine_impl ()
    {
      delete[] d_prev_samples;
      delete d_energy_detection;
      delete d_cav_engine;
    }

    int
    detection_engine_impl::work (int noutput_items,
                                 gr_vector_const_void_star &input_items,
                                 gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];

      switch (d_method)
        {
        case ENERGY_DETECTION:
          d_energy_detection->energy_detection_init (in, noutput_items);
          message_out_print (d_energy_detection->d_psd, d_fft_size);
          break;
        case COVARIANCE:
          /* Save number of samples + smoothing factor for processing */
          gr_complex tmp_input[d_num_samples + d_smoothing_factor];
          size_t last_index;
          size_t d_rep_cnt = 0;
          size_t d_num_full_packets = noutput_items / d_num_samples;
          int detected = 0;

          std::memcpy (&tmp_input[0], d_prev_samples,
                       d_smoothing_factor * sizeof(gr_complex));
          while (d_rep_cnt < d_num_full_packets) {
            std::memcpy (&tmp_input[d_smoothing_factor],
                         &in[d_rep_cnt * d_num_samples],
                         d_num_samples * sizeof(gr_complex));

            /* Run CAV on new samples */
            detected = 
                  d_cav_engine->covariance_absolute_values_engine (tmp_input);
            /*
            if (detected == 1) {
              ORBITSENSE_DEBUG("Signal_detected!");
            }*/

            message_out_print(detected);
            
            /* Save last smoothing factor samples for next run*/
            last_index = (d_rep_cnt + 1) * d_num_samples
                         - d_smoothing_factor;
            std::memcpy (&tmp_input[0], &in[last_index],
                         d_smoothing_factor * sizeof(gr_complex));

            d_rep_cnt++;
          }
          /* Save previous samples for next run */
          std::memcpy (d_prev_samples, &in[last_index],
                       d_smoothing_factor * sizeof(gr_complex));
          break;
        }
      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

    void
    detection_engine_impl::message_out_print (float *vector, int vector_len)
    {
      pmt::pmt_t vector_combined = pmt::init_f32vector (vector_len, vector);
      message_port_pub (pmt::mp ("data_out"), vector_combined);
    }

    void
    detection_engine_impl::message_out_print (bool in)
    {
      pmt::pmt_t detected = pmt::from_bool(in);
      message_port_pub (pmt::mp ("data_out"), detected);
    }

  } /* namespace orbitsense */
} /* namespace gr */
