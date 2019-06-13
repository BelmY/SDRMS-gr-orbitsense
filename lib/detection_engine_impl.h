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

#ifndef INCLUDED_ORBITSENSE_DETECTION_ENGINE_IMPL_H
#define INCLUDED_ORBITSENSE_DETECTION_ENGINE_IMPL_H

#include <orbitsense/detection_engine.h>
#include <gnuradio/fft/fft.h>
#include <math.h>
#include <boost/fiber/condition_variable.hpp>
#include <boost/fiber/all.hpp>
#include <volk/volk.h>
#include <orbitsense/log.h>

namespace gr
{
  namespace orbitsense
  {

    class detection_engine_impl : public detection_engine
    {
    private:

      /* The FFT size */
      const size_t d_fft_size;

      /* Method of sensing engine */
      uint8_t d_method;

      /* Threshold for checking if signal exists */
      float d_energy_thresh_dB;

      /* flag for calculating the noise floor */
      bool d_nf_est;

      /* manual noise floor value in dB */
      float d_noise_floor_val;

      /* sec for estimating noise floor */
      float d_noise_floor_time;

      /* Operational sampling rate */
      const double d_sampling_rate;

      uint8_t d_window;

      uint8_t d_mode;

      /* Auxiliary FFT vector */
      fft::fft_complex *d_fft;

      /* Auxiliary vector */
      float *d_data_combined;

      /* Auxiliary vector to store FFT shift */
      float *d_shift_tmp;

      /* Auxiliary vector to store FFT shift for each input */
      gr_complex *d_fftshift;

      /* window */
      float *d_window_values;

      /* Noise floor values */
      float *d_noise_floor_vec_dB;

      /* PSD variable */
      float *d_psd;

      /* Threshold vector */
      float *d_energy_thresh_vec_dB;

      /* number of loops for all the samples */
      size_t d_num_full_packets;

      size_t d_rep_cnt;

      size_t d_noise_floor_est_cnt;

    public:
      detection_engine_impl (const size_t fft_size, uint8_t method,
                             float energy_thresh_dB, uint8_t nf_est,
                             float noise_floor_val, float noise_floor_time,
                             const double sampling_rate, uint8_t window);
      ~detection_engine_impl ();

      // Where all the action really happens
      int
      work (int noutput_items, gr_vector_const_void_star &input_items,
            gr_vector_void_star &output_items);

      void
      create_window_values ();

      void
      psd_estimation (const gr_complex *in);

      void
      noise_floor_estimation (const gr_complex *in);

      void
      energy_detector(const gr_complex *in);

      void
      message_out_print(float *vector, int vector_len);
    };

  } // namespace orbitsense
} // namespace gr

#endif /* INCLUDED_ORBITSENSE_DETECTION_ENGINE_IMPL_H */

