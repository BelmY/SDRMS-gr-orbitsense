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
#include <math.h>
#include <boost/fiber/condition_variable.hpp>
#include <boost/fiber/all.hpp>
#include <volk/volk.h>
#include <orbitsense/log.h>
#include <energy_detection.h>

namespace gr
{
  namespace orbitsense
  {

    class detection_engine_impl : public detection_engine
    {
    private:

      /* Define static pi */
      const double pi = 3.141592653589793;

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

      /* Duty cycle class instance */
      energy_detection* d_energy_detection;

      /* Covariance matrix for detection */
      gr_complex **d_covariance_matrix;
	
      /* Smoothing factor previous samples array */
      gr_complex *d_prev_samples;

      /* Number of samples used for covariance based detection */
      const size_t d_num_samples;

      /* Smoothing factor for covariance based detection */
      uint8_t d_smoothing_factor;
	
      /* Probability of false alarm for covariance based detection */
      float d_false_alarm_probability;
	
      /* Threshold based on false alarm probability */
      double d_threshold;

    public:
      detection_engine_impl (const size_t fft_size, uint8_t method,
                             float energy_thresh_dB, uint8_t nf_est,
                             float noise_floor_val, float noise_floor_time,
                             const double sampling_rate, uint8_t window,
			     const size_t num_samples, uint8_t smoothing_factor,
			     float false_alarm);
      ~detection_engine_impl ();

      // Where all the action really happens
      int
      work (int noutput_items, gr_vector_const_void_star &input_items,
            gr_vector_void_star &output_items);
      
      void
      compute_covariance_matrix(const gr_complex *in);

      gr_complex
      compute_autocorrelations(const gr_complex *in, size_t lamda);

      void
      compute_correlations(gr_complex **matrix, uint8_t smoothing_factor,
			 double *thres1, double *thres2);

      void
      compute_threshold(const float probability_false_alarm, size_t num_samples,
		        uint8_t smoothing_factor, double *threshold);
    
      message_out_print (float *vector, int vector_len);

    };

  } // namespace orbitsense
} // namespace gr

#endif /* INCLUDED_ORBITSENSE_DETECTION_ENGINE_IMPL_H */

