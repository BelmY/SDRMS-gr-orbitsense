/*
 * energy_detection.h
 *
 *  Created on: Jun 24, 2019
 *      Author: nkaramolegos
 */

#ifndef LIB_ENERGY_DETECTION_H_
#define LIB_ENERGY_DETECTION_H_

#include <volk/volk.h>
#include <iostream>
#include <string.h>
#include <gnuradio/fft/fft.h>

namespace gr
{
  namespace orbitsense
  {

    class energy_detection
    {
    private:
      /* The FFT size */
      const size_t d_fft_size;

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

      size_t d_noise_floor_est_cnt;

      /* number of loops for all the samples */
      size_t d_num_full_packets;

      size_t d_rep_cnt;

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

      /* Threshold vector */
      float *d_energy_thresh_vec_dB;

    public:

      /* PSD variable */
      float *d_psd;

      energy_detection (const size_t fft_size, float energy_thresh_dB,
                        bool nf_est, float noise_floor_val,
                        float noise_floor_time, const double sampling_rate,
                        uint8_t window);

      ~energy_detection ();

      void
      create_window_values ();

      void
      psd_estimation (const gr_complex *in);

      void
      noise_floor_estimation (const gr_complex *in);

      void
      energy_detector (const gr_complex *in);

      void
      energy_detection_init (const gr_complex *in, int noutput_items);

    };

  } /* namespace orbitsense */
} /* namespace gr */

#endif /* LIB_ENERGY_DETECTION_H_ */
