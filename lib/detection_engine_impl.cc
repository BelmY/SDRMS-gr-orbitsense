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
#include "detection_engine_impl.h"

namespace gr
{
  namespace orbitsense
  {

    detection_engine::sptr
    detection_engine::make (const size_t fft_size, uint8_t method,
                            float energy_thresh_dB, uint8_t nf_est,
                            float noise_floor_val, float noise_floor_time,
                            const double sampling_rate, uint8_t window)
    {
      return gnuradio::get_initial_sptr (
          new detection_engine_impl (fft_size, method, energy_thresh_dB, nf_est,
                                     noise_floor_val, noise_floor_time,
                                     sampling_rate, window));
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
                                                  uint8_t window) :
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
            d_num_full_packets (0),
            d_rep_cnt (0),
            d_noise_floor_est_cnt (
                floor (d_noise_floor_time / ((d_fft_size / d_sampling_rate))))
    {
      message_port_register_out (pmt::mp ("data_out"));

      /* Process in a per-FFT basis */
      set_output_multiple (d_fft_size);

      d_fft = new fft::fft_complex (d_fft_size, true, 1);

      d_data_combined = (float*) volk_malloc (d_fft_size * sizeof(float), 32);

      d_shift_tmp = (float*) volk_malloc ((d_fft_size) * sizeof(float), 32);

      d_fftshift = (gr_complex*) volk_malloc ((d_fft_size) * sizeof(gr_complex),
                                              32);

      d_window_values = (float*) volk_malloc (d_fft_size * sizeof(float), 32);

      d_noise_floor_vec_dB = (float*) volk_malloc ((d_fft_size) * sizeof(float),
                                                   32);

      d_psd = (float*) volk_malloc ((d_fft_size) * sizeof(float), 32);

      d_energy_thresh_vec_dB = (float*) volk_malloc ((d_fft_size) * sizeof(float), 32);

      memset (d_energy_thresh_vec_dB, d_energy_thresh_dB, d_fft_size * sizeof(float));

      // create the value for each window
      create_window_values ();

      std::cout << d_nf_est << std::endl;

      if (d_nf_est == true) {
        d_mode = NOISE_FLOOR_ESTIMATION;
        ORBITSENSE_DEBUG("NF will be estimated\n");
      }
      else {
        // set manual NF
        for (int i = 0; i < d_fft_size; i++) {
          d_noise_floor_vec_dB[i] = d_noise_floor_val;
        }
        ORBITSENSE_DEBUG(
            "There is no NF estimation, setting manually at %2f db",
            d_noise_floor_val);
        d_mode = SPECTRUM_DETECTION;
      }

    }

    /*
     * Our virtual destructor.
     */
    detection_engine_impl::~detection_engine_impl ()
    {
      delete d_fft;
      volk_free (d_window_values);
      volk_free (d_data_combined);
      volk_free (d_fftshift);
      volk_free (d_shift_tmp);
      volk_free (d_noise_floor_vec_dB);
      volk_free (d_psd);
    }

    int
    detection_engine_impl::work (int noutput_items,
                                 gr_vector_const_void_star &input_items,
                                 gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];

      // for now we don't check if there is signal
      d_num_full_packets = noutput_items / d_fft_size;
      d_rep_cnt = 0;

      switch (d_mode)
        {
        case NOISE_FLOOR_ESTIMATION:
          noise_floor_estimation (in);
          for (int i = 0; i < d_fft_size; i++) {
            ORBITSENSE_DEBUG("%2f, ", d_noise_floor_vec_dB[i]);
          }
          // add to NF the threshhold
          volk_32f_x2_add_32f (d_noise_floor_vec_dB, d_energy_thresh_vec_dB, d_noise_floor_vec_dB, d_fft_size);
          break;
        case SPECTRUM_DETECTION:
          // check the selected method. You can split it in two classes.
          switch (d_method)
            {
            case ENERGY_DETECTION:

              while (d_rep_cnt < d_num_full_packets) {

                energy_detector (in);
                message_out_print (d_psd, d_fft_size);
                d_rep_cnt++;
              }

              break;

            case COVARIANCE:

              break;
            }

          break;
        }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

    void
    detection_engine_impl::create_window_values ()
    {
      switch (d_window)
        {
        case FLAT_TOP:
          for (size_t i = 0; i < d_fft_size; i++) {
            d_window_values[i] = 0.215578950
                - 0.416631580 * cos ((2 * M_PI * i) / (d_fft_size - 1))
                + 0.2772631580 * cos ((4 * M_PI * i) / (d_fft_size - 1))
                - 0.0835789470 * cos ((6 * M_PI * i) / (d_fft_size - 1))
                + 0.0069473680 * cos ((8 * M_PI * i) / (d_fft_size - 1));
            // set NF to something very small
            d_noise_floor_vec_dB[i] = -400;
          }
          break;

        case BLACKMANN_HARRIS:
          /* Calculation of Blackmann-Harris window */
          for (size_t i = 0; i < d_fft_size; i++) {
            d_window_values[i] = 0.35875
                - 0.48829 * cos ((2 * M_PI * i) / (d_fft_size - 1))
                + 0.14128 * cos ((4 * M_PI * i) / (d_fft_size - 1))
                - 0.01168 * cos ((6 * M_PI * i) / (d_fft_size - 1));
            // set NF to something very small
            d_noise_floor_vec_dB[i] = -400;
          }
          break;
        case NONE:
          for (size_t i = 0; i < d_fft_size; i++) {
            d_window_values[i] = 1;
            // set NF to something very small
            d_noise_floor_vec_dB[i] = -400;
          }
          break;
        }
    }

    void
    detection_engine_impl::psd_estimation (const gr_complex *in)
    {
      /* Apply window */
      volk_32fc_32f_multiply_32fc (&d_fft->get_inbuf ()[0], in, d_window_values,
                                   d_fft_size);

      /* Perform fft */
      d_fft->execute ();
      /* Perform fftshift  */
      memcpy (&d_fftshift[0], &d_fft->get_outbuf ()[(d_fft_size / 2)],
              sizeof(gr_complex) * (d_fft_size / 2));
      memcpy (&d_fftshift[(d_fft_size / 2)], &d_fft->get_outbuf ()[0],
              sizeof(gr_complex) * (d_fft_size / 2));

      /* Calculate psd in dB. We can change it to linear if it is faster*/
      volk_32fc_s32f_x2_power_spectral_density_32f (d_psd, d_fftshift,
                                                    d_fft_size, 1, d_fft_size);
    }

    void
    detection_engine_impl::noise_floor_estimation (const gr_complex *in)
    {

      d_rep_cnt = 0;
      while (d_rep_cnt < d_num_full_packets) {

        if (d_noise_floor_est_cnt > 0) {

          /* estimate PSD */
          psd_estimation (in);

          /* Get maximum for each sub-carrier*/
          volk_32f_x2_max_32f (d_noise_floor_vec_dB, d_noise_floor_vec_dB,
                               d_psd, d_fft_size);
          d_noise_floor_est_cnt--;
        }
        else {
          d_mode = SPECTRUM_DETECTION;
          break;
        }
        d_rep_cnt++;
      }
    }

    void
    detection_engine_impl::energy_detector (const gr_complex *in)
    {
      psd_estimation (in + d_rep_cnt * d_fft_size);

      volk_32f_x2_subtract_32f (d_psd, d_noise_floor_vec_dB, d_psd, d_fft_size);

      /* Get the signs */
      for (int i = 0; i < d_fft_size; i++) {

        d_psd[i] = d_psd[i] < 0;
      }

    }

    void
    detection_engine_impl::message_out_print (float *vector, int vector_len)
    {
      pmt::pmt_t vector_combined = pmt::init_f32vector (vector_len, vector);
      message_port_pub (pmt::mp ("data_out"), vector_combined);
    }

  } /* namespace orbitsense */
} /* namespace gr */

