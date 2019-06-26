/*
 * energy_detection.cpp
 *
 *  Created on: Jun 24, 2019
 *      Author: nkaramolegos
 */

#include "energy_detection.h"
#include <orbitsense/log.h>
#include <orbitsense/api.h>

namespace gr
{
  namespace orbitsense
  {

    energy_detection::energy_detection (const size_t fft_size,
                                        float energy_thresh_dB, bool nf_est,
                                        float noise_floor_val,
                                        float noise_floor_time,
                                        const double sampling_rate,
                                        uint8_t window) :
            d_fft_size (fft_size),
            d_energy_thresh_dB (energy_thresh_dB),
            d_nf_est (nf_est),
            d_noise_floor_val (noise_floor_val),
            d_noise_floor_time (noise_floor_time),
            d_sampling_rate (sampling_rate),
            d_window (window),
            d_noise_floor_est_cnt (
                floor (d_noise_floor_time / ((d_fft_size / d_sampling_rate))))
    {
      ORBITSENSE_DEBUG("Energy detection Method\n");

      d_fft = new fft::fft_complex (d_fft_size, true, 1);

      d_data_combined = (float*) volk_malloc (d_fft_size * sizeof(float), 32);

      d_shift_tmp = (float*) volk_malloc ((d_fft_size) * sizeof(float), 32);

      d_fftshift = (gr_complex*) volk_malloc ((d_fft_size) * sizeof(gr_complex),
                                              32);

      d_window_values = (float*) volk_malloc (d_fft_size * sizeof(float), 32);

      d_noise_floor_vec_dB = (float*) volk_malloc ((d_fft_size) * sizeof(float),
                                                   32);

      d_psd = (float*) volk_malloc ((d_fft_size) * sizeof(float), 32);

      d_energy_thresh_vec_dB = (float*) volk_malloc (
          (d_fft_size) * sizeof(float), 32);

      memset (d_energy_thresh_vec_dB, d_energy_thresh_dB,
              d_fft_size * sizeof(float));

      // create the value for each window
      create_window_values ();

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

    energy_detection::~energy_detection ()
    {
      delete d_fft;
      volk_free (d_window_values);
      volk_free (d_data_combined);
      volk_free (d_fftshift);
      volk_free (d_shift_tmp);
      volk_free (d_noise_floor_vec_dB);
      volk_free (d_psd);
    }

    void
    energy_detection::create_window_values ()
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
        }
    }

    void
    energy_detection::psd_estimation (const gr_complex *in)
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
    energy_detection::energy_detector (const gr_complex *in)
    {
      psd_estimation (in + d_rep_cnt * d_fft_size);

      volk_32f_x2_subtract_32f (d_psd, d_noise_floor_vec_dB, d_psd, d_fft_size);

      /* Get the signs */
      for (int i = 0; i < d_fft_size; i++) {

        d_psd[i] = d_psd[i] < 0;
      }

    }

    void
    energy_detection::noise_floor_estimation (const gr_complex *in)
    {

      /* estimate PSD */
      psd_estimation (in);

      /* Get maximum for each sub-carrier*/
      volk_32f_x2_max_32f (d_noise_floor_vec_dB, d_noise_floor_vec_dB, d_psd,
                           d_fft_size);
    }

    void
    energy_detection::energy_detection_init (const gr_complex *in,
                                             int noutput_items)
    {

      d_rep_cnt = 0;
      d_num_full_packets = noutput_items / d_fft_size;

      switch (d_mode)
        {
        case NOISE_FLOOR_ESTIMATION:

          while (d_rep_cnt < d_num_full_packets) {

            if (d_noise_floor_est_cnt > 0) {

              noise_floor_estimation (in);

              d_noise_floor_est_cnt--;
            }
            else {
              // add to NF the threshhold
              volk_32f_x2_add_32f (d_noise_floor_vec_dB, d_energy_thresh_vec_dB,
                                   d_noise_floor_vec_dB, d_fft_size);
              ORBITSENSE_DEBUG("Noise Floor estimation is done!\n");
              d_mode = SPECTRUM_DETECTION;
              break;
            }
            d_rep_cnt++;
          }
          break;
        case SPECTRUM_DETECTION:
          // check the selected method. You can split it in two classes.

          while (d_rep_cnt < d_num_full_packets) {

            energy_detector (in);
            //message_out_print (d_psd, d_fft_size);
            d_rep_cnt++;
          }

          break;
        }
    }

  } /* namespace orbitsense */
} /* namespace gr */
