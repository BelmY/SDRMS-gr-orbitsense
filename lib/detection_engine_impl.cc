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
            d_window (window)
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
          break;
        }

    }

    /*
     * Our virtual destructor.
     */
    detection_engine_impl::~detection_engine_impl ()
    {
      delete d_energy_detection;
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
          //twra tipwnei kathe noutput_items des ti tha kaneis apo edw kai pera
          message_out_print (d_energy_detection->d_psd, d_fft_size);
          break;
        case COVARIANCE:
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

  } /* namespace orbitsense */
} /* namespace gr */

