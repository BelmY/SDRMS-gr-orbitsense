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

#ifndef INCLUDED_ORBITSENSE_DETECTION_ENGINE_H
#define INCLUDED_ORBITSENSE_DETECTION_ENGINE_H

#include <orbitsense/api.h>
#include <gnuradio/sync_block.h>

namespace gr
{
  namespace orbitsense
  {

    /*!
     * \brief <+description of block+>
     * \ingroup orbitsense
     *
     */
    class ORBITSENSE_API detection_engine : virtual public gr::sync_block
    {
    public:
      typedef boost::shared_ptr<detection_engine> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of orbitsense::detection_engine.
       *
       * To avoid accidental use of raw pointers, orbitsense::detection_engine's
       * constructor is in a private implementation
       * class. orbitsense::detection_engine::make is the public interface for
       * creating new instances.
       */
      static sptr
      make (const size_t fft_size, uint8_t method, float energy_thresh_dB,
            uint8_t nf_est, float noise_floor_val, float d_noise_floor_time, const double sampling_rate,
            uint8_t window, const size_t num_samples, uint8_t smoothing_factor, float false_alarm);
    };

  } // namespace orbitsense
} // namespace gr

#endif /* INCLUDED_ORBITSENSE_DETECTION_ENGINE_H */

