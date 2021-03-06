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

#ifndef INCLUDED_COVARIANCE_ABSOLUTE_VALUES_H
#define INCLUDED_COVARIANCE_ABSOLUTE_VALUES_H

#include <complex>
#include <vector>

namespace cav
{
  class CAV
  {
  private:

    /* Define pi */
    const double pi = 3.141592653589793;

    /* Number of samples used for covariance based detection */
    const size_t d_num_samples;

    /* Smoothing factor for covariance based detection */
    uint8_t d_smoothing_factor;

    /* Probability of false alarm for covariance based detection */
    float d_false_alarm_probability;

    /* SNR estimation flag */
    uint8_t d_est_snr;

    /* Sum of main diagonal of covariance matrix */
    double d_t1;

    /* Sum of other elements of covariance matrix */
    double d_t2;

    /* Threshold based on false alarm probability */
    double d_threshold;

    /* Covariance matrix for detection */
    std::complex<float>* d_covariance_matrix;

    /* Number of samples plus smoothing factor samples array */
    std::complex<float>** d_samples;

    /* Number of saved samples */
    unsigned int d_samples_index;

    /* Unfiltered samples plus history samples array */
    std::complex<float>* d_samples_unfilt;

    /* Estimated noise power */
    double d_estimated_noise_power;

    /* Times of noise estimation */
    unsigned int d_estimated_noise_times;

    /* Estimated signal + noise power */
    double d_estimated_signal_noise_power;

    void
    compute_covariance_matrix (const std::complex<float> *samples);

    std::complex<float>
    compute_autocorrelations (const std::complex<float> *samples, 
                              uint8_t lamda);

    void
    compute_correlations (void);

    void
    compute_threshold (void);

    std::vector<int>
    unfiltered_save_run (const std::complex<float> *samples, size_t nitems);

    void
    estimate_noise_rms (std::complex<float> *in, size_t nitems);

    void
    estimate_signal_noise_rms (std::complex<float>* in, size_t nitems);

    void
    print_snr_db (void);

    void
    run_cav (int i);

  public:

    CAV (const size_t num_samples, const uint8_t smoothing_factor,
         const double false_alarm_probability, const uint8_t est_snr);

    ~CAV ();

    std::vector<int> 
    covariance_absolute_values_engine (const std::complex<float> *samples, size_t nitems); 
  };
} /* namespase cav */

#endif /* INCLUDED_COVARIANCE_ABSOLUTE_VALUES_H */ 
