#include "covariance_absolute_values.h"
#include <boost/math/special_functions/erf.hpp>
#include <volk/volk.h>
#include <cmath> 
#include <iostream>

namespace cav 
{
  CAV::CAV (const size_t num_samples, const uint8_t smoothing_factor,
            const double false_alarm_probability, const uint8_t est_snr) :
    d_num_samples (num_samples),
    d_smoothing_factor (smoothing_factor),
    d_false_alarm_probability (false_alarm_probability),
    d_est_snr (est_snr),
    d_t1 (0.0),
    d_t2 (0.0),
    d_threshold (0.0),
    d_samples_index (0),
    d_estimated_noise_power (0.0),
    d_estimated_noise_times (0),
    d_estimated_signal_noise_power (0.0)
  {
    d_samples = new std::complex<float>* [1];
    d_samples[0] = new std::complex<float> [num_samples + smoothing_factor];
    d_samples_unfilt = new std::complex<float> [d_num_samples];
    d_covariance_matrix = new std::complex<float> [smoothing_factor]; 
    compute_threshold();
  }

  CAV::~CAV ()
  {
    delete[] d_samples[0];
    delete[] d_samples;
    delete[] d_covariance_matrix;
  }

  void
  CAV::compute_covariance_matrix (const std::complex<float> *samples)
  {
    for (uint8_t i = 0; i < d_smoothing_factor; i++) {
      d_covariance_matrix[i] = compute_autocorrelations(samples, i);
    }
  }

  std::complex<float>
  CAV::compute_autocorrelations (const std::complex<float> *samples,
                                 uint8_t lamda)
  {
    std::complex<float> sum = std::complex<float> (0, 0);
    volk_32fc_x2_conjugate_dot_prod_32fc (&sum, &samples[d_smoothing_factor],
                                          &samples[d_smoothing_factor - lamda],
                                          d_num_samples);

    return std::complex<float> (sum.real() / d_num_samples,
                                sum.imag() / d_num_samples);
  }

  void
  CAV::compute_correlations (void)
  {
    d_t1 = 0.0;
    d_t2 = 0.0; 

    for (uint8_t i = 1; i < d_smoothing_factor; i++) {
      d_t1 += (d_smoothing_factor - i) * std::abs(d_covariance_matrix[i]);
    }
    d_t1 = std::abs(d_covariance_matrix[0]) + (2.0 * d_t1) / d_smoothing_factor;
    d_t2 = std::abs(d_covariance_matrix[0]);
  }

  void
  CAV::compute_threshold (void)
  {
    
    const float pfa = d_false_alarm_probability;
    /* Q^-1(x) = sqrt(2) * inverse complementary error function(2x) */
    d_threshold = (1.0 +
              ((d_smoothing_factor - 1.0) * std::sqrt (2.0 / d_num_samples * pi)))  
              / (1.0 -
                (std::sqrt (2.0) * boost::math::erfc_inv (2.0 * pfa)
                  * std::sqrt (2.0 / d_num_samples)));
  }

  std::vector<int>
  CAV::unfiltered_save_run (const std::complex<float> *samples, size_t nitems)
  {
    /* If less samples that needed */
    if (d_samples_index + nitems < d_num_samples) {
      memcpy(&d_samples[0][d_smoothing_factor + d_samples_index], samples,
             sizeof(std::complex<float>) * nitems);
      d_samples_index += nitems;
      return std::vector<int> (1, 2);
    }
    /* If enough samples */
    else if (d_samples_index + nitems == d_num_samples) {
      memcpy(&d_samples[0][d_smoothing_factor + d_samples_index], samples,
             sizeof(std::complex<float>) * nitems);
      
      compute_covariance_matrix (d_samples[0]);
      compute_correlations ();
      
      /* Copy smoothing factor samples for next run */
      memcpy(d_samples[0], &d_samples[d_num_samples - d_smoothing_factor],
             sizeof(std::complex<float>) * d_smoothing_factor);

      d_samples_index = 0;
    }
    /* If enough samples plus extra */
    else {
      int extra_samples = d_samples_index + nitems - d_num_samples;
      /* Save needed samples to fill buffer */
      memcpy(&d_samples[0][d_smoothing_factor + nitems - extra_samples],
             samples, sizeof(std::complex<float>) * (nitems - extra_samples));
      
      compute_covariance_matrix (d_samples[0]);
      compute_correlations ();

      /* Copy smoothing factor samples for next run */
      memcpy(d_samples[0], &d_samples[d_num_samples - d_smoothing_factor],
             sizeof(std::complex<float>) * d_smoothing_factor);

      /* Copy extra samples for next run */
      memcpy(&d_samples[0][d_smoothing_factor], &samples[nitems - extra_samples],
             sizeof(std::complex<float>) * extra_samples);
    }
    
    if (d_t1 / d_t2 > d_threshold) {
      if (d_est_snr && d_estimated_noise_times < 50) {
        std::cout << "Noise power estimation not ready" << std::endl;
      }
      if (d_est_snr && d_estimated_noise_times >= 50) {
        print_snr_db ();
      }
      return std::vector<int> (1, 1);
    }
    else {
      estimate_noise_rms (&d_samples[0][d_smoothing_factor], d_num_samples);
      return std::vector<int> (1, 0);
    }
  }

  void
  CAV::estimate_noise_rms (std::complex<float>* in, size_t nitems)
  {
    double magnitude = 0.0;
    double alpha = 0.0001;
    double beta = 1 - alpha;
    
    for (int i = 0; i < nitems; i++) {
      magnitude = std::pow(std::abs(in[i]), 2);
      d_estimated_noise_power = beta * d_estimated_noise_power + alpha * magnitude;
    }
    d_estimated_noise_times = std::min(++d_estimated_noise_times, 50u);
    if (d_estimated_noise_times == 50u && d_est_snr) {
      std::cout << "Noise power estimation ready!" << std::endl;
    }
  }

  void
  CAV::estimate_signal_noise_rms (std::complex<float>* in, size_t nitems)
  {
    double magnitude = 0.0;
    double alpha = 0.0001;
    double beta = 1 - alpha;
    for (int i = 0; i < nitems; i++) {
      magnitude = std::pow(std::abs(in[i]), 2);
      d_estimated_signal_noise_power = beta * d_estimated_signal_noise_power + 
                                        alpha * magnitude;
    }
  }

  void
  CAV::print_snr_db (void)
  {
    estimate_signal_noise_rms (&d_samples[0][d_smoothing_factor], d_num_samples);
    std::cout << "SNR: " << 10*log10(
                    (d_estimated_signal_noise_power - d_estimated_noise_power)/
                     d_estimated_noise_power) << " dB" << std::endl;
  }
  
  std::vector<int> 
  CAV::covariance_absolute_values_engine (const std::complex<float> *samples,
                                          size_t nitems)
  {
    return unfiltered_save_run (samples, nitems);
  }

} /* namespace cav */
