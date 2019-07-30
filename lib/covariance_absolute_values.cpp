#include "covariance_absolute_values.h"
#include <boost/math/special_functions/erf.hpp>
#include <volk/volk.h>
#include <cmath> 
#include <iostream>

namespace cav 
{
  CAV::CAV (const size_t num_samples, const uint8_t smoothing_factor,
            const double false_alarm_probability) :
    d_num_samples (num_samples),
    d_smoothing_factor (smoothing_factor),
    d_false_alarm_probability (false_alarm_probability),
    d_t1 (0.0),
    d_t2 (0.0),
    d_threshold (0.0)
  {
    d_samples = new std::complex<float>[num_samples + smoothing_factor];
    d_covariance_matrix = new std::complex<float>[smoothing_factor]; 
    compute_threshold();
  }

  CAV::~CAV ()
  {
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

  int
  CAV::covariance_absolute_values_engine (const std::complex<float> *samples)
  {
    compute_covariance_matrix (samples);
    compute_correlations();

    if (d_t1 / d_t2 > d_threshold) {
      return 1;
    }
    else {
      return 0;
    }
  }

} /* namespace cav */
