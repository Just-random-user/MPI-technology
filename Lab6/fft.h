#pragma once
#include <vector>
#include <complex>
// v2: use fft + swap order to implement the ifft
namespace fft
{
	
	std::vector<std::complex<double>> fft(const std::vector<std::complex<double>> &vec);

	
	std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>> &vec);

	std::vector<int> round(const std::vector<std::complex<double>> &vec);

	
	std::vector<double> real(const std::vector<std::complex<double>> &vec);
} //namespace fft