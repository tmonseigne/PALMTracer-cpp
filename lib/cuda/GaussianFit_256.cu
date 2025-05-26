#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cmath>

// extern "C" Fonction appelable depuis un CPP
// __device__ Fonction CUDA
// __global__ Kernel CUDA

constexpr size_t GRID_SIZE = 256;


__device__ double Psf(const double x, const double y, const double i0, const double sigma, const double offset)
{
	return i0 * exp(-(x * x + y * y) / (2 * (sigma * sigma))) + offset;
}

__device__ double Psf(const double x, const double y, const double i0, const double sigmaX, const double sigmaY, const double theta, const double offset)
{
	const auto CosTheta = cos(theta), SinTheta = sin(theta);
	const auto Factor1  = (x * CosTheta - y * SinTheta) / sigmaX,
	           Factor2  = (x * SinTheta + y * CosTheta) / sigmaY;

	return i0 * exp(-(Factor1 * Factor1 + Factor2 * Factor2)) + offset;
}

__device__ double Psf(const double* p, const double x, const double y, const int mode, const double sigma, const double theta)
{
	switch (mode) {
		case 0: return Psf(x - p[1], y - p[2], p[0], sigma, p[3]);
		case 1: return Psf(x - p[1], y - p[2], p[0], p[3], p[4]);
		case 2: return Psf(x - p[1], y - p[2], p[0], p[3], p[4], theta, p[5]);
		case 3: return Psf(x - p[1], y - p[2], p[0], p[3], p[4], p[5], p[6]);
		default: return 0;
	}
}


__global__ void ComputeHatDelta(const double* dat, const size_t* t, double* hat, double* delta, const double* p,
                                const int mode, const double sigma, const double theta, const size_t n)
{
	const auto Idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (Idx >= n) { return; }
	const double Val = Psf(p, static_cast<double>(t[Idx * 2]), static_cast<double>(t[Idx * 2 + 1]), mode, sigma, theta);
	hat[Idx]         = Val;
	delta[Idx]       = dat[Idx] - Val;
}

extern "C"
void LaunchComputeHatDelta(const double* dat, const size_t* t, double* hat, double* delta,
                           const double* p, const int mode, const double sigma, const double theta, const size_t n)
{
	double *d_dat, *d_hat, *d_delta, *d_p;
	size_t* d_t;

	cudaMalloc(&d_dat, n * sizeof(double));
	cudaMalloc(&d_hat, n * sizeof(double));
	cudaMalloc(&d_delta, n * sizeof(double));
	cudaMalloc(&d_t, 2 * n * sizeof(size_t));

	switch (mode) {
		case 0: cudaMalloc(&d_p, 4 * sizeof(double));
			cudaMemcpy(d_p, p, 4 * sizeof(double), cudaMemcpyHostToDevice);
			break;
		case 1: cudaMalloc(&d_p, 5 * sizeof(double));
			cudaMemcpy(d_p, p, 5 * sizeof(double), cudaMemcpyHostToDevice);
			break;
		case 2: cudaMalloc(&d_p, 6 * sizeof(double));
			cudaMemcpy(d_p, p, 6 * sizeof(double), cudaMemcpyHostToDevice);
			break;
		case 3: cudaMalloc(&d_p, 7 * sizeof(double));
			cudaMemcpy(d_p, p, 7 * sizeof(double), cudaMemcpyHostToDevice);
			break;
		default: break;
	}

	cudaMemcpy(d_dat, dat, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_t, t, 2 * n * sizeof(size_t), cudaMemcpyHostToDevice);

	dim3 blockSize(GRID_SIZE);
	dim3 gridSize((n + blockSize.x - 1) / blockSize.x);

	ComputeHatDelta<<< gridSize, blockSize >>>(d_dat, d_t, d_hat, d_delta, d_p, mode, sigma, theta, n);

	cudaMemcpy(hat, d_hat, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(delta, d_delta, n * sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(d_dat);
	cudaFree(d_hat);
	cudaFree(d_delta);
	cudaFree(d_p);
	cudaFree(d_t);
}
