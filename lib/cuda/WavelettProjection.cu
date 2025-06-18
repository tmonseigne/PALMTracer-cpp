#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <vector>

// extern "C" Fonction appelable depuis un CPP
// __device__ Fonction CUDA
// __global__ Kernel CUDA

constexpr double COEFF_H0 = 3.0 / 8.0, COEFF_H1 = 1.0 / 4.0, COEFF_H2 = 1.0 / 16.0;


__global__ void convertToDouble(const uint16_t* in, double* out, const size_t size)
{
	const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= size) { return; }
	out[idx] = static_cast<double>(in[idx]);
}

__global__ void WaveletHorizontal(const double* in, double* tmp, const size_t width, const size_t height, const size_t plane)
{
	const size_t X = blockIdx.x * blockDim.x + threadIdx.x;
	const size_t Y = blockIdx.y * blockDim.y + threadIdx.y;
	if (X >= width || Y >= height) { return; }
	const int Step       = 1 << plane; // 2^plane
	const int DoubleStep = Step << 1;  // 2*2^plane

	const auto Idx = Y * width + X;

	// Lecture sécurisée avec bord clamping
	auto get = [&](const size_t xx) -> double {
		const size_t Xi = max(static_cast<size_t>(0), min(width - 1, xx));
		return in[Y * width + Xi];
	};

	tmp[Idx] = COEFF_H0 * in[Idx] + COEFF_H1 * (get(X - Step) + get(X + Step)) + COEFF_H2 * (get(X - DoubleStep) + get(X + DoubleStep));
}

__global__ void WaveletVertical(const double* tmp, double* out, const size_t width, const size_t height, const size_t plane)
{
	const size_t X = blockIdx.x * blockDim.x + threadIdx.x;
	const size_t Y = blockIdx.y * blockDim.y + threadIdx.y;
	if (X >= width || Y >= height) { return; }

	const int Step       = 1 << plane;
	const int DoubleStep = Step << 1;

	const auto Idx = Y * width + X;

	// Lecture sécurisée avec bord clamping
	auto get = [&](const size_t yy) -> double {
		const size_t Yi = max(static_cast<size_t>(0), min(height - 1, yy));
		return tmp[Yi * width + X];
	};

	out[Idx] = COEFF_H0 * tmp[Idx] + COEFF_H1 * (get(Y - Step) + get(Y + Step)) + COEFF_H2 * (get(Y - DoubleStep) + get(Y + DoubleStep));
}


__global__ void Wavelet_substract(const double* in, double* out, const size_t size)
{
	const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= size) { return; }
	out[idx] -= in[idx];
}

//Parsing
extern "C"
void WaveletProjections(const uint16_t* image, const size_t height, const size_t width, const size_t nbPlane, std::vector<std::vector<double>>& decomp)
{
	//copie de l'image dans decomp[0] en passant en double
	//std::transform(image, image + size_, wDecomp_[0].begin(), [](const uint16_t val) { return static_cast<double>(val); });
	const size_t Size      = height * width;
	const size_t TotalSize = Size * nbPlane;

	// Allocation GPU
	uint16_t* image_d;
	double *decomp_d, *tmp_d;
	cudaMalloc(&image_d, Size * sizeof(uint16_t));
	cudaMalloc(&decomp_d, TotalSize * sizeof(double));
	cudaMalloc(&tmp_d, Size * sizeof(double));

	// Transfert de l'image
	cudaMemcpy(image_d, image, Size * sizeof(uint16_t), cudaMemcpyHostToDevice);

	// Appel du kernel principal CUDA répartition en 256 threads et force à passer sur totue l'image si l'image d'est pas une puissance de 2.
	dim3 blockSize(16, 16);
	dim3 gridSize((width + 15) / 16, (height + 15) / 16);
	convertToDouble<<<gridSize, blockSize>>>(image_d, decomp_d, Size);
	cudaDeviceSynchronize(); // Attendre que decomp_d soit rempli

	for (size_t i = 0; i < nbPlane - 1; ++i) {
		WaveletHorizontal<<<gridSize, blockSize>>>(decomp_d + i * Size, tmp_d, width, height, i);
		cudaDeviceSynchronize(); // Attendre que tmp soit rempli
		WaveletVertical<<<gridSize, blockSize>>>(tmp_d, decomp_d + (i + 1) * Size, width, height, i);
		cudaDeviceSynchronize(); // Attendre que out soit prêt avant de le modifier
		Wavelet_substract<<<gridSize, blockSize>>>(decomp_d + (i + 1) * Size, decomp_d + i * Size, Size);
		cudaDeviceSynchronize(); // Attendre que la décomposition soit fini
	}

	cudaDeviceSynchronize(); // S’assurer que tout le GPU a fini AVANT transfert
	// Récupération des résultats
	for (size_t p = 0; p < nbPlane; ++p) { cudaMemcpy(decomp[p].data(), decomp_d + p * Size, Size * sizeof(double), cudaMemcpyDeviceToHost); }

	// Libération mémoire
	cudaFree(image_d);
	cudaFree(decomp_d);
	cudaFree(tmp_d);
}
