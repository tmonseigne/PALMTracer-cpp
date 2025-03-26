#ifndef _DECOMPOSITION_KERNEL_H_
#define _DECOMPOSITION_KERNEL_H_

#include <helper_cuda.h>
#include <helper_timer.h>

#include "../parameters.h"

__global__ void GPUuSrtToFloat(unsigned short* input, float* output, unsigned short WIDTH, unsigned short HEIGHT);


extern "C"
void test(unsigned short* input, float* output, unsigned short numImages, unsigned short width, unsigned short height)
{
	dim3  grid( numImages * (unsigned int)ceil( (float)width / (float)WIDTH_PORTION ), (unsigned int)ceil( (float)height / (float)HEIGHT_PORTION ), 1);
    dim3  threads(WIDTH_PORTION, 1, 1);

	GPUuSrtToFloat<<<grid, threads>>>(input, output, width, height);
}

#define BLOCKDIM_X 8
#define BLOCKDIM_Y 8

__device__ int offset(int i, int j, unsigned short WIDTH, unsigned short HEIGHT)
{
	if(i < 0)       i = 0;
	if(i >= WIDTH)  i = WIDTH-1;
	if(j < 0)       j = 0;
	if(j >= HEIGHT) j = HEIGHT-1;
	return( j * WIDTH + i );
}

 //Trous2D en Y axe
__global__ void GPU_Trous2D_Y(unsigned short step, float* input, float* output, float* ptr, unsigned short WIDTH, unsigned short HEIGHT)
{

	const int ix = blockDim.x * blockIdx.x + threadIdx.x;
    const int iy = blockDim.y * blockIdx.y + threadIdx.y;

    if(ix < WIDTH && iy < HEIGHT)
	{
		if (iy < 2 * step)
			output[ iy * WIDTH + ix ] = 
				  coeffH0 * input[ iy * WIDTH + ix ]
				+ coeffH1 * (input[offset(ix, iy - step, WIDTH, HEIGHT)]     + input[offset(ix, iy + step, WIDTH, HEIGHT)])
				+ coeffH2 * (input[offset(ix, iy - 2 * step, WIDTH, HEIGHT)] + input[offset(ix, iy + 2 * step, WIDTH, HEIGHT)]);

		if (iy >= 2 * step && iy < HEIGHT - 2 * step)	
			output[ iy * WIDTH + ix ] = 
				  coeffH0 * input[ iy * WIDTH + ix ]
				+ coeffH1 * (input[(iy - step) * WIDTH + ix]     + input[(iy + step) * WIDTH + ix])
				+ coeffH2 * (input[(iy - 2 * step) * WIDTH + ix] + input[(iy + 2 * step) * WIDTH + ix]);

		if(iy >= HEIGHT - 2 * step)
			output[ iy * WIDTH + ix ] = 
				  coeffH0 * input[ iy * WIDTH + ix ]
				+ coeffH1 * (input[offset(ix, iy - step, WIDTH, HEIGHT)]     + input[offset(ix, iy + step, WIDTH, HEIGHT)])
				+ coeffH2 * (input[offset(ix, iy - 2 * step, WIDTH, HEIGHT)] + input[offset(ix, iy + 2 * step, WIDTH, HEIGHT)]);		
			
		ptr[  iy * WIDTH + ix  ] -= output[  iy * WIDTH + ix  ];	
	
	}
	
}

// Trous2D en X axe
__global__ void GPU_Trous2D_X(unsigned short step, float* input, float* output, unsigned short WIDTH, unsigned short HEIGHT)
{

	const int ix = blockDim.x * blockIdx.x + threadIdx.x;
    const int iy = blockDim.y * blockIdx.y + threadIdx.y;

    if(ix < WIDTH && iy < HEIGHT)
	{
		if (ix < 2 * step)
			output[ iy * WIDTH + ix ] = 
				  coeffH0 * input[ iy * WIDTH + ix ]
				+ coeffH1 * (input[offset(ix - step, iy, WIDTH, HEIGHT)]     + input[offset(ix + step, iy, WIDTH, HEIGHT)])
				+ coeffH2 * (input[offset(ix - 2 * step, iy, WIDTH, HEIGHT)] + input[offset(ix + 2 * step, iy, WIDTH, HEIGHT)]);

		if (ix >= 2 * step && ix < WIDTH - 2 * step)	
			output[ iy * WIDTH + ix ] = 
				  coeffH0 * input[ iy * WIDTH + ix ]
				+ coeffH1 * (input[iy * WIDTH + ix - step]     + input[iy * WIDTH + ix + step])
				+ coeffH2 * (input[iy * WIDTH + ix - 2 * step] + input[iy * WIDTH + ix + 2 * step]);

		if(ix >= WIDTH - 2 * step)
			output[ iy * WIDTH + ix ] = 
				  coeffH0 * input[ iy * WIDTH + ix ]
				+ coeffH1 * (input[offset(ix - step, iy, WIDTH, HEIGHT)]     + input[offset(ix + step, iy, WIDTH, HEIGHT)])
				+ coeffH2 * (input[offset(ix - 2 * step, iy, WIDTH, HEIGHT)] + input[offset(ix + 2 * step, iy, WIDTH, HEIGHT)]);		
	}

}


// projectionV02D en Y axe
__global__ void GPU_projectionV02D_Y(float* input, float* output, unsigned short WIDTH, unsigned short HEIGHT)
{
	// Chargement des données sur shared memory
    //__shared__ float sData [(HEIGHT_PORTION + 2) * WIDTH_PORTION];  /* 2 = templateSize/2 et l'autre 2 c'est pour les voisins haut-bas*/ 
	extern __shared__ float sData [];
    // numero d'image = blcokIdx.x / ( WIDTH / WIDTH_PORTION )

    unsigned short y0, yn;
	unsigned int index;
	unsigned short idImage = /* numero de l'image */ (blockIdx.x / (int)ceilf((float)WIDTH / (float)WIDTH_PORTION ));
	unsigned int imageSize = WIDTH * HEIGHT;
	unsigned short bx = blockIdx.x  % (int)ceilf((float)WIDTH / (float)WIDTH_PORTION);
	unsigned short H;

    y0 = 0;
    yn = HEIGHT_PORTION + 2;
	H = HEIGHT_PORTION;

	// Condition : width >= 32 and height >= 32 
    if ( blockIdx.y % (int)ceilf((float)HEIGHT / (float)HEIGHT_PORTION) == 0 ) 
	{
		y0 = 1;
		//if (HEIGHT % HEIGHT_PORTION == 0) H = HEIGHT_PORTION;
		//else H = (HEIGHT % HEIGHT_PORTION);
	}
	if ( (blockIdx.y+1) % (int)ceilf((float)HEIGHT / (float)HEIGHT_PORTION) == 0 ) 
	{
		if (HEIGHT % HEIGHT_PORTION == 0)
		{
			yn = yn - 1;
			H = HEIGHT_PORTION;
		}
		else 
		{
			yn = 1 + (HEIGHT % HEIGHT_PORTION);
			H = (HEIGHT % HEIGHT_PORTION);
		}
	}
	
	// Initialization of shared memory
    for(unsigned short y = y0; y < yn; y++)
    {
		index = ( threadIdx.x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + (y - 1) ) * WIDTH;
		if (( ( threadIdx.x + bx * WIDTH_PORTION) < WIDTH) && ( (blockIdx.y * HEIGHT_PORTION + (y - 1) ) < HEIGHT))
			sData[threadIdx.x + y * WIDTH_PORTION] = input[ idImage * imageSize + index ];
    }
   
    __syncthreads();

	for(unsigned short y = 0; y < H; y++)
	{
		index = ( threadIdx.x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + y ) * WIDTH;
		if (( ( threadIdx.x + bx * WIDTH_PORTION) < WIDTH) && ( (blockIdx.y * HEIGHT_PORTION + y ) < HEIGHT))
			output[ idImage * imageSize + index ] = 
				   _coeffH0 * sData[threadIdx.x + (y + 1) * WIDTH_PORTION]
				 + _coeffH1 * (sData[threadIdx.x + __min(y + 2, yn - 1)* WIDTH_PORTION] + sData[threadIdx.x + __max(y, y0)* WIDTH_PORTION]);
	}
}

 //projectionV02D en X axe
__global__ void GPU_projectionV02D_X(float* input, float* output, unsigned short WIDTH, unsigned short HEIGHT)
{
    // Chargement des données sur shared memory
    //__shared__ float sData [(WIDTH_PORTION + 2) * HEIGHT_PORTION];  /* 2 = templateSize/2 et l'autre 2 c'est pour les voisins gauche-droite */ 
	extern __shared__ float sData [];
    // numero d'image = blcokIdx.x / ( WIDTH / WIDTH_PORTION )

    unsigned short x0, xn;
	unsigned int index;
	unsigned short idImage = /* numero de l'image */ (blockIdx.x / (int)ceilf((float)WIDTH / (float)WIDTH_PORTION ));
	unsigned int imageSize = WIDTH * HEIGHT;
	unsigned short bx = blockIdx.x  % (int)ceilf((float)WIDTH / (float)WIDTH_PORTION);
	unsigned short W;

	x0 = 0;
    xn = WIDTH_PORTION + 2;
	W = WIDTH_PORTION;

	// Condition : width >= 32 and height >= 32 
    if ( bx == 0 ) 
	{
		x0 = 1;
		//if (WIDTH % WIDTH_PORTION == 0) W = WIDTH_PORTION;
		//else W = (WIDTH % WIDTH_PORTION);
	}
	if ( (blockIdx.x+1) % (int)ceilf((float)WIDTH / (float)WIDTH_PORTION) == 0 ) 
	{
		if (WIDTH % WIDTH_PORTION == 0) 
		{
			xn = xn - 1;
			W = WIDTH_PORTION;
		}
		else 
		{
			xn = 1 + (WIDTH % WIDTH_PORTION);
			W = (WIDTH % WIDTH_PORTION);
		}
	}

	// Initialization of shared memory
    for(unsigned short x = x0; x < xn; x++)
    {
		index = ( (x - 1) + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + threadIdx.y ) * WIDTH;
		if (( ( (x - 1) + bx * WIDTH_PORTION) < WIDTH) && ( (blockIdx.y * HEIGHT_PORTION + threadIdx.y ) < HEIGHT))
			sData[ x * HEIGHT_PORTION + threadIdx.y] = input[ idImage * imageSize +  index ];
    }
   
    __syncthreads();

	for(unsigned short x = 0; x < W; x++)
	{
		index  = ( x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + threadIdx.y ) * WIDTH;
		if (( ( x + bx * WIDTH_PORTION) < WIDTH) && ( (blockIdx.y * HEIGHT_PORTION + threadIdx.y ) < HEIGHT))
			output[ idImage * imageSize + index ] = 
				   _coeffH0 * sData[(x + 1) * HEIGHT_PORTION + threadIdx.y]
				 + _coeffH1 * (sData[__min(x + 2, xn - 1) * HEIGHT_PORTION + threadIdx.y] + sData[__max(x, x0) * HEIGHT_PORTION + threadIdx.y]);
	}
}

// Conversion : unsigned short TO float
__global__ void GPUuSrtToFloat(unsigned short* input, float* output, unsigned short WIDTH, unsigned short HEIGHT)
{
	unsigned int index;
	unsigned short idImage = /* numero de l'image */ (blockIdx.x / (int)ceilf((float)WIDTH / (float)WIDTH_PORTION ));
	unsigned int imageSize = WIDTH * HEIGHT;
	unsigned short bx = blockIdx.x  % (int)ceilf((float)WIDTH / (float)WIDTH_PORTION);

	for(unsigned short y = 0; y < HEIGHT_PORTION; y++)
	{
		index = ( ( threadIdx.x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + y ) * WIDTH );
		if ( ( ( threadIdx.x + bx * WIDTH_PORTION) < WIDTH) && ( (blockIdx.y * HEIGHT_PORTION + y ) < HEIGHT))
			output[ idImage * imageSize + index ] = (float)input[ idImage * imageSize + index ];
	}

}

extern "C"
void cudaProjectionV02DGPU(unsigned short* input, float* output, unsigned int numImages,unsigned short width, unsigned short height, bool realProj)
{ 
	//cutilSafeCall( cudaMemcpyToSymbol("HEIGHT", &height, sizeof(unsigned short), 0, cudaMemcpyHostToDevice) );
 //   cutilSafeCall( cudaMemcpyToSymbol("WIDTH", &width, sizeof(unsigned short), 0, cudaMemcpyHostToDevice) );
	
	// For debug
	//unsigned short a = 0;
	//cutilSafeCall( cudaMemcpyFromSymbol (&a, "WIDTH", sizeof(unsigned short), 0, cudaMemcpyDeviceToHost) );
	//printf("WIDTH : %d\n", a);
	//cutilSafeCall( cudaMemcpyFromSymbol (&a, "HEIGHT", sizeof(unsigned short), 0, cudaMemcpyDeviceToHost) );
	//printf("HEIGHT : %d\n", a);

	//cutilSafeCall( cudaMemcpy( &WIDTH, &width, sizeof(unsigned short), cudaMemcpyHostToDevice) );
	//cutilSafeCall( cudaMemcpy( &HEIGHT, &height, sizeof(unsigned short), cudaMemcpyHostToDevice) );
	
	// A completer
    dim3  grid( numImages * (unsigned int)ceil( (float)width / (float)WIDTH_PORTION ), (unsigned int)ceil( (float)height / (float)HEIGHT_PORTION ), 1);
    dim3  threads(WIDTH_PORTION, 1, 1);

	GPUuSrtToFloat<<<grid, threads>>>(input, output, width, height);

    // A completer
	if (realProj == true)
	{	
		float* buff;
		checkCudaErrors( cudaMalloc( (void**) &buff, numImages * width * height * sizeof(float) ) );

		// setup execution parameters
		grid.x = numImages * (unsigned int)ceil( (float)width / (float)WIDTH_PORTION );
		grid.y = (unsigned int)ceil( (float)height / (float)HEIGHT_PORTION );
		grid.z = 1;

		threads.x = 1;
		threads.y = HEIGHT_PORTION;
		threads.z = 1;

		// execute the kernels
		GPU_projectionV02D_X<<< grid, threads, (WIDTH_PORTION + 2) * HEIGHT_PORTION * sizeof(float) >>>(output, buff, width, height);

		threads.x = WIDTH_PORTION;
		threads.y = 1;

		// execute the kernels
		GPU_projectionV02D_Y<<< grid, threads, (HEIGHT_PORTION + 2) * WIDTH_PORTION * sizeof(float) >>>(buff, output, width, height);

		// cleanup memory
		checkCudaErrors(cudaFree(buff));
	}
}
	

extern "C" 
void cudaTrous2DGPU(float* input, float* output, unsigned int numImages, unsigned short width, unsigned short height, unsigned short step)
{
	float* buff;
	checkCudaErrors( cudaMalloc( (void**) &buff, numImages * width * height * sizeof(float) ) );
	
    dim3 threads(BLOCKDIM_X, BLOCKDIM_Y);
    dim3 grid((unsigned int)ceil((float)width / (float)BLOCKDIM_X), (unsigned int)ceil((float)height / (float)BLOCKDIM_Y));
	
	GPU_Trous2D_X<<< grid, threads >>>(step, input, buff, width, height);
	
	GPU_Trous2D_Y<<< grid, threads >>>(step, buff, output, input, width, height);
	
	// cleanup memory
    checkCudaErrors(cudaFree(buff));

}
#endif // #ifndef _DECOMPOSITION_KERNEL_H_