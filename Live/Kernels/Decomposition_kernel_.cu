#ifndef _DECOMPOSITION_KERNEL_H_
#define _DECOMPOSITION_KERNEL_H_

//#include <math_constants.h>
#include <helper_cuda.h>
#include <helper_timer.h>

#include "../parameters.h"

__global__ void GPUuSrtToFloat(unsigned short* input, float* output, unsigned short WIDTH, unsigned short HEIGHT);


extern "C"
void convertUSrtToFloat(unsigned short* input, float* output, unsigned short numImages, unsigned short width, unsigned short height)
{
	dim3  grid( numImages * (unsigned int)ceil( (float)width / (float)WIDTH_PORTION ), (unsigned int)ceil( (float)height / (float)HEIGHT_PORTION ), 1);
    dim3  threads(WIDTH_PORTION, 1, 1);

	GPUuSrtToFloat<<<grid, threads>>>(input, output, width, height);
}


 //Trous2D en Y axe
__global__ void GPU_Trous2D_Y(unsigned short step, float* input, float* output, float* ptr, unsigned short WIDTH, unsigned short HEIGHT)
{
	// Chargement des données sur shared memory
    //__shared__ float sData [(HEIGHT_PORTION + 2 * 2 * step) * WIDTH_PORTION];  /* 2 = templateSize/2 et l'autre 2 c'est pour les voisins haut-bas*/ 
	extern __shared__ float sData [];
    // numero d'image = blcokIdx.x / ( WIDTH / WIDTH_PORTION )

    unsigned short y0, yn;
	unsigned int index;
	unsigned short idImage = /* numero de l'image */ (blockIdx.x / (int)ceilf((float)WIDTH / (float)WIDTH_PORTION ));
	unsigned int imageSize = WIDTH * HEIGHT;
	unsigned short bx = blockIdx.x  % (int)ceilf((float)WIDTH / (float)WIDTH_PORTION);
	unsigned short H;

    y0 = 0;
    yn = HEIGHT_PORTION + 2 * 2 * step;
	H = HEIGHT_PORTION;

	// Condition : width >= 32 and height >= 32 
    if ( blockIdx.y % (int)ceilf((float)HEIGHT / (float)HEIGHT_PORTION) == 0 ) 
	{
		y0 = 2 * step;
		//if (HEIGHT % HEIGHT_PORTION == 0) H = HEIGHT_PORTION;
		//else H = (HEIGHT % HEIGHT_PORTION);
	}
	if ( (blockIdx.y+1) % (int)ceilf((float)HEIGHT / (float)HEIGHT_PORTION) == 0 ) 
	{
		if (HEIGHT % HEIGHT_PORTION == 0)
		{
			yn = yn - 2 * step;
			H = HEIGHT_PORTION;
		}
		else 
		{
			yn = 2 * step + (HEIGHT % HEIGHT_PORTION);
			H = (HEIGHT % HEIGHT_PORTION);
		}
	}

	// Initialization of shared memory
    for(unsigned short y = y0; y < yn; y++)
	{
		index = ( threadIdx.x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + (y - 2 * step) ) * WIDTH;
		if ( (( threadIdx.x + bx * WIDTH_PORTION) < WIDTH) && ((blockIdx.y * HEIGHT_PORTION + (y - 2 * step)) < HEIGHT) )
			sData[threadIdx.x + y * WIDTH_PORTION] = input[ idImage * imageSize + index ];
	}
   
    __syncthreads();

	for(unsigned short y = 0; y < H; y++)
	{
		index  = ( threadIdx.x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + y ) * WIDTH;
		if ( (( threadIdx.x + bx * WIDTH_PORTION) < WIDTH) && ((blockIdx.y * HEIGHT_PORTION + y ) < HEIGHT) )
			output[ idImage * imageSize + index ] =
				   coeffH0 * sData[threadIdx.x + (y + 2 * step) * WIDTH_PORTION]
				 + coeffH1 * (sData[threadIdx.x + __max(y + step, y0) * WIDTH_PORTION] + sData[threadIdx.x + __min(y + 3 * step, yn - 1) * WIDTH_PORTION])
				 + coeffH2 * (sData[threadIdx.x + __max(y, y0) * WIDTH_PORTION] + sData[threadIdx.x + __min(y + 4 * step, yn - 1) * WIDTH_PORTION]);
	}


	for(unsigned short y = 0; y < H; y++)
	{
		index = ( threadIdx.x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + y ) * WIDTH;
		if ( (( threadIdx.x + bx * WIDTH_PORTION) < WIDTH) && ((blockIdx.y * HEIGHT_PORTION + y ) < HEIGHT) )
			ptr[ idImage * imageSize + index ] -= output[ idImage * imageSize + index ];
	}
}

// Trous2D en X axe
__global__ void GPU_Trous2D_X(unsigned short step, float* input, float* output, unsigned short WIDTH, unsigned short HEIGHT)
{
    // Chargement des données sur shared memory
    //__shared__ float sData [(WIDTH_PORTION + 2 * 2 * step) * HEIGHT_PORTION];  /* 2 = templateSize/2 et l'autre 2 c'est pour les voisins gauche-droite */ 
	extern __shared__ float sData [];
    // numero d'image = blcokIdx.x / ( WIDTH / WIDTH_PORTION )

    unsigned short x0, xn;
	unsigned int index;
	unsigned short idImage = (blockIdx.x / (int)ceilf((float)WIDTH / (float)WIDTH_PORTION )); /* numero de l'image */ 
	unsigned int imageSize = WIDTH * HEIGHT;
	unsigned short bx = blockIdx.x  % (int)ceilf((float)WIDTH / (float)WIDTH_PORTION);
	unsigned short W;

	x0 = 0;
    xn = WIDTH_PORTION + 2 * 2 * step;
	W = WIDTH_PORTION;

	// Condition : width >= 32 and height >= 32 
    if ( bx == 0 ) 
	{
		x0 = 2 * step;
		//if (WIDTH % WIDTH_PORTION == 0) W = WIDTH_PORTION;
		//else W = (WIDTH % WIDTH_PORTION);
	}
	if ( (blockIdx.x+1) % (int)ceilf((float)WIDTH / (float)WIDTH_PORTION) == 0 ) 
	{
		if (WIDTH % WIDTH_PORTION == 0) 
		{
			xn = xn - 2 * step;
			W = WIDTH_PORTION;
		}
		else 
		{
			xn = 2 * step + (WIDTH % WIDTH_PORTION)/* + 1*/;
			W = (WIDTH % WIDTH_PORTION);
		}
	}

	// Initialization of shared memory
    for(unsigned short x = x0; x < xn; x++)
	{
		index = ( (x - 2 * step) + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + threadIdx.y ) * WIDTH ;
		if (( ( (x - 2 * step) + bx * WIDTH_PORTION) < WIDTH) && ( (blockIdx.y * HEIGHT_PORTION + threadIdx.y ) < HEIGHT))
			sData[ x * HEIGHT_PORTION + threadIdx.y] = input[ idImage * imageSize + index ];
	}
   
    __syncthreads();
	
	for(unsigned short x = 0; x < W; x++)
	{
		index = ( x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + threadIdx.y ) * WIDTH;
		if (( ( x + bx * WIDTH_PORTION) < WIDTH) && ( (blockIdx.y * HEIGHT_PORTION + threadIdx.y ) < HEIGHT))
			output[ idImage * imageSize + index ] = 
				   coeffH0 * sData[(x + 2 * step) * HEIGHT_PORTION + threadIdx.y]
				 + coeffH1 * (sData[__max(x + step, x0) * HEIGHT_PORTION + threadIdx.y] + sData[__min(x + 3 * step, xn - 1) * HEIGHT_PORTION + threadIdx.y])
				 + coeffH2 * (sData[__max(x, x0) * HEIGHT_PORTION + threadIdx.y] + sData[__min(x + 4 * step, xn - 1) * HEIGHT_PORTION + threadIdx.y]);		
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
	//checkCudaErrors( cudaMemcpy( &WIDTH, &width, sizeof(unsigned short), cudaMemcpyHostToDevice) );
	//checkCudaErrors( cudaMemcpy( &HEIGHT, &height, sizeof(unsigned short), cudaMemcpyHostToDevice) );

	float* buff;
	checkCudaErrors( cudaMalloc( (void**) &buff, numImages * width * height * sizeof(float) ) );

	// execute the kernels
    dim3  grid( numImages * (unsigned int)ceil( (float)width / (float)WIDTH_PORTION ), (unsigned int)ceil( (float)height / (float)HEIGHT_PORTION ), 1);
    dim3  threads(1, HEIGHT_PORTION, 1);
    GPU_Trous2D_X<<< grid, threads, (WIDTH_PORTION + 2 * 2 * step) * HEIGHT_PORTION * sizeof(float) >>>(step, input, buff, width, height);

	// execute the kernels
	threads.x = WIDTH_PORTION;
	threads.y = 1;
    GPU_Trous2D_Y<<< grid, threads, (HEIGHT_PORTION + 2 * 2 * step) * WIDTH_PORTION * sizeof(float) >>>(step, buff, output, input, width, height);

	// cleanup memory
    checkCudaErrors(cudaFree(buff));
}

#endif // #ifndef _DECOMPOSITION_KERNEL_H_