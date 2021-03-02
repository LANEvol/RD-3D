#include <iostream>
#include <stdio.h>
#include <set>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <ctime>

// includes, cuda
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

// Utilities and timing functions
#include <cuda.h>    // includes cuda.h and cuda_runtime_api.h
#include <cuda_runtime_api.h>
//#include <timer.h>               // timing functions

// CUDA helper functions
//#include <helper_cuda.h>         // helper functions for CUDA error check 
#include <curand.h>
#include <map>
#include <utility>



using namespace std;

GLuint gl_u, gl_v, gl_points;

double *d_u_old, *d_u_new, *d_v_old, *d_v_new, *d_w_old, *d_w_new;
float  *border, *d_border, *d_z, *z, *Deff, *d_Deff;

double *u_new, *v_new, *w_new, *u_old, *v_old, *w_old, *centers, *pars;

int  *closest_el, *neighbours;
bool *initialColors;
bool periodic, constantDiffusion, correction, growth, incrementE;

map<int, vector<int> >  contributions;
map<int, vector<int> >  nodesPerHexa;

int nodesSize, a, b, c, centersSize, nodesToSave, networkType;
int  saveFrequency, iteration, saveIteration, quasiStaticStep;
int stop, saving;
double h, dt, Du, Dv, Dw, U, V, W, cu, cV, cw, c1, c2, c3, c4, c5, c6, c7, c8, c9, error, hz,P, bt, border_thickness_in_elements,sigma, minDiffusionFactor;
double initialGreen[3];
double initialBlack[3];
double initialUniform[3];
double sx, sy;
char saveFolder[200];


dim3 DimBlock;
dim3 DimGrid;
dim3 DimBlockSim;
dim3 DimGridSim;



__constant__ int dimensions[3];
__constant__ double parameters3var[20];
texture <float, 1, cudaReadModeElementType> borderTexture;
texture <float, 1, cudaReadModeElementType> zTexture;
texture <float, 1, cudaReadModeElementType> DeffTexture;



#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


__device__ double reaction_u(double u1, double v1, double w1)
{
	double U = parameters3var[2] * v1 + parameters3var[3] * w1 + parameters3var[4];
	if(U > parameters3var[17])
	{
		return parameters3var[17];
	}
	else if(U < 0)
	{
		return 0;
	}
	return U;
}

__device__ double reaction_v(double u1, double v1, double w1)
{
	double U = parameters3var[5] * u1 + parameters3var[6] * w1 + parameters3var[7];
	if(U > parameters3var[18])
	{
		return parameters3var[18];
	}
	else if(U < 0)
	{
		return 0;
	}
	return U;
}

__device__ double reaction_w(double u1, double v1, double w1)
{
	double U = parameters3var[8] * u1 + parameters3var[9] * v1 + parameters3var[10];
	if(U > parameters3var[19])
	{
		return parameters3var[19];
	}
	else if(U < 0)
	{
		return 0;
	}
	return U;
}



__global__ void mainFunction3D3varBorder_P(double *u_old, double *u_new, double *v_old, double *v_new, double *w_old, double *w_new)
{
	__shared__ double u[16][16][4];
	__shared__ double v[16][16][4];
	__shared__ double w[16][16][4];
	int myIndX = blockIdx.x * (blockDim.x - 2) + threadIdx.x - 1;
	int myIndY = blockIdx.y * (blockDim.y - 2) + threadIdx.y - 1;
	int myIndZ = blockIdx.z * (blockDim.z - 2) + threadIdx.z - 1;
	int myIndXLocal = threadIdx.x;
	int myIndYLocal = threadIdx.y;
	int myIndZLocal = threadIdx.z;
	int index = 0;
	bool go = false;


	if(myIndX < (dimensions[0] + 1) && myIndY < (dimensions[1] + 1) && myIndZ < (dimensions[2] + 1))
	{
		if(myIndX >=0 && myIndX <= dimensions[0] - 1  && myIndY >= 0 && myIndY <= dimensions[1] - 1 && myIndZ >= 0 && myIndZ <= dimensions[2] - 1)
		{
			go = true;
		}
		if(myIndX < 0) myIndX = dimensions[0] - 1;
		if(myIndY < 0) myIndY = dimensions[1] - 1;
		if(myIndZ < 0) myIndZ = 0; 
		if(myIndX == dimensions[0]) myIndX = 0;
		if(myIndY == dimensions[1]) myIndY  = 0;
		if(myIndZ == dimensions[2]) myIndZ  = dimensions[2] - 1;
		//global index in 1D array
		index = myIndZ * (dimensions[0] * dimensions[1]) + myIndX * dimensions[1] + myIndY;

		u[myIndXLocal][myIndYLocal][myIndZLocal] = u_old[index];
		v[myIndXLocal][myIndYLocal][myIndZLocal] = v_old[index];
		w[myIndXLocal][myIndYLocal][myIndZLocal] = w_old[index];
	}
	//solve for the indices out of domain
	__syncthreads();
	if(go && myIndXLocal > 0 && myIndXLocal < 15 && myIndYLocal > 0 && myIndYLocal < 15 && myIndZLocal > 0 && myIndZLocal < 3)
	{

		float xl = tex1Dfetch(borderTexture, 6 * index);
		float xr = tex1Dfetch(borderTexture, 6 * index + 1);
		float yl = tex1Dfetch(borderTexture, 6 * index + 2);
		float yr = tex1Dfetch(borderTexture, 6 * index + 3);
		float zl = tex1Dfetch(borderTexture, 6 * index + 4);
		float zr = tex1Dfetch(borderTexture, 6 * index + 5);
		
		//check if the node has to be simulated
		if(xl < -1) return;

		double Laplacian_u = 0;
		double Laplacian_v = 0;
		double Laplacian_w = 0;
		
		Laplacian_u += xl * (u[myIndXLocal - 1][myIndYLocal][myIndZLocal] - u[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_v += xl * (v[myIndXLocal - 1][myIndYLocal][myIndZLocal] - v[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_w += xl * (w[myIndXLocal - 1][myIndYLocal][myIndZLocal] - w[myIndXLocal][myIndYLocal][myIndZLocal]);	
		
		Laplacian_u += xr * (u[myIndXLocal + 1][myIndYLocal][myIndZLocal] - u[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_v += xr * (v[myIndXLocal + 1][myIndYLocal][myIndZLocal] - v[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_w += xr * (w[myIndXLocal + 1][myIndYLocal][myIndZLocal] - w[myIndXLocal][myIndYLocal][myIndZLocal]);		
		
		Laplacian_u += yl * (u[myIndXLocal][myIndYLocal - 1][myIndZLocal] - u[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_v += yl * (v[myIndXLocal][myIndYLocal - 1][myIndZLocal] - v[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_w += yl * (w[myIndXLocal][myIndYLocal - 1][myIndZLocal] - w[myIndXLocal][myIndYLocal][myIndZLocal]);		
		
		Laplacian_u += yr * (u[myIndXLocal][myIndYLocal + 1][myIndZLocal] - u[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_v += yr * (v[myIndXLocal][myIndYLocal + 1][myIndZLocal] - v[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_w += yr * (w[myIndXLocal][myIndYLocal + 1][myIndZLocal] - w[myIndXLocal][myIndYLocal][myIndZLocal]);			
	
		Laplacian_u += zl * (u[myIndXLocal][myIndYLocal][myIndZLocal - 1] - u[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_v += zl * (v[myIndXLocal][myIndYLocal][myIndZLocal - 1] - v[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_w += zl * (w[myIndXLocal][myIndYLocal][myIndZLocal - 1] - w[myIndXLocal][myIndYLocal][myIndZLocal]);

		Laplacian_u += zr * (u[myIndXLocal][myIndYLocal][myIndZLocal + 1] - u[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_v += zr * (v[myIndXLocal][myIndYLocal][myIndZLocal + 1] - v[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_w += zr * (w[myIndXLocal][myIndYLocal][myIndZLocal + 1] - w[myIndXLocal][myIndYLocal][myIndZLocal]);
		
		Laplacian_u = u[myIndXLocal][myIndYLocal][myIndZLocal] + parameters3var[1] * (Laplacian_u * parameters3var[14]/parameters3var[0] + reaction_u(u[myIndXLocal][myIndYLocal][myIndZLocal], v[myIndXLocal][myIndYLocal][myIndZLocal], w[myIndXLocal][myIndYLocal][myIndZLocal]) - parameters3var[11] * u[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_v = v[myIndXLocal][myIndYLocal][myIndZLocal] + parameters3var[1] * (Laplacian_v * parameters3var[15]/parameters3var[0] + reaction_v(u[myIndXLocal][myIndYLocal][myIndZLocal], v[myIndXLocal][myIndYLocal][myIndZLocal], w[myIndXLocal][myIndYLocal][myIndZLocal]) - parameters3var[12] * v[myIndXLocal][myIndYLocal][myIndZLocal]);
		Laplacian_w = w[myIndXLocal][myIndYLocal][myIndZLocal] + parameters3var[1] * (Laplacian_w * parameters3var[16]/parameters3var[0] + reaction_w(u[myIndXLocal][myIndYLocal][myIndZLocal], v[myIndXLocal][myIndYLocal][myIndZLocal], w[myIndXLocal][myIndYLocal][myIndZLocal]) - parameters3var[13] * w[myIndXLocal][myIndYLocal][myIndZLocal]);

		u_new[index] = Laplacian_u > 0 ? Laplacian_u : 0;
		v_new[index] = Laplacian_v > 0 ? Laplacian_v : 0;
		w_new[index] = Laplacian_w > 0 ? Laplacian_w : 0;
	}
}

//it works in 2D only for a squared lattice, otherwise call 3D method
__global__ void mainFunction2D3varBorder_P(double *u_old, double *u_new, double *v_old, double *v_new, double *w_old, double *w_new)
{
	__shared__ double u[32][32];
	__shared__ double v[32][32];
	__shared__ double w[32][32];
	int myIndX = blockIdx.x * (blockDim.x - 2) + threadIdx.x - 1;
	int myIndY = blockIdx.y * (blockDim.y - 2) + threadIdx.y - 1;
	int myIndXLocal = threadIdx.x;
	int myIndYLocal = threadIdx.y;
	int index = 0;
	bool go = false;


	if(myIndX < (dimensions[0] + 1) && myIndY < (dimensions[1] + 1))
	{
		if(myIndX >=0 && myIndX <= dimensions[0] - 1  && myIndY >= 0 && myIndY <= dimensions[1] - 1)
		{
			go = true;
		}
		if(myIndX < 0) myIndX = 0;
		if(myIndY < 0) myIndY = 0;
		if(myIndX == dimensions[0]) myIndX = dimensions[0] - 1;
		if(myIndY == dimensions[1]) myIndY  = dimensions[1] - 1;
		//global index in 1D array
		index = myIndX * dimensions[1] + myIndY;

		u[myIndXLocal][myIndYLocal] = u_old[index];
		v[myIndXLocal][myIndYLocal] = v_old[index];
		w[myIndXLocal][myIndYLocal] = w_old[index];
	}
	//solve for the indices out of domain
	__syncthreads();
	if(go && myIndXLocal > 0 && myIndXLocal < 31 && myIndYLocal > 0 && myIndYLocal < 31)
	{


		float xl = tex1Dfetch(borderTexture, 6 * index);
		float xr = tex1Dfetch(borderTexture, 6 * index + 1);
		float yl = tex1Dfetch(borderTexture, 6 * index + 2);
		float yr = tex1Dfetch(borderTexture, 6 * index + 3);
		double Laplacian_u = 0;
		double Laplacian_v = 0;
		double Laplacian_w = 0;
		
		if(xl < -1) return;
		
		Laplacian_u += xl * (u[myIndXLocal - 1][myIndYLocal] - u[myIndXLocal][myIndYLocal]);
		Laplacian_v += xl * (v[myIndXLocal - 1][myIndYLocal] - v[myIndXLocal][myIndYLocal]);
		Laplacian_w += xl * (w[myIndXLocal - 1][myIndYLocal] - w[myIndXLocal][myIndYLocal]);	
		
		Laplacian_u += xr * (u[myIndXLocal + 1][myIndYLocal] - u[myIndXLocal][myIndYLocal]);
		Laplacian_v += xr * (v[myIndXLocal + 1][myIndYLocal] - v[myIndXLocal][myIndYLocal]);
		Laplacian_w += xr * (w[myIndXLocal + 1][myIndYLocal] - w[myIndXLocal][myIndYLocal]);		
		
		Laplacian_u += yl * (u[myIndXLocal][myIndYLocal - 1] - u[myIndXLocal][myIndYLocal]);
		Laplacian_v += yl * (v[myIndXLocal][myIndYLocal - 1] - v[myIndXLocal][myIndYLocal]);
		Laplacian_w += yl * (w[myIndXLocal][myIndYLocal - 1] - w[myIndXLocal][myIndYLocal]);		
		
		Laplacian_u += yr * (u[myIndXLocal][myIndYLocal + 1] - u[myIndXLocal][myIndYLocal]);
		Laplacian_v += yr * (v[myIndXLocal][myIndYLocal + 1] - v[myIndXLocal][myIndYLocal]);
		Laplacian_w += yr * (w[myIndXLocal][myIndYLocal + 1] - w[myIndXLocal][myIndYLocal]);			
		
		Laplacian_u = u[myIndXLocal][myIndYLocal] + parameters3var[1] * (Laplacian_u * parameters3var[14]/parameters3var[0] + reaction_u(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[11] * u[myIndXLocal][myIndYLocal]);
		Laplacian_v = v[myIndXLocal][myIndYLocal] + parameters3var[1] * (Laplacian_v * parameters3var[15]/parameters3var[0] + reaction_v(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[12] * v[myIndXLocal][myIndYLocal]);
		Laplacian_w = w[myIndXLocal][myIndYLocal] + parameters3var[1] * (Laplacian_w * parameters3var[16]/parameters3var[0] + reaction_w(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[13] * w[myIndXLocal][myIndYLocal]);

		u_new[index] = Laplacian_u > 0 ? Laplacian_u : 0;
		v_new[index] = Laplacian_v > 0 ? Laplacian_v : 0;
		w_new[index] = Laplacian_w > 0 ? Laplacian_w : 0;
	}
}

//it works in 2D only for a squared lattice
//the D texture contains the heights at each point
//the arrays u,v, and w contain the concentrations is 2D (multiplied by the height of the system)
__global__ void mainFunction2D3varBorder_D(double *u_old, double *u_new, double *v_old, double *v_new, double *w_old, double *w_new)
{
	__shared__ double u[16][16];
	__shared__ double v[16][16];
	__shared__ double w[16][16];
	int myIndX = blockIdx.x * (blockDim.x - 2) + threadIdx.x - 1;
	int myIndY = blockIdx.y * (blockDim.y - 2) + threadIdx.y - 1;
	int myIndXLocal = threadIdx.x;
	int myIndYLocal = threadIdx.y;
	int index = 0;
	bool go = false;


	if(myIndX < (dimensions[0] + 1) && myIndY < (dimensions[1] + 1))
	{
		if(myIndX >=0 && myIndX <= dimensions[0] - 1  && myIndY >= 0 && myIndY <= dimensions[1] - 1)
		{
			go = true;
		}
		if(myIndX < 0) myIndX = 0;
		if(myIndY < 0) myIndY = 0;
		if(myIndX == dimensions[0]) myIndX = dimensions[0] - 1;
		if(myIndY == dimensions[1]) myIndY  = dimensions[1] - 1;
		//global index in 1D array
		index = myIndX * dimensions[1] + myIndY;

		u[myIndXLocal][myIndYLocal] = u_old[index];
		v[myIndXLocal][myIndYLocal] = v_old[index];
		w[myIndXLocal][myIndYLocal] = w_old[index];
	}
	//solve for the indices out of domain
	__syncthreads();
	if(go && myIndXLocal > 0 && myIndXLocal < 15 && myIndYLocal > 0 && myIndYLocal < 15)
	{

		float z = tex1Dfetch(zTexture, 3 * index);
		float dzdx = tex1Dfetch(zTexture, 3 * index + 1);
		float dzdy = tex1Dfetch(zTexture, 3 * index + 2);

		
		double Laplacian_u = z * (u[myIndXLocal - 1][myIndYLocal] + u[myIndXLocal + 1][myIndYLocal] + u[myIndXLocal][myIndYLocal - 1] + u[myIndXLocal][myIndYLocal + 1] - 4 * u[myIndXLocal][myIndYLocal])/(parameters3var[0]*parameters3var[0]);
		double Laplacian_v = z * (v[myIndXLocal - 1][myIndYLocal] + v[myIndXLocal + 1][myIndYLocal] + v[myIndXLocal][myIndYLocal - 1] + v[myIndXLocal][myIndYLocal + 1] - 4 * v[myIndXLocal][myIndYLocal])/(parameters3var[0]*parameters3var[0]);
		double Laplacian_w = z * (w[myIndXLocal - 1][myIndYLocal] + w[myIndXLocal + 1][myIndYLocal] + w[myIndXLocal][myIndYLocal - 1] + w[myIndXLocal][myIndYLocal + 1] - 4 * w[myIndXLocal][myIndYLocal])/(parameters3var[0]*parameters3var[0]);	
		
		Laplacian_u += dzdx * (u[myIndXLocal + 1][myIndYLocal] - u[myIndXLocal - 1][myIndYLocal])/(2*parameters3var[0]);
		Laplacian_v += dzdx * (v[myIndXLocal + 1][myIndYLocal] - v[myIndXLocal - 1][myIndYLocal])/(2*parameters3var[0]);		
		Laplacian_w += dzdx * (w[myIndXLocal + 1][myIndYLocal] - w[myIndXLocal - 1][myIndYLocal])/(2*parameters3var[0]);		
		
		Laplacian_u += dzdy * (u[myIndXLocal][myIndYLocal + 1] - u[myIndXLocal][myIndYLocal - 1])/(2*parameters3var[0]);
		Laplacian_v += dzdy * (v[myIndXLocal][myIndYLocal + 1] - v[myIndXLocal][myIndYLocal - 1])/(2*parameters3var[0]);		
		Laplacian_w += dzdy * (w[myIndXLocal][myIndYLocal + 1] - w[myIndXLocal][myIndYLocal - 1])/(2*parameters3var[0]);		

		Laplacian_u *= parameters3var[14];
		Laplacian_v *= parameters3var[15];
		Laplacian_w *= parameters3var[16];

		
		Laplacian_u = u[myIndXLocal][myIndYLocal] + parameters3var[1]/z * (Laplacian_u + (reaction_u(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[11] * u[myIndXLocal][myIndYLocal])*z);
		Laplacian_v = v[myIndXLocal][myIndYLocal] + parameters3var[1]/z * (Laplacian_v + (reaction_v(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[12] * v[myIndXLocal][myIndYLocal])*z);
		Laplacian_w = w[myIndXLocal][myIndYLocal] + parameters3var[1]/z * (Laplacian_w + (reaction_w(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[13] * w[myIndXLocal][myIndYLocal])*z);
	
		u_new[index] = Laplacian_u;
		v_new[index] = Laplacian_v;
		w_new[index] = Laplacian_w;
	}
}

//it works in 2D only for a squared lattice
//the D texture contains the heights at each point
//the arrays u,v, and w contain the concentrations is 2D (multiplied by the height of the system)
__global__ void mainFunction2D3varBorder_D_Correction(double *u_old, double *u_new, double *v_old, double *v_new, double *w_old, double *w_new)
{
	__shared__ double u[16][16];
	__shared__ double v[16][16];
	__shared__ double w[16][16];
	int myIndX = blockIdx.x * (blockDim.x - 2) + threadIdx.x - 1;
	int myIndY = blockIdx.y * (blockDim.y - 2) + threadIdx.y - 1;
	int myIndXLocal = threadIdx.x;
	int myIndYLocal = threadIdx.y;
	int index = 0;
	bool go = false;


	if(myIndX < (dimensions[0] + 1) && myIndY < (dimensions[1] + 1))
	{
		if(myIndX >=0 && myIndX <= dimensions[0] - 1  && myIndY >= 0 && myIndY <= dimensions[1] - 1)
		{
			go = true;
		}
		if(myIndX < 0) myIndX = 0;
		if(myIndY < 0) myIndY = 0;
		if(myIndX == dimensions[0]) myIndX = dimensions[0] - 1;
		if(myIndY == dimensions[1]) myIndY  = dimensions[1] - 1;
		//global index in 1D array
		index = myIndX * dimensions[1] + myIndY;

		u[myIndXLocal][myIndYLocal] = u_old[index];
		v[myIndXLocal][myIndYLocal] = v_old[index];
		w[myIndXLocal][myIndYLocal] = w_old[index];
	}
	//solve for the indices out of domain
	__syncthreads();
	if(go && myIndXLocal > 0 && myIndXLocal < 15 && myIndYLocal > 0 && myIndYLocal < 15)
	{

		float z = tex1Dfetch(zTexture, 3 * index);
		float dzdx = tex1Dfetch(zTexture, 3 * index + 1);
		float dzdy = tex1Dfetch(zTexture, 3 * index + 2);
		float Dx = tex1Dfetch(DeffTexture, 2 * index);
		float Dy = tex1Dfetch(DeffTexture, 2 * index + 1);

		
		double Laplacian_u = Dx * (u[myIndXLocal - 1][myIndYLocal] + u[myIndXLocal + 1][myIndYLocal] - 2 * u[myIndXLocal][myIndYLocal]);
		double Laplacian_v = Dx * (v[myIndXLocal - 1][myIndYLocal] + v[myIndXLocal + 1][myIndYLocal] - 2 * v[myIndXLocal][myIndYLocal]);
		double Laplacian_w = Dx * (w[myIndXLocal - 1][myIndYLocal] + w[myIndXLocal + 1][myIndYLocal] - 2 * w[myIndXLocal][myIndYLocal]);
		
		Laplacian_u += Dy * (u[myIndXLocal][myIndYLocal - 1] + u[myIndXLocal][myIndYLocal + 1] - 2 * u[myIndXLocal][myIndYLocal]);
		Laplacian_v += Dy * (v[myIndXLocal][myIndYLocal - 1] + v[myIndXLocal][myIndYLocal + 1] - 2 * v[myIndXLocal][myIndYLocal]);
		Laplacian_w += Dy * (w[myIndXLocal][myIndYLocal - 1] + w[myIndXLocal][myIndYLocal + 1] - 2 * w[myIndXLocal][myIndYLocal]);
		
		Laplacian_u *= z/(parameters3var[0] * parameters3var[0]);
		Laplacian_v *= z/(parameters3var[0] * parameters3var[0]);
		Laplacian_w *= z/(parameters3var[0] * parameters3var[0]);
						
		Laplacian_u += dzdx * (u[myIndXLocal + 1][myIndYLocal] - u[myIndXLocal - 1][myIndYLocal])/(2*parameters3var[0]);
		Laplacian_v += dzdx * (v[myIndXLocal + 1][myIndYLocal] - v[myIndXLocal - 1][myIndYLocal])/(2*parameters3var[0]);		
		Laplacian_w += dzdx * (w[myIndXLocal + 1][myIndYLocal] - w[myIndXLocal - 1][myIndYLocal])/(2*parameters3var[0]);		
		
		Laplacian_u += dzdy * (u[myIndXLocal][myIndYLocal + 1] - u[myIndXLocal][myIndYLocal - 1])/(2*parameters3var[0]);
		Laplacian_v += dzdy * (v[myIndXLocal][myIndYLocal + 1] - v[myIndXLocal][myIndYLocal - 1])/(2*parameters3var[0]);		
		Laplacian_w += dzdy * (w[myIndXLocal][myIndYLocal + 1] - w[myIndXLocal][myIndYLocal - 1])/(2*parameters3var[0]);		

		Laplacian_u *= parameters3var[14];
		Laplacian_v *= parameters3var[15];
		Laplacian_w *= parameters3var[16];

		
		Laplacian_u = u[myIndXLocal][myIndYLocal] + parameters3var[1]/z * (Laplacian_u + (reaction_u(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[11] * u[myIndXLocal][myIndYLocal])*z);
		Laplacian_v = v[myIndXLocal][myIndYLocal] + parameters3var[1]/z * (Laplacian_v + (reaction_v(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[12] * v[myIndXLocal][myIndYLocal])*z);
		Laplacian_w = w[myIndXLocal][myIndYLocal] + parameters3var[1]/z * (Laplacian_w + (reaction_w(u[myIndXLocal][myIndYLocal], v[myIndXLocal][myIndYLocal], w[myIndXLocal][myIndYLocal]) - parameters3var[13] * w[myIndXLocal][myIndYLocal])*z);
	
		u_new[index] = Laplacian_u;
		v_new[index] = Laplacian_v;
		w_new[index] = Laplacian_w;
	}
}

void initCuda(bool random)
{
	//dimensions of the system
	int *dim = (int*)malloc(3 * sizeof(int));
	dim[0] = a;
	dim[1] = b; 
	dim[2] = c;
	cudaMemcpyToSymbol(dimensions, dim, 3 * sizeof(int));
	free(dim);

	pars = (double*)malloc(20 * sizeof(double));
	pars[0] = h; pars[1] = dt; pars[2] = c1; pars[3] = c2; pars[4] = c3; pars[5] = c4;
	pars[6] = c5; pars[7] = c6; pars[8] = c7;
	pars[9] = c8; pars[10] = c9; pars[11] = cu;
	pars[12] = cV; pars[13] = cw; pars[14] = Du;
	pars[15] = Dv; pars[16] = Dw; pars[17] = U;
	pars[18] = V; pars[19] = W;
	cudaMemcpyToSymbol(parameters3var, pars, 20 * sizeof(double));

	//set block dimensions (32,32, 1) for 2D and (16,16,4) for 3D
	if(c == 1 && constantDiffusion) //we simulate the P value at the border
	{
		DimBlock = dim3(32,32,1);
		DimGridSim = dim3(a/30 + 1 * (bool)(a % 30), b/30 + 1 * (bool)(b % 30), 1);
	}
	else
	{
		if(c > 1)
		{	
			DimBlock = dim3(16,16,4);
			DimGridSim = dim3(a/14 + 1 * (bool)(a % 14), b/14 + 1 * (bool)(b % 14), c/2 + 1 * (bool)(c % 2));				
		}
		else
		{
			DimBlock = dim3(16,16,1);
			DimGridSim = dim3(a/14 + 1 * (bool)(a % 14), b/14 + 1 * (bool)(b % 14),1);			
		}
	}

	cudaMalloc((void**) &d_u_old, nodesSize*sizeof(double));
	cudaMalloc((void**) &d_u_new, nodesSize*sizeof(double));
	cudaMalloc((void**) &d_v_new, nodesSize*sizeof(double));
	cudaMalloc((void**) &d_v_old, nodesSize*sizeof(double));
	cudaMalloc((void**) &d_w_new, nodesSize*sizeof(double));
	cudaMalloc((void**) &d_w_old, nodesSize*sizeof(double));
	cudaMalloc((void**) &d_border, 6 * nodesSize*sizeof(float));
	if(!random)
	{
		cudaMemcpy(d_u_old, u_old, nodesSize*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(d_v_old, v_old, nodesSize*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(d_w_old, w_old, nodesSize*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(d_u_new, u_new, nodesSize*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(d_v_new, v_new, nodesSize*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(d_w_new, w_new, nodesSize*sizeof(double), cudaMemcpyHostToDevice);
	}
	cudaMemcpy(d_border, border, 6*nodesSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaBindTexture(0, borderTexture ,d_border,6*nodesSize*sizeof(float));

	if(!constantDiffusion)
	{
		cudaMalloc((void**) &d_z, 3 * nodesSize*sizeof(float));
		cudaMemcpy(d_z, z, 3 * nodesSize*sizeof(float), cudaMemcpyHostToDevice);
		cudaBindTexture(0, zTexture ,d_z,3*nodesSize*sizeof(float));
		if(correction)
		{
			cudaMalloc((void**) &d_Deff, 2 * nodesSize*sizeof(float));
			cudaMemcpy(d_Deff, Deff, 2 * nodesSize*sizeof(float), cudaMemcpyHostToDevice);
			cudaBindTexture(0, DeffTexture ,d_Deff,2*nodesSize*sizeof(float));			
		}
	}
	cout << "Initialization done" << endl;
	cudaDeviceSynchronize();


}

void cudaIteration3var_P(int number)
{
	if(number % 2 == 0)
	{
		//cout << "Starting old new" << endl;
		mainFunction3D3varBorder_P<<<DimGridSim, DimBlock>>>(d_u_old, d_u_new, d_v_old, d_v_new, d_w_old, d_w_new);
	}
	else
	{
		//cout << "Starting new old" << endl;
		mainFunction3D3varBorder_P<<<DimGridSim, DimBlock>>>(d_u_new, d_u_old, d_v_new, d_v_old, d_w_new, d_w_old);
	}
	cudaDeviceSynchronize();
	//gpuErrchk( cudaPeekAtLastError() );
}

void cudaIteration3var_P2D(int number)
{
	if(number % 2 == 0)
	{
		//cout << "Starting old new" << endl;
		mainFunction2D3varBorder_P<<<DimGridSim, DimBlock>>>(d_u_old, d_u_new, d_v_old, d_v_new, d_w_old, d_w_new);
	}
	else
	{
		//cout << "Starting new old" << endl;
		mainFunction2D3varBorder_P<<<DimGridSim, DimBlock>>>(d_u_new, d_u_old, d_v_new, d_v_old, d_w_new, d_w_old);
	}
	cudaDeviceSynchronize();
}

void cudaIteration3var_D2D(int number)
{
	if(number % 2 == 0)
	{
		//cout << "Starting old new" << endl;
		mainFunction2D3varBorder_D<<<DimGridSim, DimBlock>>>(d_u_old, d_u_new, d_v_old, d_v_new, d_w_old, d_w_new);
	}
	else
	{
		//cout << "Starting new old" << endl;
		mainFunction2D3varBorder_D<<<DimGridSim, DimBlock>>>(d_u_new, d_u_old, d_v_new, d_v_old, d_w_new, d_w_old);
	}
	cudaDeviceSynchronize();

	gpuErrchk( cudaPeekAtLastError() );
}

void cudaIteration3var_D2D_Correction(int number)
{
	if(number % 2 == 0)
	{
		//cout << "Starting old new" << endl;
		mainFunction2D3varBorder_D_Correction<<<DimGridSim, DimBlock>>>(d_u_old, d_u_new, d_v_old, d_v_new, d_w_old, d_w_new);
	}
	else
	{
		//cout << "Starting new old" << endl;
		mainFunction2D3varBorder_D_Correction<<<DimGridSim, DimBlock>>>(d_u_new, d_u_old, d_v_new, d_v_old, d_w_new, d_w_old);
	}
	cudaDeviceSynchronize();

	gpuErrchk( cudaPeekAtLastError() );
}

void copyToHost()
{
	if(iteration % 2 == 1)
	{
		//cout << "Copy new -> new" << endl;
		cudaMemcpy(u_new, d_u_new, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(v_new, d_v_new, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(w_new, d_w_new, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(u_old, d_u_old, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(v_old, d_v_old, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(w_old, d_w_old, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
	}
	else
	{
		//cout << "Copy old -> new" << endl;
		cudaMemcpy(u_new, d_u_old, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(v_new, d_v_old, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(w_new, d_w_old, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);	
		cudaMemcpy(u_old, d_u_new, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(v_old, d_v_new, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(w_old, d_w_new, nodesSize*sizeof(double), cudaMemcpyDeviceToHost);		
	}
}

void changeTimeStep(double e)
{
	//increase time step for 1.5
	if(e < 0.000005)
	{
		cout << "Increasing the time step" << endl;
		pars[1] = pars[1] * 1.5;
		cudaMemcpyToSymbol(parameters3var, pars, 20 * sizeof(double));
	}
	else if(e > 0.00009)
	{
		cout << "Decrease the time step" << endl;
		pars[1] = pars[1] * 0.5;
		cudaMemcpyToSymbol(parameters3var, pars, 20 * sizeof(double));
	}
	
}

void updateMeshSpacing(double newSpacing)
{
	h = newSpacing;
	pars[0] = newSpacing;
	cudaMemcpyToSymbol(parameters3var, pars, 20 * sizeof(double));
	cout << "Mesh spacing updated to: " << newSpacing << endl;
}

void updateSpaceDependantDiffusion()
{
	cudaMemcpy(d_z, z, 3 * nodesSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaBindTexture(0, zTexture ,d_z,3*nodesSize*sizeof(float));

	cudaMemcpy(d_Deff, Deff, 2 * nodesSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaBindTexture(0, DeffTexture ,d_Deff,2*nodesSize*sizeof(float));			

}

void freeMemory()
{
	cudaFree(d_u_old);
	cudaFree(d_u_new);
	cudaFree(d_v_old);
	cudaFree(d_v_new);
	cudaFree(d_w_old);
	cudaFree(d_w_new);
	cudaUnbindTexture(borderTexture);

	free(u_new);
	free(v_new);
	free(w_new);
}
