#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <mutex>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>

#include "loading.h"

pthread_barrier_t youCanCopyBarrier, youCanSaveBarrier;

//return true if simulation needs to be stopped
bool checkStoppingCriterion(double *e)
{
	double max = -1;
	for (int i = 0; i < a * b * c; ++i)
	{
		double e = abs(u_new[i] - u_old[i]) + abs(v_new[i] - v_old[i]) + abs(w_new[i] - w_old[i]);
		if (e > max && border[6 * i] > -2)
		{
			max = e;
		}
	}
	cout << iteration << ": " << max << endl;
	(*e) = max;
	if (max >= error)
	{
		return false;
	}
	return true;
}

//Outoput statistics at the end of the simulation
void getStats()
{
	double blackPixelsSum[3] = {0, 0, 0};
	double greenPixelsSum[3] = {0, 0, 0};
	int blackCount = 0;
	int greenCount = 0;
	for (int i = 0; i < nodesSize; ++i)
	{
		//color is black
		if (u_new[i] > v_new[i])
		{
			blackPixelsSum[0] += u_new[i];
			blackPixelsSum[1] += v_new[i];
			blackPixelsSum[2] += w_new[i];
			blackCount++;
		}
		else
		{
			greenPixelsSum[0] += u_new[i];
			greenPixelsSum[1] += v_new[i];
			greenPixelsSum[2] += w_new[i];
			greenCount++;
		}
	}
	cout << "Average concentrations black:" << endl;
	cout << blackPixelsSum[0] / blackCount << endl;
	cout << blackPixelsSum[1] / blackCount << endl;
	cout << blackPixelsSum[2] / blackCount << endl;
	cout << "Average concentrations green:" << endl;
	cout << greenPixelsSum[0] / greenCount << endl;
	cout << greenPixelsSum[1] / greenCount << endl;
	cout << greenPixelsSum[2] / greenCount << endl;
}

//save the average color in z direction as an image
void save_image_projection(char *fileName)
{
	char temp[200];
	strcpy(temp, fileName);
	strcat(temp, ".ppm");
	FILE *out = fopen(temp, "wb");
	fprintf(out, "P6 %d %d 255\n", a, b);
	for (int i = 0; i < b; ++i)
		for (int j = 0; j < a; ++j)
		{
			double sum_u = 0;
			double sum_v = 0;
			double height = 0;
			if (c > 1)
			{
				for (int k = 0; k < c; k++)
				{
					int index = k * (a * b) + j * b + i;
					if (border[6 * index] >= -1)
					{
						sum_u += u_new[index];
						sum_v += v_new[index];
						height += 1;
					}
				}
			}
			else
			{
				sum_u = u_new[j * b + i];
				sum_v = v_new[j * b + i];
				height = 1;
			}
			sum_u /= height;
			sum_v /= height;
			int G = ((sum_v - sum_u) / (sum_v + sum_u) + 0.5) * 255;
			if (G < 0)
				G = 0;
			if (G > 255)
				G = 255;
			putc(0, out);
			putc(G, out);
			putc(0, out);
		}
	fclose(out);
}

//save simulation on height k
void saveImage(char *fileName, int k)
{
	char temp[200];
	if (k >= c)
		k = c - 1;
	strcpy(temp, fileName);
	strcat(temp, ".ppm");
	FILE *out = fopen(temp, "wb");
	fprintf(out, "P6 %d %d 255\n", a, b);
	for (int i = 0; i < b; ++i)
		for (int j = 0; j < a; ++j)
		{
			int index = k * (a * b) + j * b + i;
			int R = (!(abs(border[6 * index] - 1 / h) > 0.0001) != !(abs(border[6 * index + 1] - 1 / h) > 0.0001)) && networkType != 4 ? 255 : 0;
			int G = 0;
			int B = 0;
			if (R == 0)
				R = (!(abs(border[6 * index + 2] - 1 / h) > 0.0001) != !(abs(border[6 * index + 3] - 1 / h) > 0.0001)) && networkType != 4 ? 255 : 0;
			if (border[6 * index] < -1)
			{
				R = 128;
				G = 128;
				B = 128;
			}
			else if (R == 0 && B == 0)
			{
				G = ((v_new[index] - u_new[index]) / (v_new[index] + u_new[index]) + 0.5) * 255;
				if (G < 0)
					G = 0;
				if (G > 255)
					G = 255;
			}
			putc(R, out);
			putc(G, out);
			putc(B, out);
		}
	fclose(out);
}

//save current simulation state
//if row is negative, save a s ply file
void save(char *saveFile, int row)
{
	if (row < 0)
	{
		char temp[200];
		strcpy(temp, saveFile);
		strcat(temp, ".ply");
		ofstream ply(temp);
		if (!ply)
		{
			cout << "Unable to create saving file " << saveFile << endl;
		}
		ply.precision(15);
		ply.setf(ios::fixed);
		ply.setf(ios::showpoint);
		ply << "ply \nformat ascii 1.0 \nelement vertex ";
		ply << nodesToSave;
		ply << "\nproperty double x \nproperty double y \nproperty double z";
		ply << "\nproperty uchar red  \nproperty uchar green \nproperty uchar blue\n";
		ply << "end_header\n";
		for (int i = 0; i < a; ++i)
			for (int j = 0; j < b; ++j)
				for (int k = 0; k < c; ++k)
				{
					int index = k * (a * b) + i * b + j;
					if (border[6 * index] < -1)
					{
						continue;
					}
					int value = ((v_new[index] - u_new[index]) / (v_new[index] + u_new[index]) + 0.5) * 255;
					if (value < 0)
						value = 0;
					if (value > 255)
						value = 255;
					ply << i * h;
					ply << " ";
					ply << j * h;
					ply << " ";
					ply << k * h;
					ply << " ";
					int R = 0;
					int G = 0;
					int B = 0;
					if (R == 0 && B == 0)
					{
						G = value;
					}
					ply << R;
					ply << " ";
					ply << G;
					ply << " ";
					ply << B;
					ply << endl;
				}
		ply.close();
	}
	else
	{
		saveImage(saveFile, row);
	}
}

//save current simulation state
void savePoints()
{
	bool first_iteration = true;
	while (!stop)
	{
		pthread_barrier_wait(&youCanSaveBarrier);
		cout << "Saving: " << saveIteration << endl;
		char saveFile[50];
		sprintf(saveFile, "%s/%d", saveFolder, saveIteration * saveFrequency);
		if (saveIteration == 0) //save in 3D
		{
			save(saveFile, -1);
		}
		double error;
		bool returnValue = checkStoppingCriterion(&error);
		cout << "Current: " << saveIteration * saveFrequency << endl;
		if(networkType == 0)
			save(saveFile, 0);
		else
			save(saveFile, c-1);
		if (returnValue && !first_iteration)
		{
			if (!growth)
			{
				stop = true;
				cout << "Saving thread barrier stop" << endl;
			}
			else
			{
				incrementE = true;
				first_iteration = true;
				cout << "Apply growth" << endl;
			}
			save(saveFile, -1);
		}
		else
		{
			first_iteration = false;
		}
		pthread_barrier_wait(&youCanCopyBarrier);
	}
}

void saveHexaStats()
{
	char saveFile[50];
	sprintf(saveFile, "%s/final_stats.txt", saveFolder);
	ofstream f(saveFile);
	for (int i = 0; i < centersSize; ++i)
	{
		for (int j = 0; j < 6; ++j)
		{
			f << neighbours[6 * i + j] << " ";
		}
		f << nodesPerHexa[i].size();
		for (int j = 0; j < nodesPerHexa[i].size(); ++j)
		{
			f << " " << u_new[nodesPerHexa[i][j]] << " " << v_new[nodesPerHexa[i][j]] << " " << w_new[nodesPerHexa[i][j]];
		}
		f << endl;
	}
	f.close();
}

void initHostArrays()
{
	u_new = (double *)malloc(nodesSize * sizeof(double));
	v_new = (double *)malloc(nodesSize * sizeof(double));
	w_new = (double *)malloc(nodesSize * sizeof(double));

	u_old = (double *)malloc(nodesSize * sizeof(double));
	v_old = (double *)malloc(nodesSize * sizeof(double));
	w_old = (double *)malloc(nodesSize * sizeof(double));
}

float epsilon_value(int step, float increment_factor, float initial_epsilon)
{
	return initial_epsilon + step * increment_factor;
}

int main(int argc, char **argv)
{
	if (argc < 5)
	{
		cout << "Argument 1: network type (1 - hexagonal prisms, 0 - gaussian bumps)" << endl;
		cout << "Argument 2: initial conditions if network type is 1 (1 - uniformly colored scales, 0 - randomized uniform steady state)" << endl;
		cout << "Argument 2: gaussian bumps sigma (integer) if network type is 0" << endl;
		cout << "Argument 3: total height of the simulation domain" << endl;
		cout << "Argument 4: domain thickness between prisms/gaussian bumps" << endl;
		cout << "Argument 5: mesh spacing epsilon" << endl;
		exit(-1);
	}
	h = atof(argv[5]);
	periodic = true;
	bool initialCond = atoi(argv[2]);
	border_thickness_in_elements = 1 / h;
	cout << "Element spacing: " << h << endl;
	cout << "Border thickness in elements: " << border_thickness_in_elements << endl;
	bt = 0.5 * border_thickness_in_elements + 1; //border thickess (for border thickness of x elements put 0.5*x + 1)
	dt = 0.012 * (h * h);
	error = 0.0000001;
	saveFrequency = 10000;
	constantDiffusion = true;
	//load network
	int res = 200;
	cout << "Resolution and spacing" << endl;
	cout << res << " " << h << endl;
	networkType = atoi(argv[1]);
	int maxZ = atoi(argv[3]);
	int sigma = -1;
	if(networkType == 0)
		sigma = atoi(argv[2]);
	int type = 1;
	if (networkType == 1) //hexagonal prisms
	{
		int freeZ = atoi(argv[4]);
		P = 0;
		c = maxZ;
		type = loadRegularNetwork("../networks/hexa10by10.txt", 1, 1, res);
		if (bt < 1.5)
		{
			getPixels(type, maxZ, freeZ); //one element at the border
			nodesToSave = nodesSize;
		}
		else
		{
			nodesToSave = getPixelsP0btMoreThanOne(type, maxZ, freeZ);
		}
		sprintf(saveFolder, "../output/prisms_%d_%d_%d_%s", initialCond, freeZ, maxZ, argv[5]);
		int result = mkdir(saveFolder, S_IRUSR | S_IWUSR | S_IXUSR);
		nodesSize = a * b * c;
		cout << "Nodes size" << nodesSize << endl;
	}
	else if (networkType == 0)  //gaussian bumps
	{
		if(sigma < 0)
			cout << "Wrong sigma value" << endl;
		int freeZ = atoi(argv[4]);
		char gaussNetworkFile[200];
		sprintf(gaussNetworkFile, "../networks/gauss10by10_%d_%d_%d.txt", freeZ, maxZ, sigma);
		P = 0;
		periodic = false;
		nodesToSave = getPixelsGaussianBumps(gaussNetworkFile);
		nodesSize = a * b * c;
		sprintf(saveFolder, "../output/gauss_%d_%d_%d_%s", freeZ, maxZ, sigma, argv[5]);
		int result = mkdir(saveFolder, S_IRUSR | S_IWUSR | S_IXUSR);
	}
	else
	{
		cout << "Wrong network type" << endl;
		exit(-1);
	}
	//just in case
	correctBordersToTheNonsimulatedNodes();
	if (type == 1 && periodic)
		hexaNetworkCorrection();

	cout << "Dimensions:" << a << " " << b << " " << c << " " << endl;
	cout << "Number of nodes: " << nodesSize << endl;
	cout << "Number of nodes to save: " << nodesToSave << endl;

	//simulation parameters
	c1 = -0.04;
	c2 = -0.056;
	c3 = 0.382;
	c4 = -0.05;
	c5 = 0;
	c6 = 0.25;
	c7 = 0.016;
	c8 = -0.03;
	c9 = 0.24;
	cu = 0.02;
	cV = 0.025;
	cw = 0.06;

	Du = 1.125;
	Dv = 1.125;
	Dw = 12 * Du;

	cout << "Diffusion coeff " << Du << endl;

	U = 0.5;
	V = 0.5;
	W = 0.5;

	if(initialCond && networkType == 1)
	{
		initialGreen[0] = 1.2;
		initialGreen[1] = 6.6;
		initialGreen[2] = 2.3;
		initialBlack[0] = 5.3;
		initialBlack[1] = 0.92;
		initialBlack[2] = 4;
	}
	else
	{
		initialUniform[0] = 3.47561;
		initialUniform[1] = 3.0488;
		initialUniform[2] = 3.40243;
		initialUniform[0] = 3.47561;
		initialUniform[1] = 3.0488;
		initialUniform[2] = 3.40243;
	}

	double randomIntervalSize = 0.1;
	initHostArrays();
		//can be replaced by loadInitialConcentrations to continue a simulation
		//bool ret = loadInitialConcentrations("/home/lane/PROJECTS/RD/code/concentrations_7.txt");
	if(networkType == 1)
		bool ret = loadInitialColors("../networks/hexa10by10_c.txt", randomIntervalSize, type, -1);
	else
		createRandomInitialConditions(randomIntervalSize);

	cout << "Initializing cuda" << endl;
	initCuda(false);
	iteration = 0;
	stop = false;
	saving = true;

	pthread_barrier_init(&youCanCopyBarrier, NULL, 2);
	pthread_barrier_init(&youCanSaveBarrier, NULL, 2);
	thread savingThread(savePoints);
	copyToHost();
	saveIteration = 0;
	pthread_barrier_wait(&youCanSaveBarrier);

	cout << "Starting the simulation" << endl;

	if (c == 1)
	{
		cout << "Running 2D version" << endl;
		while (1)
		{
			cout << "Thread 0: " << iteration << endl;
			for (int i = 0; i < saveFrequency; ++i)
			{
				cudaIteration3var_P2D(iteration);
				iteration++;
			}
			pthread_barrier_wait(&youCanCopyBarrier);
			if (stop)
			{
				cout << "Done" << endl;
				break;
			}
			copyToHost();
			saveIteration++;
			pthread_barrier_wait(&youCanSaveBarrier);
		}
	}
	else
	{
		cout << "Running 3D version" << endl;
		while (1)
		{
			cout << "Thread 0: " << iteration << endl;
			for (int i = 0; i < saveFrequency; ++i)
			{
				cudaIteration3var_P(iteration);
				iteration++;
			}
			pthread_barrier_wait(&youCanCopyBarrier);
			if (stop)
			{
				cout << "Done" << endl;
				break;
			}
			copyToHost();
			saveIteration++;
			pthread_barrier_wait(&youCanSaveBarrier);
		}
	}
	savingThread.join();
	return 0;
}