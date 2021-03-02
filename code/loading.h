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


#include "kernel.cu"


#define PI 3.14159265
#define sin60 0.866025
#define sin30 0.5

/*
-----used only for 2D simulations-----
2*i contains the diffusion coefficient
2*i + 1 contains its derivative in x direction
2*i + 2 contains its derivative in y direction
*/

//meshes have the same number of elements
//we are saving just the height


int loadDiffusionFromPly(const char *top_mesh, const char *bottom_mesh)
{
	z = (float*)malloc(3*a*b*sizeof(float));
	bool equal = false;
	if(strcmp(top_mesh, bottom_mesh) == 0) equal = true;
	for(int i = 0; i < a * b; ++i)
	{
		z[3 * i] = h;
	}
	ifstream f(top_mesh);
	ifstream bo(bottom_mesh);
	if(!f || !bo)
	{
		cout << "Diffusion ply file doesn't exist" << endl;
		return -1;
	}
	string line;
	for(int i = 0; i < 9; ++i) getline(f,line);
	for(int i = 0; i < 9; ++i) getline(bo,line);
	float temp, temp_b;
	float min = 1000;
	float max = -1000;
	int a1, b1;
	bool zero_detected = false;
	for(int i = 0; i < a * b; ++i)
	{
		f >> temp;
		bo >> temp_b;
		a1 = int(std::round(temp));
		if(temp != temp_b)
		{
			cout << "The files are not sincronized" << endl;
			return -1;
		}
		f >> temp;
		bo >> temp_b;
		if(temp != temp_b)
		{
			cout << "The files are not sincronized" << endl;
			return -1;
		}
		b1 = int(std::round(temp));
		f >> temp;
		bo >> temp_b;
		//cout << temp << " " << temp_b << endl;
		
		if(equal) z[3*(a1 * b + b1)] = temp * h;
		else z[3*(a1 * b + b1)] = (temp - temp_b)* h;
		if(z[3*(a1 * b + b1)] == 0) zero_detected = true;
		if(z[3*(a1 * b + b1)] < min) min = z[3*(a1 * b + b1)];
		if(z[3*(a1 * b + b1)] > max)
		{
			max = z[3*(a1 * b + b1)];
		} 
	}
	cout << "Minimum height: " << min << endl;
	cout << "Maximum height: " <<  max << endl;
	//update all the derivatives and save file
	char temp1[200];
	sprintf(temp1, "%s/diffusion_%f_%f.ppm", saveFolder, min , max);
	FILE * out = fopen (temp1, "wb");
	fprintf(out, "P6 %d %d 255\n", a, b);
	for(int i = 0; i < b; ++i)
		for(int j = 0; j < a; ++j) 
 		{
			int index = j * b + i;
			float value = (z[3*index] - min)/(max - min);
			if(zero_detected) z[3*index] += h/2;
 			int R = int(value * 255);
 			int G = 0;
 			int B = 0;
 			
 			putc(R, out);
 			putc(G, out);
 			putc(B, out);
		}
	fclose(out);
	//calculate the derivatives in the x and y direction
	max = -1000;
	for(int i = 0; i < a; ++i)
		for(int j = 0; j < b; ++j)
		{
			int index = i * b + j;
			int xb = (i - 1) * b + j;
			int xa = (i + 1) * b + j;
			int yb = i * b + j - 1;
			int ya = i * b + j + 1;
			//x direction
			if(i == 0)
			{
				z[3 * index + 1] = (z[3 * xa] - z[3 * index])/h;
			}
			else if(i == a - 1)
			{
				z[3 * index + 1] = (z[3 * index] - z[3 * xb])/h;		
			}
			else
			{
				z[3 * index + 1] = (z[3 * xa] - z[3 * xb])/(2 * h);				
			}
			//y direction
			if(j == 0)
			{
				z[3 * index + 2] = (z[3 * ya] - z[3 * index])/h;
			}
			else if(j == b - 1)
			{
				z[3 * index + 2] = (z[3 * index] - z[3 * yb])/h;		
			}
			else
			{
				z[3 * index + 2] = (z[3 * ya] - z[3 * yb])/(2 * h);				
			}
			if(abs(z[3 * index + 2]) > max) max = abs(z[3 * index + 2]);
			if(abs(z[3 * index + 1]) > max) max = abs(z[3 * index + 1]);
		}	
	cout << "Maximum der norm: " <<  max << endl;
	sprintf(temp1, "%s/%f_diffusion_derivatives.ppm", saveFolder,max);
	FILE * out1 = fopen (temp1, "wb");
	fprintf(out1, "P6 %d %d 255\n", a, b);
	for(int i = 0; i < b; ++i)
		for(int j = 0; j < a; ++j) 
 		{
			int index = j * b + i;
			float value = (z[3*index + 2] - min)/(max - min);
		
 			int R = int(value * 255);
 			int G = 0;
 			int B = 0;
 			
 			putc(R, out1);
 			putc(G, out1);
 			putc(B, out1);
		}
	fclose(out1);
	cout << "Done" << endl;
	return 1;

}

int loadDiffusionFromPlyCorrection(const char *top_mesh, const char *bottom_mesh, bool first)
{
	if(first)
	{
		z = (float*)malloc(3*a*b*sizeof(float));	
		Deff = (float*)malloc(2*a*b*sizeof(float));
	}
	float* midline = (float*)malloc(3*a*b*sizeof(float));
	bool equal = false;
	if(strcmp(top_mesh, bottom_mesh) == 0)
	{
		equal = true;
		cout << "The meshes are equal" << endl;
	}
	else
	{
		cout << "The meshes are not equal" << endl;
	}
	for(int i = 0; i < a * b; ++i)
	{
		z[3 * i] = h;
	}
	ifstream f(top_mesh);
	ifstream bo(bottom_mesh);
	if(!f || !bo)
	{
		cout << "Diffusion ply file doesn't exist" << endl;
		return -1;
	}
	string line;
	for(int i = 0; i < 9; ++i) getline(f,line);
	for(int i = 0; i < 9; ++i) getline(bo,line);
	float temp, temp_b;
	float min = 1000;
	float max = -1000;
	int a1, b1;
	bool zero_detected = false;
	for(int i = 0; i < a * b; ++i)
	{
		f >> temp;
		bo >> temp_b;
		a1 = int(std::round(temp));
		if(temp != temp_b)
		{
			cout << "The files are not sincronized" << endl;
			return -1;
		}
		f >> temp;
		bo >> temp_b;
		if(temp != temp_b)
		{
			cout << "The files are not sincronized" << endl;
			return -1;
		}
		b1 = int(std::round(temp));
		f >> temp;
		bo >> temp_b;
		if(equal)
		{
			z[3*(a1 * b + b1)] = temp * h;
			midline[3 * (a1 * b + b1)] = temp/2 * h; //because we just need the derivatives	
		}
		else
		{
			z[3*(a1 * b + b1)] = (temp - temp_b) * h;
			midline[3 * (a1 * b + b1)] = (temp + temp_b)/2 * h; //because we just need the derivatives
		}
		if(z[3*(a1 * b + b1)] == 0)
		{
			zero_detected = true;
			cout << "Zero detected" << endl;
		}	
		if(z[3*(a1 * b + b1)] < min) min = z[3*(a1 * b + b1)];
		if(z[3*(a1 * b + b1)] > max)
		{
			max = z[3*(a1 * b + b1)];
		} 
	}	
	cout << "Minimum height: " << min << endl;
	cout << "Maximum height: " <<  max << endl;
	//update all the derivatives and save file
	char temp1[200];
	sprintf(temp1, "%s/diffusion_%f_%f.ppm", saveFolder, max, min);
	FILE * out = fopen (temp1, "wb");
	fprintf(out, "P6 %d %d 255\n", a, b);
	for(int i = 0; i < b; ++i)
		for(int j = 0; j < a; ++j) 
 		{
			int index = j * b + i;
			float value = (z[3*index] - min)/(max - min);
			if(zero_detected)
			{
				z[3*index] += h/2;
				midline[3*index] += h/4;
			}
 			int R = int(value * 255);
 			int G = 0;
 			int B = 0;
 			
 			putc(R, out);
 			putc(G, out);
 			putc(B, out);
		}
	fclose(out);
	//calculate the derivatives in the x and y direction
	max = -1000;
	for(int i = 0; i < a; ++i)
		for(int j = 0; j < b; ++j)
		{
			int index = i * b + j;
			int xb = (i - 1) * b + j;
			int xa = (i + 1) * b + j;
			int yb = i * b + j - 1;
			int ya = i * b + j + 1;
			//x direction
			if(i == 0)
			{
				z[3 * index + 1] = (z[3 * xa] - z[3 * index])/(h);
				midline[3 * index + 1] = (midline[3 * xa] - midline[3 * index])/(h);
			}
			else if(i == a - 1)
			{
				z[3 * index + 1] = (z[3 * index] - z[3 * xb])/(h);
				midline[3 * index + 1] = (midline[3 * index] - midline[3 * xb])/(h);			
			}
			else
			{
				z[3 * index + 1] = (z[3 * xa] - z[3 * xb])/(2 * h);	
				midline[3 * index + 1] = (midline[3 * xa] - midline[3 * xb])/(2 * h);				
			}
			//y direction
			if(j == 0)
			{
				z[3 * index + 2] = (z[3 * ya] - z[3 * index])/(h);
				midline[3 * index + 2] = (midline[3 * ya] - midline[3 * index])/(h);
			}
			else if(j == b - 1) 
			{
				z[3 * index + 2] = (z[3 * index] - z[3 * yb])/(h);
				midline[3 * index + 2] = (midline[3 * index] - midline[3 * yb])/(h);		
			}
			else
			{
				z[3 * index + 2] = (z[3 * ya] - z[3 * yb])/(2 * h);	
				midline[3 * index + 2] = (midline[3 * ya] - midline[3 * yb])/(2 * h);			
			}

			Deff[2 * index] = 1.0/(1 + midline[3 * index + 1]*midline[3 * index + 1] + 1.0/12.0 * z[3 * index + 1]*z[3 * index + 1]);
			Deff[2 * index + 1] = 1.0/(1 + midline[3 * index + 2]*midline[3 * index + 2] + 1.0/12.0 * z[3 * index + 2]*z[3 * index + 2]);				
			if(index == 500)
			{
			cout << 1.0/(1 + midline[3 * index + 2]*midline[3 * index + 2] + 1.0/12.0 * z[3 * index + 2]*z[3 * index + 2]) << endl;
			cout << 1.0/(1 + 1.0/3.0 * z[3 * index + 2]*z[3 * index + 2]) << endl;
			}
			if(abs(z[3 * index + 1]) > max) max = abs(z[3 * index + 1]);
			if(abs(z[3 * index + 2]) > max) max = abs(z[3 * index + 2]);
		}	
	for(int i = 0; i < a; ++i)
		for(int j = 0; j < b; ++j)
		{
			int index = i * b + j;
			int xb = (i - 1) * b + j;
			int xa = (i + 1) * b + j;
			int yb = i * b + j - 1;
			int ya = i * b + j + 1;
			float xDeff = 0;
			float yDeff = 0;
			//x direction
			if(i == 0)
			{
				xDeff = (Deff[2 * xa] - Deff[2 * index])/h;
			}
			else if(i == a - 1)
			{
				xDeff = (Deff[2 * index] - Deff[2 * xb])/h;		
			}
			else
			{
				xDeff = (Deff[2 * xa] - Deff[2 * xb])/(2 * h);				
			}
			//y direction
			if(j == 0)
			{
				yDeff = (Deff[2 * ya + 1] - Deff[2 * index + 1])/h;
			}
			else if(j == b - 1)
			{
				yDeff = (Deff[2 * index + 1] - Deff[2 * yb + 1])/h;		
			}
			else
			{
				yDeff = (Deff[2 * ya + 1] - Deff[2 * yb + 1])/(2 * h);				
			}		
			z[3 * index + 1] = z[3 * index + 1] * Deff[2 * index] + z[3 * index] * xDeff;
			z[3 * index + 2] = z[3 * index + 2] * Deff[2 * index + 1] + z[3 * index] * yDeff;
		}	
	
	cout << "Maximum grad norm: " <<  max << endl;
	sprintf(temp1, "%s/%f_diffusion_derivatives.ppm", saveFolder,max);
	FILE * out1 = fopen (temp1, "wb");
	fprintf(out1, "P6 %d %d 255\n", a, b);
	for(int i = 0; i < b; ++i)
		for(int j = 0; j < a; ++j) 
 		{
			int index = j * b + i;
			float value = (z[3*index + 2] - min)/(max - min);
		
 			int R = int(value * 255);
 			int G = 0;
 			int B = 0;
 			
 			putc(R, out1);
 			putc(G, out1);
 			putc(B, out1);
		}
	fclose(out1);
	free(midline);
	cout << "Done" << endl;
	return 1;

}

int returnClosestCenter(int x, int y)
{
	double mini = 10000000000;
	int minIndex = -1;
	for(int i = 0 ; i < centersSize; ++i)
	{
		double distance = ((x - centers[2 * i])/sx)*((x - centers[2 * i])/sx) + ((y - centers[2 * i + 1])/sy)*((y - centers[2 * i + 1])/sy);
		if(distance < mini)
		{
			mini = distance;
			minIndex = i;
		}
	}
	return minIndex;
}

void correct_D_and_z(float old_h, float new_h)
{
	float* midline = (float*)malloc(3*a*b*sizeof(float));
	for(int i = 0; i < a * b; ++i)
	{
		z[3 * i] *= new_h/old_h; 
		midline[3 * i] = z[3 * i]/2;
	}
	for(int i = 0; i < a; ++i)
		for(int j = 0; j < b; ++j)
		{
			int index = i * b + j;
			int xb = (i - 1) * b + j;
			int xa = (i + 1) * b + j;
			int yb = i * b + j - 1;
			int ya = i * b + j + 1;
			//x direction
			if(i == 0)
			{
				z[3 * index + 1] = (z[3 * xa] - z[3 * index])/(h);
				midline[3 * index + 1] = (midline[3 * xa] - midline[3 * index])/(h);
			}
			else if(i == a - 1)
			{
				z[3 * index + 1] = (z[3 * index] - z[3 * xb])/(h);
				midline[3 * index + 1] = (midline[3 * index] - midline[3 * xb])/(h);			
			}
			else
			{
				z[3 * index + 1] = (z[3 * xa] - z[3 * xb])/(2 * h);	
				midline[3 * index + 1] = (midline[3 * xa] - midline[3 * xb])/(2 * h);				
			}
			//y direction
			if(j == 0)
			{
				z[3 * index + 2] = (z[3 * ya] - z[3 * index])/(h);
				midline[3 * index + 2] = (midline[3 * ya] - midline[3 * index])/(h);
			}
			else if(j == b - 1) 
			{
				z[3 * index + 2] = (z[3 * index] - z[3 * yb])/(h);
				midline[3 * index + 2] = (midline[3 * index] - midline[3 * yb])/(h);		
			}
			else
			{
				z[3 * index + 2] = (z[3 * ya] - z[3 * yb])/(2 * h);	
				midline[3 * index + 2] = (midline[3 * ya] - midline[3 * yb])/(2 * h);			
			}

			Deff[2 * index] = 1.0/(1 + midline[3 * index + 1]*midline[3 * index + 1] + 1.0/12.0 * z[3 * index + 1]*z[3 * index + 1]);
			Deff[2 * index + 1] = 1.0/(1 + midline[3 * index + 2]*midline[3 * index + 2] + 1.0/12.0 * z[3 * index + 2]*z[3 * index + 2]);				
		}	
		for(int i = 0; i < a; ++i)
			for(int j = 0; j < b; ++j)
			{
				int index = i * b + j;
				int xb = (i - 1) * b + j;
				int xa = (i + 1) * b + j;
				int yb = i * b + j - 1;
				int ya = i * b + j + 1;
				float xDeff = 0;
				float yDeff = 0;
				//x direction
				if(i == 0)
				{
					xDeff = (Deff[2 * xa] - Deff[2 * index])/h;
				}
				else if(i == a - 1)
				{
					xDeff = (Deff[2 * index] - Deff[2 * xb])/h;		
				}
				else
				{
					xDeff = (Deff[2 * xa] - Deff[2 * xb])/(2 * h);				
				}
				//y direction
				if(j == 0)
				{
					yDeff = (Deff[2 * ya + 1] - Deff[2 * index + 1])/h;
				}
				else if(j == b - 1)
				{
					yDeff = (Deff[2 * index + 1] - Deff[2 * yb + 1])/h;		
				}
				else
				{
					yDeff = (Deff[2 * ya + 1] - Deff[2 * yb + 1])/(2 * h);				
				}		
				z[3 * index + 1] = z[3 * index + 1] * Deff[2 * index] + z[3 * index] * xDeff;
				z[3 * index + 2] = z[3 * index + 2] * Deff[2 * index + 1] + z[3 * index] * yDeff;
			}	
	free(midline);
}

void correctBordersToTheNonsimulatedNodes()
{
	float *newBorder = (float*)malloc(a*b*c*6*sizeof(float));
	memcpy(newBorder, border, a*b*c*6*sizeof(float));
	for(int i = 0;i < a*b*c; ++i)
	{
			if(border[6*i] == -2)
			{
				int z = i/(a * b);
				int temp = i - z * (a * b);
				int x = (temp / b);
				int y = (temp % b);	
				//left
				int ind = (x > 0) ? z * (a * b) + (x - 1) * b + y : z * (a * b) + (a - 1) * b + y;
				if(border[6 * ind + 1] > 0) newBorder[6 * ind + 1] = 0;
				//right
				ind = (x < a-1) ? z * (a * b) + (x + 1) * b + y : z * (a * b) + y;
				if(border[6 * ind] > 0) newBorder[6 * ind] = 0;
				//bottom
				ind = (y > 0) ? z * (a * b) + x * b + y - 1 : z * (a * b) + x * b + b-1;
				if(border[6 * ind + 3] > 0) newBorder[6 * ind + 3] = 0;
				//top
				ind = (y < b-1) ? z * (a * b) + x * b + y + 1 : z * (a * b) + x * b;
				if(border[6 * ind + 2] > 0) newBorder[6 * ind + 2] = 0;
				//up
				if(z < c - 1)
				{
					ind = (z+1) * (a * b) + x * b + y;
					if(border[6 * ind + 4] > 0) newBorder[6 * ind + 4] = 0; 
				}
				if(z > 0)
				{
					ind = (z-1) * (a * b) + x * b + y;
					if(border[6 * ind + 5] > 0) newBorder[6 * ind + 5] = 0;					
				}				
			}
	}
	 memcpy(border, newBorder, a*b*c*6*sizeof(float));
}

void hexaNetworkCorrection()
{
	for(int i = 0; i < centersSize; ++i)
	{
		for(int j = 0; j < contributions[i].size(); ++j)
		{
			for(int k = 0; k < nodesPerHexa[contributions[i][j]].size(); ++k)
			{
				nodesPerHexa[i].push_back(nodesPerHexa[contributions[i][j]][k]); 
				closest_el[nodesPerHexa[contributions[i][j]][k]] = i;
			}
			nodesPerHexa[contributions[i][j]].clear();
		}
	}
}


int loadRegularNetwork(const char *fileName, double scalex, double scaley, int xdim)
{
	ifstream f(fileName);
	if(!f)
	{
		cout << "Network file doesn't exist" << endl;
		return -1;
	}
	int type, size, index, contrNum, el, temp;
	f >> type;
	f >> sx;
	f >> sy;
	f >> size;
	centers = (double*)malloc(2 * size * sizeof(double));
	double maxX = -1;
	double maxY = -1;
	centersSize = size;
	cout << type << endl;
	//load rectangles
	if(type == 0)
	{
		for(int i = 0; i < size; ++i)
		{
			f >> index;
			f >> centers[2 * i];
			centers[2 * i] += 0.5;
			f >> centers[2 * i + 1];
			centers[2 * i + 1] += 0.5;
			f >> temp;
			if(centers[2 * i] > maxX) maxX = centers[2 * i];
			if(centers[2 * i + 1] > maxY) maxY = centers[2 * i + 1];
			for(int j = 0 ; j < 8; ++j)
			{
				f >> temp; 
			}
			f >> contrNum;
			for(int j = 0; j < contrNum; ++j)
			{
				f >> el;
			}
		}
		maxX -= 0.5;
		maxY -= 0.5;
		//scale all centers
		double scale_factor = double(xdim)/(maxX  * scalex);
		cout << scale_factor << endl;
		// distance between neighboursbouring centers
		for(int i = 0; i < size; ++i)
		{
			centers[2 * i] *= scale_factor * scalex - 0.00001;
			centers[2 * i + 1] *= scale_factor * scaley - 0.00001;
		}
		a = int(maxX * scale_factor * scalex);
		b = int(maxY * scale_factor * scaley);
		sx = sx * scalex * scale_factor;
		sy = sy * scaley * scale_factor;
	}
	//load hexagons
	else if(type == 1)
	{
		neighbours = (int*)malloc(6 * size * sizeof(int));
		map<int, int> indexCorr;
		for(int i = 0; i < size; ++i)
		{
			f >> index;
			indexCorr[index] = i;
			f >> centers[2 * i] >> centers[2 * i + 1];
			f >> temp;
			if(centers[2 * i] > maxX) maxX = centers[2 * i];
			if(centers[2 * i + 1] > maxY) maxY = centers[2 * i + 1];
			for(int j = 0 ; j < 6; ++j)
			{
				f >> neighbours[6 * i + j]; 
			}
			f >> contrNum;
			for(int j = 0; j < contrNum; ++j)
			{
				f >> el;
				contributions[i].push_back(el);
			}
		}
		for(int i = 0; i < 6 * size; ++i)
		{
			if(neighbours[i] != -1)
			{
				neighbours[i] = indexCorr[neighbours[i]];
			}
		}
		for(int i = 0; i < size; ++i)
		{
			for(int j = 0; j < contributions[i].size(); ++j)
			{
				contributions[i][j] = indexCorr[contributions[i][j]];
			}
		}
		//scale all centers
		if(sx == 1)
		{
			double scale_factor = double(xdim)/(maxX  * scalex);
			cout << scale_factor << endl;
			// distance between neighboursbouring centers
			for(int i = 0; i < size; ++i)
			{
				centers[2 * i] *= scale_factor * scalex - 0.00001;
				centers[2 * i + 1] *= scale_factor * scaley + 0.00001;
			}
			a = int(maxX * scale_factor * scalex);
			b = int(maxY * scale_factor * scaley);
			sx = sx * scalex * scale_factor;
			sy = sy * scaley * scale_factor;
		}
	}
	cout << a << " " << b << endl;
	return type;
}

//putting the border green
void randomizeInitialColorsUsingPercentage(float percBlack)
{
	cout << "Randomizing initial colors using percentage" << endl;
	srand (time(0));
	int totalNumber = 0;
	int black = 0;
	bool black_border = false;
	//if((double)rand()/(RAND_MAX) < 0.5) black_border = true;
	for(int i = 0; i < centersSize; ++i)
	{
		if((double)rand()/(RAND_MAX) < percBlack)
		{
			initialColors[i] = true;
			black++;
		}
		else
		{
			initialColors[i] = false;
		}
		if(neighbours[6*i] == -1)
		{
			if(black_border) initialColors[i] = true;
			else initialColors[i] = false;
		}
		totalNumber++;
	}
	cout << "Percentage of black scales: " << float(black)/totalNumber << endl; 
}
//closest has to be specified
bool loadInitialColors(const char *fileName, double randomIntervalSize, int type, float blackPercentage)
{
	ifstream f(fileName);
	if(!f)
	{
		cout << "Colors file doesn't exist" << endl;
		return false;
	}
	int size;
	bool temp;
	f >> size;
	initialColors = (bool*)malloc(size * sizeof(bool));
	int goodColors = 0;
	for(int i = 0; i < size; ++i)
	{
		if(centers[2 * i] >= 0 && centers[2 * i] < a && centers[2 * i + 1] >= 0 && centers[2 * i + 1] < b)
		{
			if(type == 0)
			{
				f >> initialColors[i];
			}
			else
			{
				
				f >> initialColors[i];
				f >> temp;
			}
			goodColors++;
		}
	}
	srand (time(0));
	cout << "Good colors loaded " << goodColors << endl;
	if(blackPercentage > 0) randomizeInitialColorsUsingPercentage(blackPercentage);
	for(int i = 0; i < a * b * c; ++i)
	{
		//true is black
		if(closest_el[i] < 0 || border[6 * i] < -1 || neighbours[6*closest_el[i]] == -1)
		{
			u_old[i] = initialUniform[0] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
			v_old[i] = initialUniform[1] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;;
			w_old[i] = initialUniform[2] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;;
			continue;
		}		
		if(initialColors[closest_el[i]])
		{
			u_old[i] = initialBlack[0] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
			v_old[i] = initialBlack[1] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
			w_old[i] = initialBlack[2] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
		}
		else
		{
			u_old[i] = initialGreen[0] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
			v_old[i] = initialGreen[1] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
			w_old[i] = initialGreen[2] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;			
		}			
	}
	return true;
}

bool loadInitialConcentrations(const char *fileName)
{
	cout << "Reading concentration file" << endl;
	ifstream f(fileName);
	if(!f)
	{
		cout << "Concentration file doesn't exist" << endl;
		return false;
	}
	int size;
	f >> size;
	int i,j,k;
	for(int m = 0; m < size; ++m)
	{
		f >> i >> j >> k;
		int index = k * (a * b) + i * b + j;
		f >> u_old[index] >> v_old[index] >> w_old[index];
	}
	return true;
}

bool checkIfHole(map<pair<int, int>, bool> *holes, int x, int y)
{
	if(holes->count(make_pair(x,y)) == 1) return true;
	return false;
}


int getPixelsUsingUniformHoles(int polygon, const char *fileName)
{
	border = (float*)malloc(6 * a * b * c * sizeof(float));
	closest_el = (int*)malloc(a * b * c * sizeof(int));
	int numberOfNodesToSimulate = 0;
	for(int i = 0 ; i < 6 * a * b * c; ++i)
	{
		border[i] = 1/h;
	}
	//load uniform holes
	map<pair<int, int>, bool> holes; 
	ifstream f(fileName);
	if(!f)
	{
		cout << "Holes file doesn't exist" << endl;
		return -1;
	}
	int size, hx, hy;
	f >> size;
	for(int i = 0; i < size; ++i)
	{
		f >> hx;
		f >> hy;
		holes[make_pair(hx, hy)] = true;
	}
	//done loading holes
	if (polygon == 0)
	{	
		double maxDistance = float(border_thickness_in_elements)/2;
		bool dont_simulate = true;
		if(maxDistance == 0)
		{
			dont_simulate = false;
			maxDistance = 0.5;
		}
		for(int i = 0; i < a * b * c; ++i)
		{
			int z = i/(a * b);
			int temp = i - z * (a * b);
			int x = (temp / b);
			int y = (temp % b);	
			int closest = returnClosestCenter(x,y);
			closest_el[i] = closest;
			//left side
			double d = abs(x - centers[2 * closest] + sx/2);
			if(d <= maxDistance) 
			{
				if(dont_simulate && !checkIfHole(&holes, y, z))
				{
					border[6 * i] = -2;
				}
				else if(!dont_simulate && !checkIfHole(&holes, y, z))
				{
					border[6 * i] = 0;
				}
				
			}
			//right side
			d = abs(x - centers[2 * closest] - sx/2);
			if(d <= maxDistance) 
			{
				if(dont_simulate && !checkIfHole(&holes, y, z))
				{
					border[6 * i] = -2;
				}
				else if(!dont_simulate && !checkIfHole(&holes, y, z))
				{
					border[6 * i + 1] = 0;
				}
				
			}	
			//top
			d = abs(y - centers[2 * closest + 1] - sy/2);
			if(d <= maxDistance) 
			{
				if(dont_simulate && !checkIfHole(&holes, x, z))
				{
					border[6 * i] = -2;
				}
				else if(!dont_simulate && !checkIfHole(&holes, x, z))
				{
					border[6 * i + 3] = 0;
				}
				
			}	
			//bottom
			d = abs(y - centers[2 * closest + 1] + sy/2);
			if(d <= maxDistance) 
			{
				if(dont_simulate && !checkIfHole(&holes, x, z))
				{
					border[6 * i] = -2;
				}
				else if(!dont_simulate && !checkIfHole(&holes, x, z))
				{
					border[6 * i + 2] = 0;
				}
				
			}	
			if(z == 0)
			{
				border[6 * i + 4] = 0;
			}
			if(z == c - 1)
			{
				border[6 * i + 5] = 0;
			}
			if(x < maxDistance * 1.1 && !periodic)
			{
				if(dont_simulate) border[6 * i] = -2;
				else
				{
					border[6 * i] = 0;
				}
			}	
			if(y < maxDistance * 1.1 && !periodic)
			{
				if(dont_simulate) border[6 * i] = -2;
				else
				{
					border[6 * i + 2] = 0;
				}
			}		
			if(abs(x - a + 1) <  maxDistance * 1.1 && !periodic)
			{
				if(dont_simulate) border[6 * i] = -2;
				else
				{
					border[6 * i + 1] = 0;
				}
			}	
			if(abs(y - b + 1) < maxDistance * 1.1 && !periodic)
			{
				if(dont_simulate) border[6 * i] = -2;
				else
				{
					border[6 * i + 3] = 0;
				}
			}	
			if(border[6 * i] != -2)
			{
				numberOfNodesToSimulate++;	
				nodesPerHexa[closest].push_back(i);
			}					
		}
	}
	return numberOfNodesToSimulate;	
} 


int getPixelsGaussianBumps(const char *file)
{
	ifstream f(file);
	int numberOfNodesTosimulate = 0;
	if(!f)
	{
		cout << "Bumps file doesn't exist" << endl;
		return -1;
	}
	f >> a >> b >> c;
	border = (float*)malloc(6 * a * b * c * sizeof(float));
	closest_el = (int*)malloc(a * b * c * sizeof(int));
	for(int i = 0 ; i < 6 * a * b * c; ++i)
	{
		border[i] = 1/h;
	}
	int x,y,z, temp;	
	for(int k = 0; k < a * b * c; ++k)
	{
		f >> x >> y >> z;
		int i = z * (a * b) + x * b + y;
		for(int j = 0; j < 6; ++j)
		{
			f >> temp;
			if(temp == 1) border[6 * i + j] = 0;
			else if(temp == -2) border[6 * i + j] = -2;
			//else border[6 * i + j] = temp;
		}
		if(border[6 * i] < -1) 
			closest_el[i] = -1;
		else
		{
			numberOfNodesTosimulate++;
			closest_el[i] = returnClosestCenter(x,y);
		}
			//use to fix a specific elements
	if(closest_el[i] == 1000000)
	{
		border[6 * i] = -1;
	}
	}
	
	return numberOfNodesTosimulate;
} 


void getPixels(int polygon, int maxZ, int lowZ)
{
	border = (float*)malloc(6 * a * b * c * sizeof(float));
	closest_el = (int*)malloc(a * b * c * sizeof(int));
	for(int i = 0 ; i < 6 * a * b * c; ++i)
	{
		border[i] = 1/h;
	}		
	if (polygon == 0)
	{	
		for(int i = 0; i < a * b * c; ++i)
		{
			int z = i/(a * b);
			int temp = i - z * (a * b);
			int x = (temp / b);
			int y = (temp % b);	
			int closest = returnClosestCenter(x,y);
			closest_el[i] = closest;
 			//check if you intersect left and top
			if(x - 1 < (centers[2 * closest] - sx/2) && x >= (centers[2 * closest] - sx/2) && z >= lowZ)
			{
				//intersects on the left
				border[6 * i] = P;
				//add value to the left neighbour
				if(x != 0)
				{
					int ind = z * (a * b) + (x - 1) * b + y;
					border[6 * ind + 1] = P; 
				}
			}
			if(x == a - 1 && !periodic)
			{
				border[6 * i + 1] = 0; 
			}
			if(x == 0 && !periodic)
			{
				border[6 * i] = 0;
			}
			if(y + 1 >= (centers[2 * closest + 1] + sy/2) && y < (centers[2 * closest + 1] + sy/2) && z >= lowZ)
			{
				//intersects on the top
				border[6 * i + 3] = P;
				if(y != b - 1)
				{
					int ind = z * (a * b) + x * b + y + 1;
					border[6 * ind + 2] = P;
				}
			}
			if(y == 0 && !periodic)
			{
				border[6 * i + 2] = 0;
			}
			if(y == b - 1 && !periodic)
			{
				border[6 * i + 3] = 0;
			}
			if(z == 0)
			{
				border[6 * i + 4] = 0;
			}
			if(z == c - 1)
			{
				border[6 * i + 5] = 0;
			}

		}
	}
	if(polygon == 1)
	{
		float hexaX = sx/sqrt(3);
		float hexaY = sy/sqrt(3);
		for(int i = 0; i < a * b * c; ++i)
		{
			int z = i/(a * b);
			int temp = i - z * (a * b);
			int x = (temp / b);
			int y = (temp % b);	
			int closest = returnClosestCenter(x,y);
			closest_el[i] = closest;
			nodesPerHexa[closest].push_back(i);
			for(int j = 0; j < 6; ++j)
			{
				float angle = PI/6 + j * PI/3;
				float angleNext = angle + PI/3;
				float x1 = centers[2 * closest] + hexaX * cos(angle);
				float y1 = centers[2 * closest + 1] + hexaY * sin(angle);
				float x2 = centers[2 * closest] + hexaX * cos(angleNext);
				float y2 = centers[2 * closest + 1] + hexaY * sin(angleNext);
				float a_line = (y2 - y1)/(x2 - x1);
				float b_line = -1 * a_line * x1 + y1;
				float mySide = y - a_line * x - b_line;
				float leftNeighborSide = y - a_line * (x - 1) - b_line;
				float topNeighborSide = (y + 1) - a_line * x - b_line;
		
				if(j == 0 || j == 3)
				{
					if(((mySide >= 0 && leftNeighborSide < 0) || (mySide <= 0 && leftNeighborSide > 0)) && z >= lowZ)
					{
						float value = sin30 * P;
						border[6 * i] = value;
						if(x != 0)
						{ 
							int ind = z * (a * b) + (x - 1) * b + y;
							border[6 * ind + 1] = value; 
						}
						else if(periodic)
						{
							int ind = z * (a * b) + (a-1) * b + y;
							border[6 * ind + 1] = value; 
						}
						else
						{
							border[6 * i] = 0;
							int ind = z * (a * b) + (a-1) * b + y;
							border[6 * ind + 1] = 0; 
						}
					}	
					if(((mySide >= 0 && topNeighborSide < 0) || (mySide <= 0 && topNeighborSide > 0)) && z >= lowZ)
					{
						float value = sin60 * P;
						border[6 * i + 3] = value; 
						if(y != b - 1)
						{
							int ind = z * (a * b) + x * b + y + 1;
							border[6 * ind + 2] = value;
						}
						else if(periodic)
						{
							int ind = z * (a * b) + x * b;
							border[6 * ind + 2] = value;
						}
						else
						{
							border[6 * i + 3] = 0;
							int ind = z * (a * b) + x * b;
							border[6 * ind + 2] = 0;
						}
					}				
				}
				else if(j == 1 || j == 4)
				{
					//horizontal line for the oblique edges
					if(((mySide >= 0 && leftNeighborSide < 0) || (mySide <= 0 && leftNeighborSide > 0)) && z >= lowZ)
					{
						float value = sin30 * P;
						border[6 * i] = value;
						if(x != 0)
						{ 
							int ind = z * (a * b) + (x - 1) * b + y;
							border[6 * ind + 1] = value; 
						}
						else if(periodic)
						{
							int ind = z * (a * b) + (a-1) * b + y;
							border[6 * ind + 1] = value; 
						}
						else
						{
							border[6 * i] = 0;
							int ind = z * (a * b) + (a-1) * b + y;
							border[6 * ind + 1] = 0; 
						}
					}
					//vertical line for the oblique edges
					if(((mySide >= 0 && topNeighborSide < 0) || (mySide <= 0 && topNeighborSide > 0)) && z >= lowZ)
					{
						float value = sin60 * P;
						border[6 * i + 3] = value; 
						if(y != b - 1)
						{
							int ind = z * (a * b) + x * b + y + 1;
							border[6 * ind + 2] = value;
						}
						else if(periodic)
						{
							int ind = z * (a * b) + x * b;
							border[6 * ind + 2] = value;
						}
						else
						{
							border[6 * i + 3] = 0;
							int ind = z * (a * b) + x * b;
							border[6 * ind + 2] = 0;
						}
					}	
					
				}
				else if(j == 2 || j == 5)
				{
				
					if(x >= x1 && (x-1) < x1  && z >= lowZ)
					{
						border[6 * i] = P; 
						if(x != 0)
						{ 
							int ind = z * (a * b) + (x - 1) * b + y;
							border[6 * ind + 1] = P; 
						}
						else if(periodic)
						{
							int ind = z * (a * b) + (a-1) * b + y;
							border[6 * ind + 1] = P; 
						}
						else
						{
							border[6 * i] = 0;
							int ind = z * (a * b) + (a-1) * b + y;
							border[6 * ind + 1] = 0; 
						}	
					}
				}
			}
			if(z == 0)
			{
				border[6 * i + 4] = 0;
			}
			if(z == c - 1)
			{
				border[6 * i + 5] = 0;
			}	
			if(x == a - 1 && !periodic)
			{
				border[6 * i + 1] = 0; 
			}
			if(x == 0 && !periodic)
			{
				border[6 * i] = 0;
			}	
			if(y == 0 && !periodic)
			{
				border[6 * i + 2] = 0;
			}
			if(y == b - 1 && !periodic)
			{
				border[6 * i + 3] = 0;
			}				
			
		}
	}
} 

//border thickness contains more than one element
int getPixelsP0btMoreThanOne(int polygon, int maxZ, int lowZ)
{
	border = (float*)malloc(6 * a * b * c * sizeof(float));
	closest_el = (int*)malloc(a * b * c * sizeof(int));
	int numberOfNodesToSimulate = 0;
	for(int i = 0 ; i < 6 * a * b * c; ++i)
	{
		border[i] = 1/h;
	}		
	if (polygon == 0)
	{	
		for(int i = 0; i < a * b * c; ++i)
		{
			int z = i/(a * b);
			int temp = i - z * (a * b);
			int x = (temp / b);
			int y = (temp % b);	
			int closest = returnClosestCenter(x,y);
			closest_el[i] = closest;
 			double dl = abs(x - centers[2 * closest] + sx/2);
 			double dto = abs(y - centers[2 * closest + 1] - sy/2);
 			//left
			if(dl <= bt)
			{
				if(z >= lowZ)
				{
					//intersects on the left
					if(dl + 1 <= bt)
					{
						border[6 * i] = -2;
					}
					else
					{
						border[6 * i] = 0;
					}
				}
				else if(z == lowZ - 1 && dl + 1 <= bt)
				{
					border[6 * i + 5] = 0;
				}
			}
			//right
			if(sx - dl <= bt)
			{
				if(z >= lowZ)
				{
					if(sx - dl + 1 <= bt)
					{
						border[6 * i] = -2;
					}
					else
					{
						border[6 * i + 1] = 0;
					}
				}	
				else if(z == lowZ - 1 && sx - dl + 1 <= bt)	
				{
					border[6 * i + 5] = 0;
				}		
			}
			//right outer border
			if(x >= a - 1 - bt && !periodic)
			{
				if(x - 1 >= a - 1 - bt)
				{
					border[6 * i] = -2;					
				}
				else
				{
					border[6 * i + 1] = 0; 
				}				
			}
			//left outer border
			if(x <= bt && !periodic)
			{
				if(x + 1 <= bt)
				{
					border[6 * i] = -2;					
				}
				else
				{
					border[6 * i] = 0;
				}
			}			
			//top
			if(dto <= bt)
			{
				if(z >= lowZ)
				{
					if(dto + 1 <= bt)
					{
						border[6 * i] = -2;					
					}
					else
					{
						border[6 * i + 3] = 0;					
					}
				}
				else if(z == lowZ - 1 && dto + 1 <= bt)
				{
					border[6 * i + 5] = 0;
				}
			}
			//bottom 
			if(sy - dto <= bt)
			{
				if(z >= lowZ)
				{
					if(sx - dto + 1 <=bt)
					{
						border[6 * i] = -2;				
					}
					else
					{
						border[6 * i + 2] = 0;				
					}
				}
				else if(z == lowZ - 1 && sx - dto + 1 <=bt)
				{
					border[6 * i + 5] = 0;					
				}
			}
			if(y <= bt && !periodic)
			{
				if(y + 1 <= bt)
				{
					border[6 * i] = -2;
				}
				else
				{
					border[6 * i + 2] = 0;
				}
			}
			if(y >= b - 1 - bt && !periodic)
			{
				if(y - 1 >= b - 1 - bt)
				{
					border[6 * i] = -2;				
				}
				else
				{
					border[6 * i + 3] = 0;
				}
			}
			if(z == 0)
			{
				border[6 * i + 4] = 0;
			}
			if(z == c - 1)
			{
				border[6 * i + 5] = 0;
			}
			if(border[6 * i] != -2) numberOfNodesToSimulate++;

		}
	}
	if(polygon == 1)
	{
		float hexaX = sx/sqrt(3);
		float hexaY = sy/sqrt(3);
		float maxDistance =  0.7071 * border_thickness_in_elements;
		for(int i = 0; i < a * b * c; ++i)
		{
			int z = i/(a * b);
			int temp = i - z * (a * b);
			int x = (temp / b);
			int y = (temp % b);	
			int closest = returnClosestCenter(x,y);
			closest_el[i] = closest;
			for(int j = 0; j < 6; ++j)
			{
				float angle = PI/6 + j * PI/3;
				float angleNext = angle + PI/3;
				float x1 = centers[2 * closest] + hexaX * cos(angle);
				float y1 = centers[2 * closest + 1] + hexaY * sin(angle);
				float x2 = centers[2 * closest] + hexaX * cos(angleNext);
				float y2 = centers[2 * closest + 1] + hexaY * sin(angleNext);
				float distance = abs((y2 - y1)*x - (x2 - x1)*y + x2*y1 - y2*x1)/sqrt((y2 - y1)*(y2 - y1) + (x2-x1)*(x2-x1));
				if(distance <= maxDistance && z >= lowZ)
				{
					border[6 * i] = -2;
					break;
				}
			}
			if(z == 0)
			{
				border[6 * i + 4] = 0;
			}
			if(z == c - 1)
			{
				border[6 * i + 5] = 0;
			}
			if(x < maxDistance && !periodic)
			{
				border[6 * i] = -2;
			}	
			if(y < maxDistance && !periodic)
			{
				border[6 * i] = -2;
			}		
			if(border[6 * i] != -2)
			{
				numberOfNodesToSimulate++;	
				nodesPerHexa[closest].push_back(i);
			}					
			
		}
		
	}
	return numberOfNodesToSimulate;
} 

bool createRandomInitialConditions(double randomIntervalSize)
{
	cout << "Creating random initial conditions using concentrations of black" << endl;
	for(int i = 0; i < a * b * c; ++i)
	{
		u_old[i] = initialBlack[0] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
		v_old[i] = initialBlack[1] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
		w_old[i] = initialBlack[2] + ((double)rand()/(RAND_MAX) * 2 - 1) * randomIntervalSize;
	}
	return true;
}

