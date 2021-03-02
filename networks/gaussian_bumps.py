import numpy as np
import random
import pickle
import sys
import mesh

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)

def closest_node_dist(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.min(dist_2)

def gaussian2D(x, y, mu, sig):
	return  np.exp(-1 * (np.power(x - mu[0], 2.) / (2 * np.power(sig, 2.)) + np.power(y - mu[1], 2.) / (2 * np.power(sig, 2.))))

def gauss_for_plotting(x, y, sig, centers, h, a):
	c = centers[closest_node([x,y], centers)] 
	return max(h + a * gaussian2D(x, y, c, sig),0)
	

def dont_simulate(x, y, centers, max_dist):
	c = closest_node_dist([x,y], centers)
	if c > max_dist:
		return True
	else:
		return False
		

gauss_vectorized = np.vectorize(gauss_for_plotting)
gauss_vectorized.excluded.add(2)
gauss_vectorized.excluded.add(3)
gauss_vectorized.excluded.add(4)
gauss_vectorized.excluded.add(5)


dont_simulate_vect = np.vectorize(dont_simulate)
dont_simulate_vect.excluded.add(2)
dont_simulate_vect.excluded.add(3)


def gaussian_hexa(res_x, res_y, scales_x, scales_y, height_bottom, sigma, max_height, ply_file, plot = False, export_centers = False):
	dist_x = 1.0 * res_x/scales_x
	centers = []
	#res_y is changed
	dist_y = dist_x /np.sqrt(3)
	for i in range(scales_y):
		start_x = dist_x/2 * ((i+1)%2) + dist_x * ((i+2)%2)
		final = scales_y - ((i+2)%2)
		for j in range(final):
			centers.append([start_x + j * dist_x, dist_x/2 + i * 1.5 * dist_y])	
	res_y = int(round(max([i[1] for i in centers]) + dist_x/2 + 1))
	maximum = gaussian2D(0,0, [0,0], sigma)
	minimum = gaussian2D(dist_x, centers[0][1], centers[0], sigma)
	a = (height_bottom - max_height)/(minimum - maximum)
	height_bottom_ = max_height - a * maximum
	X, Y = np.meshgrid(np.arange(0, res_x, 1), np.arange(0, res_y, 1))
	Z = gauss_vectorized(X, Y, sigma, centers, height_bottom_, a)
	Z[np.where(Z < height_bottom + 0.00001)] = height_bottom + 0.00001
	#save the layer in a ply file
	f = open(ply_file, "w")
	f.write("ply\n")
	f.write("format ascii 1.0\n")
	f.write("element vertex " + str(res_x * res_y) + "\n")
	f.write("property float x\n")
	f.write("property float y\n")
	f.write("property float z\n")
	f.write("element face " + str((res_x - 1) * (res_y - 1)) + "\n")
	f.write("property list uchar int vertex_indices\n")
	f.write("end_header\n")	
	for i in range(Z.shape[1]):
		for j in range(Z.shape[0]):
			f.write(str(i) + " " + str(j) + " " + str(Z[j][i]) + "\n")
	for i in range(Z.shape[1]-1):
		for j in range(Z.shape[0]-1):
			index = i * Z.shape[0] + j
			f.write("4 ")
			f.write(str(index) + " " + str(index + 1) + " " + str(index + Z.shape[0] + 1) + " " + str(index + Z.shape[0]) + "\n")
	f.close()
	exclude = dont_simulate_vect(X,Y, centers,(dist_y + 1)**2)
	Z[exclude] = -100
	return Z.transpose()

def flat_layer(res_x, res_y, height):
	return np.zeros(shape=(res_x, res_y)) + height

#predefine the dimensions in x and y directions and the number of elements
#bottom eight, top height and sigma can be set

def main():
	if len(sys.argv) < 4:
		print("Command line parameters")
		print("Argument 1: Maximum network height (integer)") 
		print("Argument 2: Domain thickness between gaussian bumps (integer)")
		print("Argument 3: Standard deviation of the bumps (integer)")
		exit()
	maxZ = int(sys.argv[1])
	spacing = int(sys.argv[2])
	sigma = int(sys.argv[3])

	name = "gauss10by10_" + str(spacing) + "_" + str(maxZ) + "_"+ str(sigma)
	dim_x = 200
	dim_y = 173
	scales_x = 10
	scales_y = 10
	top_layer = gaussian_hexa(dim_x, dim_y, scales_x, scales_y, spacing, sigma, maxZ, name + "_top_layer.ply")
	bottom_layer = flat_layer(top_layer.shape[0], top_layer.shape[1], 0)
	mesh.create_points(top_layer.shape[0], top_layer.shape[1], bottom_layer, top_layer, name + ".ply", name + ".txt")

if __name__ == "__main__":
    main()