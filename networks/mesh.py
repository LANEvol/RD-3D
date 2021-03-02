import numpy as np

def flat_layer(res_x, res_y, height):
	return np.zeros(shape=(res_x, res_y)) + height

def is_below(layer, x, y, z):
	if layer[x][y] >= z:
		return True
	else:
		return False

def is_above(layer, x, y, z):
	if layer[x][y] <= z:
		return True
	else:
		return False

#0 x, 1 y, 2 z
def intersection(bottom_layer, top_layer, x, y, z, direction):
	to_return = [False, False]
	if direction == 0:
		if x == 0:
			to_return[0] = True
		elif not is_below(top_layer, x - 1, y, z) or not is_above(bottom_layer, x - 1, y, z):
			to_return[0] = True
		elif top_layer[x-1][y] < -1:
			to_return[0] = True
		if x == top_layer.shape[0] - 1:
			to_return[1] = True
		elif not is_below(top_layer, x + 1, y, z) or not is_above(bottom_layer, x + 1, y, z):
			to_return[1] = True		
		elif top_layer[x+1][y] < -1:
			to_return[1] = True
	elif direction == 1:
		if y == 0:
			to_return[0] = True
		elif not is_below(top_layer, x, y - 1, z) or not is_above(bottom_layer, x, y - 1, z):
			to_return[0] = True
		elif top_layer[x][y-1] < -1:
			to_return[0] = True		
		if y == top_layer.shape[1] - 1:
			to_return[1] = True
		elif not is_below(top_layer, x, y + 1, z) or not is_above(bottom_layer, x, y + 1, z):
			to_return[1] = True	
		elif top_layer[x][y+1] < -1:
			to_return[1] = True
	elif direction == 2:
		if not is_below(top_layer, x, y, z - 1) or not is_above(bottom_layer, x, y, z - 1):
			to_return[0] = True
		if not is_below(top_layer, x, y, z + 1) or not is_above(bottom_layer, x, y, z + 1):
			to_return[1] = True	
	return to_return	


def save_ply_file(points, ply_file, dim_x, dim_y, dim_z):
	f = open(ply_file, "w")
	f.write("ply\n")
	f.write("format ascii 1.0\n")
	f.write("element vertex " + str(dim_x * dim_y * dim_z) + "\n")
	f.write("property float x\n")
	f.write("property float y\n")
	f.write("property float z\n")
	f.write("property uchar red\n")
	f.write("property uchar green\n")
	f.write("property uchar blue\n")
	f.write("end_header\n")
	for i in range(dim_x):
		for j in range(dim_y):
			for k in range(dim_z):
				f.write(str(i) + " " + str(j) + " " + str(k) + " ")
				r = 0
				g = 0
				b = 0
				if points[i][j][k][0] == 1 or points[i][j][k][1] == 1:
					r = 255
				if points[i][j][k][2] == 1 or points[i][j][k][3] == 1:
					b = 255
				if points[i][j][k][4] == 1 or points[i][j][k][5] == 1:
					g = 255
				f.write(str(r) + " " + str(g) + " " + str(b) + "\n")
	f.close()

def save_txt_file(points, txt_file, dim_x, dim_y, dim_z):
	f = open(txt_file, "w")
	f.write(str(dim_x) + " " + str(dim_y) + " " + str(dim_z) + "\n")
	for i in range(dim_x):
		for j in range(dim_y):
			for k in range(dim_z):
				f.write(str(i) + " " + str(j) + " " + str(k))
				for r in range(6):
					f.write(" " + str(points[i][j][k][r]))					
				f.write("\n")
	f.close()


def create_points(res_x, res_y, bottom_layer, top_layer, ply_file, txt_file):
	maximum = int(round(np.max(top_layer)))
	points = np.zeros(shape = (bottom_layer.shape[0], bottom_layer.shape[1], maximum, 6), dtype = "int16")

	for i in range(res_x):
		for j in range(res_y):
			for k in range(maximum):
				temp = np.zeros(6) - 1
				if not is_above(bottom_layer, i, j, k) or not is_below(top_layer, i, j, k):
					temp = temp - 1
					points[i][j][k] = temp
					continue
				if top_layer[i][j] < -1:
					temp = temp - 1
					points[i][j][k] = temp
					continue					
				intersections = [intersection(bottom_layer, top_layer, i, j, k, n) for n in range(3)]
				for t in range(6):
					if intersections[int(t/2)][int(t%2)]:
						temp[t] = 1
				points[i][j][k] = temp
	save_txt_file(points, txt_file, bottom_layer.shape[0], bottom_layer.shape[1], maximum)

	return points