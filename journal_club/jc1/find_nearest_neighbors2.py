'''
script to compare distances between two groups of points. Specifically, for each
point in group 1, finds the point in group 2 that is the smallest distance away
and outputs the distances and names of the points in group 2 that meet these
criteria for the points in group 1

This version assumes a csv file called 'all_cells.csv'
'''

import os

def x_y_cluster_columns(title_line):
	for name_number, name in enumerate(title_line):
		if name=='x':
			x_column=name_number
		if name=='y':
			y_column=name_number
		if name=='cluster':
			cluster_column=name_number
	return x_column, y_column, cluster_column

def modify_cluster_dict(line, cluster_dict, cluster_column):
	cell, cluster=line[0],line[cluster_column]
	if cluster not in cluster_dict:
		cluster_dict[cluster]=[]
	cluster_dict[cluster].append(cell)
	return cluster_dict

def get_coords(coord_file):
	coord_dict, cluster_dict={},{}
	for line_number, line in enumerate(open(coord_file)):
		line=line.strip().replace('"', '').split(',')
		if line_number==0:
			x_column, y_column, cluster_column=x_y_cluster_columns(line)
		else:
			coord_dict[line[0]]=[float(line[x_column]), float(line[y_column])]
			cluster_dict=modify_cluster_dict(line, cluster_dict, cluster_column)
	return coord_dict, cluster_dict

def compute_distances(cluster_dict, coord_dict, cluster1, cluster2):
	distance_dict={}
	print(len(cluster_dict[cluster1]))
	print(len(cluster_dict[cluster2]))
	for cell_number, cell1 in enumerate(cluster_dict[cluster1]):
		distance_dict[cell1]=['nonsense', float('Inf')]
		if cell_number%100==0:
			print(cell_number/len(cluster_dict[cluster1]))
		for cell2 in cluster_dict[cluster2]:
			x_distance=abs(coord_dict[cell1][0]-coord_dict[cell2][0])
			y_distance=abs(coord_dict[cell1][1]-coord_dict[cell2][1])
			total_distance=(x_distance**2+y_distance**2)**0.5
			if total_distance<distance_dict[cell1][1]:
				distance_dict[cell1]=[cell2, total_distance]
	return distance_dict
#for pair in pair_file:

def output_distances(distance_dict, first_cluster, second_cluster):
	os.makedirs('python_nearest_neighbor_outputs', exist_ok=True)
	output_path = 'python_nearest_neighbor_outputs/'+first_cluster+'_'+second_cluster+'_nearest_neighbors.csv'
	with open(output_path, 'w') as output_file:
		for cell1 in distance_dict:
			output_file.write(','.join([cell1]+list(map(str, distance_dict[cell1])))+'\n')

coord_dict, cluster_dict=get_coords('all_cells.csv')
valid_clusters=sorted(list(cluster_dict.keys()))
for first_cluster in valid_clusters:
	for second_cluster in valid_clusters:
		if first_cluster!=second_cluster:
			print(first_cluster, second_cluster)
			distance_dict=compute_distances(cluster_dict, coord_dict, first_cluster, second_cluster)
			output_distances(distance_dict, first_cluster, second_cluster)