#ifndef ALGO_H
#define ALGO_H

#include <iostream>
#include <vector>
#include<chrono>
#include "omp.h"

#include "mymesh.h"
#include "utils.h"

#include "RTree.h"

typedef RTree<int, double, 3, double> rtree_3;
typedef int ValueType;

std::vector<std::string> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue);

std::vector<std::string> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, double &elapsed_time);

std::vector<std::pair<std::string, double>> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points);

std::vector<int> helper_collision_detection_single_tissue_return_index(std::vector<Mymesh> &organ, Mymesh &tissue);

std::vector<std::pair<std::string, double>> collision_detection_volume_computation_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, GeoInfo &geo_info, int resolution, std::string method, double &elapsed_time);

std::vector<std::pair<std::string, double>> collision_detection_boolean_operation(std::vector<Mymesh> &organ, Mymesh &tissue, double &elapsed_time);

std::vector<std::pair<std::string, double>> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points, GeoInfo &geo_info, int resolution);

double voxelization_volume_computation(GeoInfo &geo_info, Mymesh &anatomical_structure, int resolution);

double voxelization_volume_computation_openmp(GeoInfo &geo_info, Mymesh &anatomical_structure, int resolution);

std::vector<std::string> collision_detection_single_tissue_parallel(std::vector<Mymesh> &organ, Mymesh &tissue);

std::vector<std::string> collision_detection_brute_force(std::vector<Mymesh> &organ, Mymesh &tissue, double &elapsed_time);

std::vector<std::string> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, rtree_3 &rtree, double &elapsed_time);

std::vector<std::string> collision_detection_with_rtree(std::vector<Mymesh> &organ, Mymesh &tissue, rtree_3 &rtree, double &elapsed_time);

bool intersection_brute_force(std::vector<Triangle> &faces1, std::vector<Triangle> &faces2);

bool MySearchCallback(ValueType id);


#endif