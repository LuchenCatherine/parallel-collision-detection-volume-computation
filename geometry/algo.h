#ifndef ALGO_H
#define ALGO_H

#include <iostream>
#include <vector>
#include<chrono>
#include "omp.h"

#include "mymesh.h"
#include "utils.h"

std::vector<std::pair<std::string, double>> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points);

std::vector<std::pair<std::string, double>> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points, GeoInfo &geo_info, int resolution);

std::vector<std::string> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue);

std::vector<std::string> collision_detection_single_tissue_parallel(std::vector<Mymesh> &organ, Mymesh &tissue);

std::vector<std::string> collision_detection_brute_force(std::vector<Mymesh> &organ, Mymesh &tissue);

bool intersection_brute_force(std::vector<Triangle> &faces1, std::vector<Triangle> &faces2);

double voxelization_volume_computation(GeoInfo &geo_info, Mymesh &anatomical_structure, int resolution);

double voxelization_volume_computation_openmp(GeoInfo &geo_info, Mymesh &anatomical_structure, int resolution);

#endif