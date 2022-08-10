#pragma once
#include <sstream>
#include <iostream>
#include<fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>
#include <math.h>
#include<boost/filesystem.hpp>

#include "mymesh.h"


struct SpatialEntity {

    SpatialEntity() = default;
    SpatialEntity(const std::string &aso, const std::string &sse, const std::string &nn, const std::string &lb, const std::string &ro, const std::string &gf): anatomical_structure_of(aso), source_spatial_entity(sse), node_name(nn), label(lb), representation_of(ro), glb_file(gf) {};

    std::string anatomical_structure_of;
    std::string source_spatial_entity;
    std::string node_name;
    std::string label;
    std::string representation_of;
    std::string glb_file;


};

struct GeoInfo {

    GeoInfo() = default;
    GeoInfo(double minx, double miny, double minz, double maxx, double maxy, double maxz, Eigen::Matrix3d rot, Eigen::Vector3d tra, Eigen::Vector3d ori): 
    min_x(minx), min_y(miny), min_z(minz), max_x(maxx), max_y(maxy), max_z(maxz), R(rot), T(tra), origin(ori) {};

    double min_x;
    double min_y;
    double min_z;

    double max_x;
    double max_y;
    double max_z;

    Eigen::Matrix3d R;
    Eigen::Vector3d T;
    Eigen::Vector3d origin;

};

//origins, meshes
void load_all_organs(const std::string &body_path, std::unordered_map<std::string, std::vector<Mymesh>> &total_body);

// generate origins from the whole model without any hole-filling.
void gen_origin(const std::string &whole_model_path, std::unordered_map<std::string, Eigen::Vector3d> &organ_origins);

//including x_scaling, x_rotation, x_translation, x_origin 
void tissue_transform(std::unordered_map<std::string, double> &params, Surface_mesh &tissue_mesh, std::vector<Point> &points, int resolution);

void tissue_transform(GeoInfo &geo_info, Surface_mesh &tissue_mesh);

GeoInfo extract_params(std::unordered_map<std::string, double> &params);

std::string organ_split(const std::string &url);

void load_ASCT_B(const std::string &file_path, std::unordered_map<std::string, std::string> &mapping, std::unordered_map<std::string, SpatialEntity> &mapping_node_spatial_entity);

void print_result(std::string &algorithm, std::vector<std::string> &result);