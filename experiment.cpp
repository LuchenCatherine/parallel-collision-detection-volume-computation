#include "algo.h"
#include "utils.h"

#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <chrono>

//global variables
std::unordered_map<std::string, Eigen::Vector3d> organ_origins;                     //origins of organs
std::unordered_map<std::string, std::string> mapping;                               //mapping from standard organ name(e.g., #VHFLeftKidney) to glb file name without suffix(e.g., VH_F_Kidney_L)
std::unordered_map<std::string, std::vector<Mymesh>> total_body;                    //mapping from organ name(glb file name) to vector of meshes of a certain organ
std::unordered_map<std::string, SpatialEntity> mapping_node_spatial_entity;         // mapping from AS to its information in asct-b table 

void init(std::unordered_map<std::string, std::vector<Mymesh>> &total_body)
{
   std::cout << "############### init start ###############\n";
   Mymesh test_tissue("/home/catherine/data/model/test_tissue.off");
   test_tissue.create_aabb_tree();
   test_tissue.extract_faces();

   for (auto &organ_item: total_body) collision_detection_single_tissue(organ_item.second, test_tissue);
   std::cout << "############### init completed ###############\n";

}

int main(int argc, char **argv)
{

    if (argc < 6)
    {
        std::cout << "Please provide the organ_origins_file_path, asct_b_file_path, body_path(model_path) server IP and port number!" << std::endl;
        return 0;
    }

    std::string organ_origins_file_path = std::string(argv[1]);
    std::string asct_b_file_path = std::string(argv[2]);
    std::string body_path = std::string(argv[3]);
    std::string server_ip = std::string(argv[4]);
    std::string port = std::string(argv[5]);

    // load origins
    gen_origin(organ_origins_file_path, organ_origins);
    load_ASCT_B(asct_b_file_path, mapping, mapping_node_spatial_entity);
    load_all_organs(body_path, total_body);

    init(total_body);

    Mymesh boolean_test_tissue("/home/catherine/Research/learningcpp/parallel-collision-detection-volume-computation/tissue_mesh.off");
    boolean_test_tissue.create_aabb_tree();
    boolean_test_tissue.extract_faces();

    double elapsed_time = 0;
    std::string organ_file_name = "VH_F_Kidney_L";
    auto x = total_body[organ_file_name];
    auto result = collision_detection_boolean_operation(total_body[organ_file_name], boolean_test_tissue, elapsed_time);
    std::cout << elapsed_time << std::endl;

}