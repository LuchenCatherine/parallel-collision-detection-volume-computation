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

    if (argc < 5)
    {
        std::cout << "Please provide the organ_origins_file_path, asct_b_file_path, body_path(model_path) and tissue path!" << std::endl;
        return 0;
    }

    std::string organ_origins_file_path = std::string(argv[1]);
    std::string asct_b_file_path = std::string(argv[2]);
    std::string body_path = std::string(argv[3]);
    std::string tissue_path = std::string(argv[4]);
    // load origins
    gen_origin(organ_origins_file_path, organ_origins);
    load_ASCT_B(asct_b_file_path, mapping, mapping_node_spatial_entity);
    load_all_organs(body_path, total_body);

    init(total_body);


    std::unordered_map<std::string, std::vector<Mymesh>> tissue_maps;  

    auto t1 = std::chrono::high_resolution_clock::now();

    for (fs::directory_entry& organ_path : fs::directory_iterator(tissue_path)) 
    {

        std::string organ_name = organ_path.path().stem().string();
           
        for (fs::directory_entry& tissue : fs::directory_iterator(organ_path)) 
        {
            std::string file_path = tissue.path().string();
            tissue_maps[organ_name].push_back(Mymesh(file_path));
        }

        for (auto &tissue: tissue_maps[organ_name]){
            tissue.create_aabb_tree();
        } 

        for (auto &tissue: tissue_maps[organ_name]) tissue.extract_faces();
    }

    auto t2 = std::chrono::high_resolution_clock::now();   
    for (auto &organ_tissues: tissue_maps)
    {
        std::string organ_file_name = organ_tissues.first;
        std::vector<Mymesh> tissue_meshes = organ_tissues.second;
        std::cout << "tissues of " << organ_file_name << std::endl;
        double elapsed_time = 0;
        for (auto &tissue: tissue_meshes)
        {
            auto result = collision_detection_boolean_operation(total_body[organ_file_name], tissue, elapsed_time);
        }
        
    }
    auto t3 = std::chrono::high_resolution_clock::now();

    for (auto &organ_tissues: tissue_maps)
    {
        std::string organ_file_name = organ_tissues.first;
        std::vector<Mymesh> tissue_meshes = organ_tissues.second;
        std::cout << "tissues of " << organ_file_name << std::endl;
        double elapsed_time = 0;
        for (auto &tissue: tissue_meshes)
        {
            auto result = collision_detection_single_tissue(total_body[organ_file_name], tissue, elapsed_time);
        }
        
    }

    auto t4 = std::chrono::high_resolution_clock::now();




    std::chrono::duration<double> duration1 = t2 - t1;
    std::chrono::duration<double> duration2 = t3 - t2;
    std::chrono::duration<double> duration3 = t4 - t3;
    
    std::cout << "loading tissues: " << duration1.count() << " boolean operation: " << duration2.count() << " collision detection: " << duration3.count() << std::endl;

}