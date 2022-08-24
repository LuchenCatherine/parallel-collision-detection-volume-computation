#include "algo.h"
#include "utils.h"

#include <cpprest/http_listener.h>
#include <cpprest/json.h>
#pragma comment(lib, "cpprest_2_10")

using namespace web;
using namespace web::http;
using namespace web::http::experimental::listener;

#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include<chrono>

//global variables
extern std::map<utility::string_t, utility::string_t> dictionary;
extern std::unordered_map<std::string, Eigen::Vector3d> organ_origins;                     //origins of organs
extern std::unordered_map<std::string, std::string> mapping;                               //mapping from standard organ name(e.g., #VHFLeftKidney) to glb file name without suffix(e.g., VH_F_Kidney_L)
extern std::unordered_map<std::string, std::vector<Mymesh>> total_body;                    //mapping from organ name(glb file name) to vector of meshes of a certain organ
extern std::unordered_map<std::string, SpatialEntity> mapping_node_spatial_entity;         // mapping from AS to its information in asct-b table 



//display json
void display_json(json::value const & jvalue, utility::string_t const & prefix);

//parse json
void parse_json(json::value const &jvalue, json::value &answer);

// construct response
void construct_response(json::value &answer, std::vector<std::pair<std::string, double>> &result, double tissue_volume, double elapsed_time);

// handle get request
void handle_get(http_request request);

// handle post request
void handle_post(http_request request);

// handle request by action
void handle_request(http_request request, std::function<void(json::value const &, json::value &)> action);

void init(std::unordered_map<std::string, std::vector<Mymesh>> &total_body);