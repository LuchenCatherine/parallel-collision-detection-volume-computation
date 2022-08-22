#include "service.h"


//display json
void display_json(
   json::value const & jvalue,
   utility::string_t const & prefix)
{
   std::cout << prefix << jvalue.serialize() << std::endl;
}

//parse json
void parse_json(json::value const &jvalue, json::value &answer)
{

   std::unordered_map<std::string, double> params;


   auto target = jvalue.at("target").as_string();
   auto refence_organ_name = organ_split(target); 
   std::cout << "target url: " << target << " target: " << refence_organ_name << " " << std::endl;
   
   // only test for kidneys, will test other organs soon.
   // if (!(refence_organ_name == "#VHFLeftKidney" || refence_organ_name == "#VHFRightKidney" || refence_organ_name == "#VHMLeftKidney" || refence_organ_name == "#VHMRightKidney"))
   // {
   //    answer[U("error_message")] = json::value::string(U("only test tissue blocks in kidneys"));
   //    return;
   // }

   // test for all organs
   if (mapping.find(refence_organ_name) == mapping.end()) 
   {
      std::cout << refence_organ_name << " doesn't exist in ASCT-B table!" << std::endl;
      return;
   }
   
   //extract parameters from json request
   params["x_dimension"] = jvalue.at("x_dimension").as_double();
   params["y_dimension"] = jvalue.at("y_dimension").as_double();
   params["z_dimension"] = jvalue.at("z_dimension").as_double();
   params["x_scaling"] = jvalue.at("x_scaling").as_double();
   params["y_scaling"] = jvalue.at("y_scaling").as_double();
   params["z_scaling"] = jvalue.at("z_scaling").as_double();
   params["x_translation"] = jvalue.at("x_translation").as_double();
   params["y_translation"] = jvalue.at("y_translation").as_double();
   params["z_translation"] = jvalue.at("z_translation").as_double();
   params["x_rotation"] = jvalue.at("x_rotation").as_double();
   params["y_rotation"] = jvalue.at("y_rotation").as_double();
   params["z_rotation"] = jvalue.at("z_rotation").as_double();
   
   auto algorithm = jvalue.at("algorithm").as_string();
   int resolution = jvalue.at("resolution").as_integer();

   std::string organ_file_name = mapping[organ_split(target)];
   std::cout << "target url: " << target << " target: " << refence_organ_name << " " << "organ file name: " << organ_file_name << std::endl;

   Eigen::Vector3d origin = organ_origins[refence_organ_name];
   params["x_origin"] = origin(0);
   params["y_origin"] = origin(1);
   params["z_origin"] = origin(2);

   
   Surface_mesh tissue_mesh;
   std::vector<Point> points; //center of voxels inside the tissue block
   tissue_transform(params, tissue_mesh, points, 10);
   // GeoInfo geo_info = extract_params(params);
   // tissue_transform(geo_info, tissue_mesh);


   Mymesh my_tissue(tissue_mesh);
   my_tissue.create_aabb_tree();
   my_tissue.extract_faces();

   //core function
   if (algorithm == "production")
   {
      std::vector<std::pair<std::string, double>> result = collision_detection_single_tissue(total_body[organ_file_name], my_tissue, points);

      // std::vector<std::pair<std::string, double>> result = collision_detection_single_tissue(total_body[organ_file_name], my_tissue, points, geo_info, resolution);
      //print result
      std::cout << "result:\nlabel         percentage" << std::endl;
      for (auto s: result) {std::cout << s.first << " " << s.second << std::endl;}

      //construct the response
      double tissue_volume = params["x_dimension"] * params["y_dimension"] * params["z_dimension"];
      for (int i = 0; i < result.size(); i++)
      {
         auto res = result[i];
         json::value AS;
         
         auto node_name = res.first;
         SpatialEntity &se = mapping_node_spatial_entity[node_name];
         AS[U("node_name")] = json::value::string(U(node_name));
         AS[U("label")] = json::value::string(U(se.label));
         AS[U("representation_of")] = json::value::string(U(se.representation_of));
         AS[U("id")] = json::value::string("http://purl.org/ccf/latest/ccf.owl" + se.source_spatial_entity + "_" + node_name);

         if (res.second < 0)  
         {
            AS[U("percentage")] = json::value(0);
            AS[U("is_closed")] = json::value(false);
         }
         else
         {
            AS[U("percentage")] = json::value(res.second);
            AS[U("volume")] = json::value(res.second * tissue_volume);
            AS[U("is_closed")] = json::value(true);
         }
         
         answer[i] = AS;
      }
   }

   else if (algorithm == "experiment_normal")
   {
      auto t1 = std::chrono::high_resolution_clock::now();

      for (auto &organ_item: total_body)
      {
         auto result = collision_detection_single_tissue(organ_item.second, my_tissue);
      }

      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration2 = t2 - t1;
      std::cout << "[test all organs] normal running time is " << duration2.count() << " seconds" << std::endl;  

   }
   else if (algorithm == "experiment_normal_single")
   {
      auto t1 = std::chrono::high_resolution_clock::now();
      auto result = collision_detection_single_tissue(total_body[organ_file_name], my_tissue);
      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration2 = t2 - t1;
      std::cout << "normal running time is " << duration2.count() << " seconds" << std::endl;  

      print_result(algorithm, result);

   }
   else if (algorithm == "experiment_parallel")
   {
      auto t1 = std::chrono::high_resolution_clock::now();
      
      for (auto &organ_item: total_body)
      {
         auto result = collision_detection_single_tissue_parallel(organ_item.second, my_tissue);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration2 = t2 - t1;
      std::cout << "[test all organs] parallel running time is " << duration2.count() << " seconds" << std::endl;  

   }
   else if (algorithm == "experiment_parallel_single")
   {
      auto t1 = std::chrono::high_resolution_clock::now();
      auto result = collision_detection_single_tissue_parallel(total_body[organ_file_name], my_tissue);
      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration2 = t2 - t1;
      std::cout << "parallel running time is " << duration2.count() << " seconds" << std::endl;  

      print_result(algorithm, result);
   }
   else if (algorithm == "brute_force")
   {
      auto t1 = std::chrono::high_resolution_clock::now();
      
      for (auto &organ_item: total_body)
      {
         auto result = collision_detection_brute_force(organ_item.second, my_tissue);
      }

      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration2 = t2 - t1;
      std::cout << "[test all organs] brute force running time is " << duration2.count() << " seconds" << std::endl;  

   }
   else if (algorithm == "brute_force_single")
   {
      auto t1 = std::chrono::high_resolution_clock::now();
      
      // for (auto &organ_item: total_body)
      // {
      //    auto result = collision_detection_brute_force(organ_item.second, my_tissue);
      //    print_result(algorithm, result);
      // }

      auto result = collision_detection_brute_force(total_body[organ_file_name], my_tissue);
      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration2 = t2 - t1;
      std::cout << "brute force running time is " << duration2.count() << " seconds" << std::endl;  

      print_result(algorithm, result);

   }

}

// handle get request
void handle_get(http_request request)
{

   std::cout << "\nhandle GET" << std::endl;

   auto answer = json::value::object();

   for (auto const & p : dictionary)
   {
      answer[p.first] = json::value::string(p.second);
   }

   display_json(json::value::null(), "R: ");
   display_json(answer, "S: ");

   request.reply(status_codes::OK, answer);
}


// handle post request
void handle_post(http_request request)
{
   // TRACE("\nhandle POST\n");
   std::cout << "\nhandle POST" << std::endl;

   handle_request(
      request,
      parse_json
   );
}


// handle request by action
void handle_request(http_request request, std::function<void(json::value const &, json::value &)> action)
{
   json::value answer;

   request
      .extract_json()
      .then([&answer, &action](pplx::task<json::value> task) {
         try
         {
            auto const & jvalue = task.get();
            display_json(jvalue, "R: ");

            if (!jvalue.is_null())
            {
               action(jvalue, answer);
            }
            
            display_json(answer, "S: ");
         }
         catch (http_exception const & e)
         {
            std::cout << e.what() << std::endl;
         }
      })
      .wait();

   if (answer != json::value::null())
      request.reply(status_codes::OK, answer);
   else
      request.reply(status_codes::OK, json::value::array());
}
