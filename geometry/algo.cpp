#include "algo.h"

std::vector<std::pair<std::string, double>> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points)
{
    auto aabbtree_tissue = tissue.get_aabb_tree();
    std::vector<std::pair<std::string, double>> result;

    std::cout << "organ size: " << organ.size() << std::endl;
    for (int i = 0; i < organ.size(); i++)
    {
        auto AS = organ[i];
        
        if (!AS.is_surface) {std::cout << "not surface: " << AS.label << std::endl; continue;}

        auto aabbtree_AS = AS.get_aabb_tree();
        if (aabbtree_AS->do_intersect(*aabbtree_tissue))
        {
            
            if (AS.is_closed) 
            {
                double percentage = AS.percentage_points_inside(points);
                result.push_back({AS.label, percentage});
            }
            else
                result.push_back({AS.label, -1.0});
        }
        else
        {

            Surface_mesh &mesh = tissue.get_raw_mesh();
            
            // the tissue block is inside the anatomical structure
            bool is_contain_1 = true;
            for (auto vd: mesh.vertices())
            {
                Point p = mesh.point(vd);
                if (!AS.point_inside(p)) 
                {
                    is_contain_1 = false; 
                    break;
                }
            }

            // the anatomical structure is wholely inside the tissue block, still use the voxel-based algorithm, can be simplified to use the volume of the anatomical structure. 
            bool is_contain_2 = true;
            Surface_mesh &AS_raw_mesh = AS.get_raw_mesh();

            for (auto vd: AS_raw_mesh.vertices())
            {
                Point p = AS_raw_mesh.point(vd);
                
                if (!tissue.point_inside(p))
                    is_contain_2 = false;
                break;
            }

            if (is_contain_1) 
                result.push_back({AS.label, 1.0});
            else if (is_contain_2)
            {
                double percentage = AS.percentage_points_inside(points);
                result.push_back({AS.label, percentage});
            }
        }
    }

    return result;
}

std::vector<std::pair<std::string, double>> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points, GeoInfo &geo_info, int resolution)
{
    auto aabbtree_tissue = tissue.get_aabb_tree();
    std::vector<std::pair<std::string, double>> result;

    std::cout << "organ size: " << organ.size() << std::endl;
    for (int i = 0; i < organ.size(); i++)
    {
        auto AS = organ[i];
        
        if (!AS.is_surface) {std::cout << "not surface: " << AS.label << std::endl; continue;}

        auto aabbtree_AS = AS.get_aabb_tree();
        if (aabbtree_AS->do_intersect(*aabbtree_tissue))
        {
            
            if (AS.is_closed) 
            {
                double t1 =  omp_get_wtime();
                double percentage = AS.percentage_points_inside(points);

                double t2 = omp_get_wtime();
                double percentage_normal = voxelization_volume_computation(geo_info, AS, resolution);
                double t3 = omp_get_wtime();
                double percentage_omp = voxelization_volume_computation_openmp(geo_info, AS, resolution);
                double t4 = omp_get_wtime();

                printf("running time %f %f %f\n", t2 - t1, t3 - t2, t4 - t3);
                result.push_back({AS.label, percentage});
                printf("coarse percentage: %f normal percentage: %f omp percentage: %f\n", percentage, percentage_normal, percentage_omp);
            }
            else
                result.push_back({AS.label, -1.0});
        }
        else
        {

            Surface_mesh &mesh = tissue.get_raw_mesh();
            
            bool is_contain = true;
            for (auto vd: mesh.vertices())
            {
                Point p = mesh.point(vd);
                if (!AS.point_inside(p)) 
                {
                    is_contain = false; 
                    break;
                }
            }

            if (is_contain) result.push_back({AS.label, 1.0});
        }
    }

    return result;
}

double voxelization_volume_computation(GeoInfo &geo_info, Mymesh &anatomical_structure, int resolution)
{

    double min_x = geo_info.min_x;
    double min_y = geo_info.min_y;
    double min_z = geo_info.min_z;
    double max_x = geo_info.max_x;
    double max_y = geo_info.max_y;
    double max_z = geo_info.max_z;

    Eigen::Matrix3d R = geo_info.R;
    Eigen::Vector3d T = geo_info.T;
    Eigen::Vector3d origin = geo_info.origin;

    // double delta_x = (max_x - min_x) / resolution, delta_y = (max_y - min_y) / resolution, delta_z = (max_z - min_z) / resolution;    
    // double center_x, center_y, center_z;

    // for (int i = 0; i < resolution; i++)
    //     for (int j = 0; j < resolution; j++)
    //         for (int k = 0; k < resolution; k++)
    //         {
    //             center_x = min_x + (i + 0.5) * delta_x;
    //             center_y = min_y + (j + 0.5) * delta_y;
    //             center_z = min_z + (k + 0.5) * delta_z; 
    //             Eigen::Vector3d vec(center_x, center_y, center_z);
    //             // Affine transformation
    //             vec = (R*vec + T)/1000.0 + origin;
    //             Point p(vec(0), vec(1), vec(2));
    //             points.push_back(p);
    //         }

    // unit in mm
    double delta = 1.0/resolution;
    long long int num_voxel_x = resolution * (max_x - min_x);
    long long int num_voxel_y = resolution * (max_y - min_y);
    long long int num_voxel_z = resolution * (max_z - min_z);
    double center_x, center_y, center_z;

    long long int sum = 0;

    for (int i = 0; i < num_voxel_x; i++)
    {
        for (int j = 0; j < num_voxel_y; j++)
        {
            for (int k = 0; k < num_voxel_z; k++)
            {
                center_x = min_x + (i + 0.5) * delta;
                center_y = min_y + (j + 0.5) * delta;
                center_z = min_z + (k + 0.5) * delta;
                Eigen::Vector3d vec(center_x, center_y, center_z);

                // affine transformation
                vec = (R*vec + T) / 1000.0 + origin;
                Point p(vec(0), vec(1), vec(2));

                if (anatomical_structure.point_inside(p))
                    sum ++;
            }
        }
    }

    double percentage = 1.0 * sum / num_voxel_x / num_voxel_y / num_voxel_z;
    return percentage;

}

double voxelization_volume_computation_openmp(GeoInfo &geo_info, Mymesh &anatomical_structure, int resolution)
{

    double min_x = geo_info.min_x;
    double min_y = geo_info.min_y;
    double min_z = geo_info.min_z;
    double max_x = geo_info.max_x;
    double max_y = geo_info.max_y;
    double max_z = geo_info.max_z;

    Eigen::Matrix3d R = geo_info.R;
    Eigen::Vector3d T = geo_info.T;
    Eigen::Vector3d origin = geo_info.origin;

    // unit in mm
    double delta = 1.0/resolution;
    long long int num_voxel_x = resolution * (max_x - min_x);
    long long int num_voxel_y = resolution * (max_y - min_y);
    long long int num_voxel_z = resolution * (max_z - min_z);
    // double center_x, center_y, center_z;

    long long int sum = 0;

    #pragma omp parallel for collapse(3) reduction(+:sum) 
    for (int i = 0; i < num_voxel_x; i++)
    {
        for (int j = 0; j < num_voxel_y; j++)
        {
            for (int k = 0; k < num_voxel_z; k++)
            {
                double center_x = min_x + (i + 0.5) * delta;
                double center_y = min_y + (j + 0.5) * delta;
                double center_z = min_z + (k + 0.5) * delta;
                Eigen::Vector3d vec(center_x, center_y, center_z);

                // affine transformation
                vec = (R*vec + T) / 1000.0 + origin;
                Point p(vec(0), vec(1), vec(2));

                if (anatomical_structure.point_inside(p))
                    sum ++;
            }
        }
    }

    double percentage = 1.0 * sum / num_voxel_x / num_voxel_y / num_voxel_z;
    return percentage;
}


std::vector<std::string> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue)
{

    auto aabbtree_tissue = tissue.get_aabb_tree();
    std::vector<std::string> result;

    auto t1 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < organ.size(); i++)
    {
        auto AS = organ[i];
        auto aabbtree_AS = AS.get_aabb_tree();

        if (aabbtree_AS->do_intersect(*aabbtree_tissue)) result.push_back(AS.label);

    }

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = t2 - t1;
    std::cout << "normal function running time: " << duration2.count() << " seconds" << std::endl;  
    return result;

}


std::vector<std::string> collision_detection_single_tissue_parallel(std::vector<Mymesh> &organ, Mymesh &tissue)
{
    auto aabbtree_tissue = tissue.get_aabb_tree();
    std::vector<std::string> result;

    // printf("Num of CPU: %d\n", omp_get_num_procs());

    auto t1 = std::chrono::high_resolution_clock::now();

    #pragma omp parallel
    {
        std::vector<std::string> result_private;

        #pragma omp for nowait
        for (int i = 0; i < organ.size(); i++)
        {
            // printf("thread num: %d\n", omp_get_thread_num());
            auto AS = organ[i];
            auto aabbtree_AS = AS.get_aabb_tree();

            if (aabbtree_AS->do_intersect(*aabbtree_tissue)) result_private.push_back(AS.label);

        }

        #pragma omp critical
        result.insert(result.end(), result_private.begin(), result_private.end());

    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = t2 - t1;
    std::cout << "parallel function running time: " << duration2.count() << " seconds" << std::endl;  

    return result;
}

std::vector<std::string> collision_detection_brute_force(std::vector<Mymesh> &organ, Mymesh &tissue)
{

    auto tissue_faces = tissue.get_faces();
    std::vector<std::string> result;

    for (int i = 0; i < organ.size(); i++)
    {
        auto &AS = organ[i];
        auto AS_faces = AS.get_faces();

        if (intersection_brute_force(*tissue_faces, *AS_faces)) result.push_back(AS.label);
    }

    return result;

}


bool intersection_brute_force(std::vector<Triangle> &faces1, std::vector<Triangle> &faces2)
{

    for (auto it1=faces1.begin(); it1 != faces1.end(); ++it1)
    {
        for (auto it2=faces2.begin(); it2 != faces2.end(); ++it2)
        {
            if (CGAL::do_intersect(*it1, *it2))
                return true;
            
        }
    }

    return false;
}


// void create_rtree(std::vector<std::vector<Mymesh>> &organ, rtree_3 &rtree)
// {
//     CGAL::Bbox_3 bb;

//     for (int i = 0; i < organ.size(); i++)
//     {
//         double min[3], max[3];

//         bb = PMP::bbox(organ[i].get_raw_mesh());
//         min[0] = bb.xmin();
//         min[1] = bb.ymin();
//         min[2] = bb.zmin();
//         max[0] = bb.xmax();
//         max[1] = bb.ymax();
//         max[2] = bb.zmax();

//         rtree.Insert(min, max, i);
//     }

// }