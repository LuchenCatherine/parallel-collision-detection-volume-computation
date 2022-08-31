# for parallel-collision-detection-volume-computation

import json
import time

import requests
import csv
from urllib.request import urlopen
import collections

server_url = 'http://localhost:12345/get-collisions'
tissue_url = 'https://ccf-api.hubmapconsortium.org/v1/hubmap/rui_locations.jsonld'


def post_request(url, rui_location, algorithm='production', resolution=1, compute_volume=True):

    tissue_id = rui_location['@id']
    tissue_id = tissue_id.split('/')[-1]
    placement = rui_location['placement']

    x_dimension = rui_location['x_dimension']
    y_dimension = rui_location['y_dimension']
    z_dimension = rui_location['z_dimension']

    target = placement['target']
    x_rotation = placement['x_rotation']
    y_rotation = placement['y_rotation']
    z_rotation = placement['z_rotation']
    x_translation = placement['x_translation']
    y_translation = placement['y_translation']
    z_translation = placement['z_translation']
    x_scaling = placement['x_scaling']
    y_scaling = placement['y_scaling']
    z_scaling = placement['z_scaling']

    json_dic = {'@id': tissue_id, 'x_dimension': x_dimension, 'y_dimension': y_dimension, 'z_dimension': z_dimension,
                'x_rotation': x_rotation, 'y_rotation': y_rotation, 'z_rotation': z_rotation,
                'x_translation': x_translation, 'y_translation': y_translation,
                'z_translation': z_translation,
                'x_scaling': x_scaling, 'y_scaling': y_scaling, 'z_scaling': z_scaling,
                'target': target, 'algorithm': algorithm, 'resolution': resolution, 'compute_volume': compute_volume}

    # print("volume: {}".format(x_dimension * y_dimension * z_dimension))
    r = requests.post(url=url, json=json_dic)
        # print(rui_location['@id'])
    return r.json()


def process_all_tissues(tissue_url, algorithm, resolution, compute_volume):
    tissue_set = set()
    total_computation_time = 0
    count = 0
    volumes = []
    with open("./rui_locations.jsonld", "r") as f:
        data = json.load(f)

        graph = data['@graph']

        for person in graph:
            samples = person['samples']
            for sample in samples:
                if sample['sample_type'] == 'Tissue Block':
                    rui_location = sample['rui_location']
                    # t1 = time.time()
                    target = rui_location['placement']['target']
                    if "#VHFLeftKidney" in target or "#VHFRightKidney" in target or "#VHMLeftKidney" in target or "#VHMRightKidney" in target:

                        tissue_id = rui_location['@id']
                        tissue_id = tissue_id.split('/')[-1]
                        # print(tissue_id)
                        tissue_set.add(tissue_id)
                        volume = rui_location['x_dimension'] * rui_location['y_dimension'] * rui_location['z_dimension']
                        # print("id: {}, volume: {}".format(tissue_id, volume))
                        volumes.append(volume)
                        r = post_request(server_url, rui_location, algorithm, resolution, compute_volume)
                        if "computation_time" in r:
                            total_computation_time += r["computation_time"]

    print(sum(volumes)/len(volumes))
    return total_computation_time


# algorithm can be "parallel" "normal" "boolean_operation"
def test_resolution(compute_volume=True):
    max_resolution = 10

    for algorithm in ["parallel", "normal", "boolean_operation"]:
        result = []
        for resolution in range(1, max_resolution):
            total_computation_time = process_all_tissues(tissue_url, algorithm, resolution, compute_volume)
            result.append(total_computation_time)

        print("{} : {}".format(algorithm, result))


# normal means only aabb
def test_index():
    for algorithm in ["brute_force", "rtree", "rtree_aabb", "normal"]:
        total_computation_time = process_all_tissues(tissue_url, algorithm, 1, False)
        print("{} costs {} to test all tissues".format(algorithm, total_computation_time))


def test_accuracy():

    max_resolution = 10

    errors = [[] for i in range(max_resolution)]
    with open("./rui_locations.jsonld", "r") as f:
        data = json.load(f)

        graph = data['@graph']
        count = 0

        for person in graph:
            samples = person['samples']
            for sample in samples:
                if sample['sample_type'] == 'Tissue Block':
                    rui_location = sample['rui_location']
                    target = rui_location['placement']['target']
                    if "#VHFLeftKidney" in target or "#VHFRightKidney" in target or "#VHMLeftKidney" in target or "#VHMRightKidney" in target:
                        count += 1
                        r2 = post_request(server_url, rui_location, "boolean_operation", 1, True)
                        for k in range(1, max_resolution):
                            r1 = post_request(server_url, rui_location, "normal", k, True)
                            if "mesh_collision_detection_result" in r2:
                                cd1 = r1["mesh_collision_detection_result"]
                                cd2 = r2["mesh_collision_detection_result"]
                                dic1 = {}
                                dic2 = {}

                                if cd2:
                                    print(cd2)
                                    for AS in cd1:
                                        dic1[AS["node_name"]] = AS["volume"]

                                    for AS in cd2:
                                        dic2[AS["node_name"]] = AS["volume"]

                                    for node in dic2:
                                        # if "renal" in node:
                                            if dic1[node] <= 0.001:
                                                error = abs(dic1[node] - dic2[node])
                                            else:
                                                error = abs(dic1[node] - dic2[node])/dic2[node]
                                            errors[k].append(error)


    print(count)
    for i in range(1, max_resolution):
        print(errors[i])
        print(len(errors[i]))
    return [sum(errors[i])/len(errors[i]) for i in range(1, max_resolution)]


if __name__ == "__main__":
    process_all_tissues(tissue_url, 'normal', 9, True)
    test_index()
    # test_resolution()
    # res = test_accuracy()
    # print(res)




