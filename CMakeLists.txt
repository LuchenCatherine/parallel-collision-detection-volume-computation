cmake_minimum_required(VERSION 3.9)
project(HTTP_SERVICE)

set(CMAKE_BUILD_TYPE "Release")

find_package(CGAL REQUIRED)
find_package(cpprestsdk REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

add_executable(server2 server2.cpp service.cpp)

file(GLOB GEOMETRY_SRC_FILES geometry/*.cpp)
add_library(geometry ${GEOMETRY_SRC_FILES})
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry)
target_link_libraries(geometry PUBLIC Eigen3::Eigen)
target_link_libraries(geometry PUBLIC CGAL::CGAL)
# target_link_libraries(geometry PUBLIC ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

target_link_libraries(geometry PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(server2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(server2 Eigen3::Eigen)
target_link_libraries(server2 OpenMP::OpenMP_CXX)
target_link_libraries(server2 CGAL::CGAL)
target_link_libraries(server2 geometry)
target_link_libraries(server2 cpprestsdk::cpprest)
# add_executable(server server.cpp)
# target_link_libraries(server PRIVATE cpprestsdk::cpprest)
# target_link_libraries (server Eigen3::Eigen)
# target_link_libraries(server PRIVATE geometry)
## target_link_libraries(server ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

add_executable(experiment experiment.cpp)
target_link_libraries(experiment Eigen3::Eigen)
target_link_libraries(experiment OpenMP::OpenMP_CXX)
target_link_libraries(experiment CGAL::CGAL)
target_link_libraries(experiment geometry)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0 COMPONENTS filesystem)




if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    target_link_libraries(server2 ${Boost_LIBRARIES})
    target_link_libraries(experiment ${Boost_LIBRARIES})
endif()


add_executable(client client.cpp)
target_link_libraries(client PRIVATE cpprestsdk::cpprest)

