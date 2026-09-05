cmake_minimum_required(VERSION 3.16)

project(bosoneW LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT REQUIRED COMPONENTS 
    Core 
    RIO 
    Hist 
    Graf 
    Gpad 
    Tree 
    TreePlayer 
    MathCore
    MathMore
    Imt 
    ROOTDataFrame 
    Physics)

file(GLOB_RECURSE SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

file(GLOB_RECURSE HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
)

add_executable(bosoneW
    main.cpp
    ${SOURCES}
    ${HEADERS}
)

target_link_libraries(bosoneW PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::MathCore
    ROOT::Imt
    ROOT::ROOTDataFrame
    ROOT::Physics
    ROOT::MathMore
)

target_compile_definitions(bosoneW PRIVATE
    DATA_DIR="${CMAKE_SOURCE_DIR}"
)

target_include_directories(bosoneW PRIVATE include)
