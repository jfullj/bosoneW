cmake_minimum_required(VERSION 3.16)

project(bosoneW LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fsanitize=address")

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

include(${ROOT_USE_FILE})

add_executable(bosoneW
    main.cpp
    src/WDecaySampler.cpp
    src/SpectrumBuilder.cpp
    src/FisherInformation.cpp
    src/HistUtils.cpp
    src/Generator.cpp
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
