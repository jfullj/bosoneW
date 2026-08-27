#include <TApplication.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <chrono>

#include <MassSensitivity.hpp>
#include <WidthSensitivity.hpp>

namespace fs = std::filesystem;

struct Parameters
{
    double W_MASS0;
    double W_MASS1;
    double W_WIDTH;
    std::size_t EVENT_COUNT;
    bool USE_ACCELERATOR_DATA;
    double INTEGRATED_LUMINOSITY;
    double CROSS_SECTION;
};
std::size_t get_event_count(Parameters p)
{
    if(p.USE_ACCELERATOR_DATA == true)
        return static_cast<std::size_t>(p.CROSS_SECTION * p.INTEGRATED_LUMINOSITY);
    return p.EVENT_COUNT;
}

Parameters read_parameters(const std::string& filename)
{
    Parameters p{};

    std::ifstream file(filename);

    if (!file)
        throw std::runtime_error("Impossibile aprire il file: " + filename);

    std::string parameter;

    p.USE_ACCELERATOR_DATA = true;
    while (file >> parameter)
    {
        if (parameter == "W_MASS0")
            file >> p.W_MASS0;
        else if (parameter == "W_MASS1")
            file >> p.W_MASS1;
        else if (parameter == "W_WIDTH")
            file >> p.W_WIDTH;
        else if(parameter == "EVENT_COUNT")
            file >> p.EVENT_COUNT, p.USE_ACCELERATOR_DATA = false;
        else if(parameter == "INTEGRATED_LUMINOSITY")
            file >> p.INTEGRATED_LUMINOSITY;
        else if(parameter == "CROSS_SECTION")
            file >> p.CROSS_SECTION;
        else
            throw std::runtime_error{ "Parametro sconosciuto: " + parameter };
    }

    return p;
}

const char* input_file{ DATA_DIR "/input/parameters.txt" };
const char* output_dir{ DATA_DIR "/results" };


int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);
    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT();

    auto start{ std::chrono::high_resolution_clock::now() };

    TApplication app("app", &argc, argv);

    fs::create_directories(output_dir);

    auto params { read_parameters(input_file) };
    auto EVENT_COUNT{ get_event_count(params) };    
    
    //Mass::junk(params.W_MASS0, params.W_MASS1, params.W_WIDTH, EVENT_COUNT);

    {
        auto pT_gen{ std::make_unique<PT_Generator>() };
        

        auto w_gen{ std::make_unique<W_Generator>(
            params.W_MASS0,
            params.W_WIDTH,
            pT_gen.get()
        )};
        WDecaySampler sampler{ w_gen.get() };
        Spectrum pdf{sampler, EVENT_COUNT};

        Width::find_best_delta_width(params.W_MASS0, params.W_WIDTH, EVENT_COUNT, 25, pdf);
    }

    auto end{ std::chrono::high_resolution_clock::now() };
    auto elapsed{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) };
    
    return 0;
}