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
        auto [sigma0, hist0] = Width0::estimate_sigma(params.W_MASS0, params.W_WIDTH, EVENT_COUNT);
        auto [sigma1, hist1] = Width1::estimate_sigma(params.W_MASS0, params.W_WIDTH, EVENT_COUNT);

        TemplateComparison tc(hist0.get(), hist1.get(), sigma0, {
            .width = 1280,
            .height = 720,
            .nominal_mass = sigma0,
            .shifted_mass = sigma1,
            .pad_height_ratio = 2,
            .gap = 0.02,
            .nominal_hist_color = kBlue + 1,
            .shifted_hist_color = kRed + 1,
            .upper_x_axis_content = "",
            .lower_x_axis_content = "p_{T}^{#mu} [GeV]",
            .upper_y_axis_content = "events / N",
            .lower_y_axis_content = "#frac{f_{m+#Deltam}}{f_{m}}",
            .title_size = 0.055,
            .label_size = 0.045,
            .title_offset_x = 1.20,
            .title_offset_y = 1.20
        });

        tc.save_as(DATA_DIR "/results/width derivatives.png");
    }

    auto end{ std::chrono::high_resolution_clock::now() };
    auto elapsed{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) };
    
    return 0;
}