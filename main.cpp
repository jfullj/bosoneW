#include <TApplication.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <SpectrumBuilder.hpp>
#include <MassSensitivityAnalyzer.hpp>
#include <HistUtils.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <chrono>
#include <sstream>


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

const std::size_t CANVAS_WIDTH = 1280;
const std::size_t CANVAS_HEIGHT = 720;

const char* output_dir{ DATA_DIR "/results" };
const char* input_file{ DATA_DIR "/input/parameters.txt" };
const char* output_file{ DATA_DIR "/results/sigma_mass.txt" };
const char* output_hist_file{ DATA_DIR "/results/template_comparison_plot.png" };

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
    double W_DELTA = params.W_MASS1 - params.W_MASS0;



    auto pT_gen{ std::make_unique<PT_Generator>() };


    WDecaySampler sampler0{
        params.W_MASS0,
        params.W_WIDTH,
        pT_gen.get()
    };
    WDecaySampler sampler1{
        params.W_MASS1,
        params.W_WIDTH,
        pT_gen.get()
    };

    auto pdf0{ SpectrumBuilder{sampler0, EVENT_COUNT}.releaseHist() };
    auto pdf1{ SpectrumBuilder{sampler1, EVENT_COUNT}.releaseHist() };

    MassSensitivityAnalyzer analyzer{
        pdf0.get(),
        pdf1.get(),
        W_DELTA
    };

    TemplateComparison tc(pdf0.get(), pdf1.get(), {
        .width = CANVAS_WIDTH,
        .height = CANVAS_HEIGHT,
        .nominal_mass = params.W_MASS0,
        .shifted_mass = params.W_MASS1,
        .pad_height_ratio = 2,
        .gap = 0.02,
        .nominal_hist_color = kBlue + 1,
        .shifted_hist_color = kRed + 1,
        .upper_x_axis_content = "",
        .lower_x_axis_content = "p_{T}^{#mu} [GeV]",
        .upper_y_axis_content = "pdf",
        .lower_y_axis_content = "ratio",
        .title_size = 0.055,
        .label_size = 0.045,
        .title_offset_x = 1.20,
        .title_offset_y = 0.70
    });

    tc.save_as(output_hist_file);

    auto ratioHist{ analyzer.releaseRatioHist() };
    double sigma{ analyzer.sigma() };

    auto end{ std::chrono::high_resolution_clock::now() };
    auto elapsed{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) };

    std::ostringstream content;
    content << "sigma_mass = " << sigma << " GeV\n\n"
            << "generated events = " << EVENT_COUNT << "\n"
            << "selected events = " << analyzer.get_selected_events_count() << "\n"
            << "acceptance mass 0 = " << (static_cast<double>(pdf0->Integral()) / EVENT_COUNT) << "\n"
            << "acceptance mass 1 = " << (static_cast<double>(pdf1->Integral()) / EVENT_COUNT) << "\n"
            << "execution time = " << elapsed.count() / 1000. << " s\n";


    std::ofstream output{ output_file };
    if (!output)
        throw std::runtime_error{ "Impossibile creare sigma_mass.txt" };

    output << content.str();
    output.close();

    std::cout << content.str();
    return 0;
}