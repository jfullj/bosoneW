#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <SpectrumBuilder.hpp>
#include <MassSensitivityAnalyzer.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <chrono>


namespace fs = std::filesystem;


struct Parameters
{
    double W_MASS0;
    double W_MASS1;
    double W_WIDTH;
    std::size_t EVENT_COUNT;
};


Parameters read_parameters(const std::string& filename)
{
    Parameters p{};

    std::ifstream file(filename);

    if (!file)
        throw std::runtime_error("Impossibile aprire il file: " + filename);

    std::string parameter;

    while (file >> parameter)
    {
        if (parameter == "W_MASS0")
            file >> p.W_MASS0;
        else if (parameter == "W_MASS1")
            file >> p.W_MASS1;
        else if (parameter == "W_WIDTH")
            file >> p.W_WIDTH;
        else if(parameter == "EVENT_COUNT")
            file >> p.EVENT_COUNT;
        else
            throw std::runtime_error{ "Parametro sconosciuto: " + parameter };
    }

    return p;
}


void save_histogram(TH1* hist, const std::string& filename)
{
    if (!hist)
        throw std::runtime_error{ "Istogramma nullo" };

    TCanvas canvas{
        "canvas",
        hist->GetTitle(),
        800,
        600
    };

    hist->Draw("HIST");

    canvas.SaveAs(filename.c_str());
}

const char* input_file{ DATA_DIR "/input/parameters.txt" };
const char* output_file{ DATA_DIR "/results/sigma_mass.txt" };


int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);
    ROOT::EnableThreadSafety();

    auto start{ std::chrono::high_resolution_clock::now() };

    TApplication app("app", &argc, argv);

    fs::create_directories("results");


    auto params { read_parameters(input_file) };


    double W_DELTA =
        params.W_MASS1 - params.W_MASS0;


    ROOT::EnableImplicitMT();


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

    auto pdf0{ SpectrumBuilder{sampler0, params.EVENT_COUNT}.releaseHist() };
    auto pdf1{ SpectrumBuilder{sampler1, params.EVENT_COUNT}.releaseHist() };

    save_histogram(pdf0.get(), "results/pdf_Wmass0.png");
    save_histogram(pdf1.get(), "results/pdf_Wmass1.png");


    MassSensitivityAnalyzer analyzer{
        pdf0.get(),
        pdf1.get(),
        W_DELTA
    };

    auto ratioHist{ analyzer.releaseRatioHist() };
    double sigma{ analyzer.sigma() };

    save_histogram(ratioHist.get(), "results/template_ratio.png");

    auto end{ std::chrono::high_resolution_clock::now() };
    auto elapsed{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) };

    std::ofstream output{ output_file };
    if (!output)
        throw std::runtime_error{ "Impossibile creare sigma_mass.txt" };

    output << "sigma_mass = " << sigma << " GeV\n\n"
           << "Execution time = " << elapsed.count() / 1000. << " s\n";

    output.close();

    std::cout << "sigma_mass = " << sigma << " GeV\n\n"
              << "Execution time = " << elapsed.count() / 1000. << " s\n";

    return 0;
}