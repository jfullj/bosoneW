#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <SpectrumBuilder.hpp>
#include <MassSensitivityAnalyzer.hpp>

#include <TH1.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>


namespace fs = std::filesystem;


struct Parameters
{
    double W_MASS0;
    double W_MASS1;
    double W_WIDTH;
};


Parameters readParameters(const std::string& filename)
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

        else
            throw std::runtime_error(
                "Parametro sconosciuto: " + parameter
            );
    }

    return p;
}


void saveHistogram(TH1* hist, const std::string& filename)
{
    if (!hist)
        throw std::runtime_error("Istogramma nullo");


    TCanvas canvas(
        "canvas",
        hist->GetTitle(),
        800,
        600
    );

    hist->Draw("HIST");

    canvas.SaveAs(filename.c_str());
}

const char* input_file{ DATA_DIR "/input/parameters.txt" };
const char* output_file{ DATA_DIR "/results/sigma_mass.txt" };


int main(int argc, char** argv)
{
    TApplication app("app", &argc, argv);

    gROOT->SetBatch(true);

    fs::create_directories("results");


    auto params { readParameters(input_file) };


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

    auto pdf0{ SpectrumBuilder{sampler0}.releaseHist() };
    auto pdf1{ SpectrumBuilder{sampler1}.releaseHist() };

    saveHistogram(pdf0.get(), "results/pdf_Wmass0.png");
    saveHistogram(pdf1.get(), "results/pdf_Wmass1.png");


    MassSensitivityAnalyzer analyzer{
        pdf0.get(),
        pdf1.get(),
        W_DELTA
    };

    auto ratioHist{ analyzer.releaseRatioHist() };
    double sigma{ analyzer.sigma() };

    saveHistogram(ratioHist.get(), "results/template_ratio.png");

    std::ofstream output{ output_file };
    if (!output)
        throw std::runtime_error{ "Impossibile creare sigma_mass.txt" };

    output
        << "W mass sensitivity\n\n"
        << "W_MASS0 = " << params.W_MASS0 << " GeV\n"
        << "W_MASS1 = " << params.W_MASS1 << " GeV\n"
        << "W_WIDTH = " << params.W_WIDTH << " GeV\n\n"
        << "delta_mass = " << W_DELTA << " GeV\n"
        << "sigma_mass = " << sigma << " GeV\n";
    output.close();

    std::cout << "sigma_mass = " << sigma << " GeV\n";

    return 0;
}