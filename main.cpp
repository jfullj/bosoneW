#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <TStyle.h>
#include <TPaveText.h>

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


const std::size_t HIST_WIDTH = 1280;
const std::size_t HIST_HEIGHT = 720;

void save_histogram(TH1* hist, const std::string& filename, bool norm_hist = true)
{
    if (!hist)
        throw std::runtime_error{"Istogramma nullo"};

    TCanvas canvas{
        "canvas",
        hist->GetTitle(),
        HIST_WIDTH,
        HIST_HEIGHT
    };

    hist->SetStats(kFALSE);
    hist->Draw("HIST");

    TPaveText stats(
        0.65, 0.72,
        0.88, 0.88,
        "NDC"
    );

    stats.SetFillColor(0);
    stats.SetBorderSize(1);
    stats.SetTextAlign(12);
    stats.SetTextFont(42);
    stats.SetTextSize(0.025);

    stats.AddText(Form("Entries    %.6g", hist->GetEntries()));
    stats.AddText(Form("Mean       %.4g", hist->GetMean()));
    stats.AddText(Form("Std Dev    %.4g", hist->GetStdDev()));
    stats.AddText(Form("Integral   %.6g", hist->Integral("width")));

    stats.Draw();

    canvas.SaveAs(filename.c_str());

    if(!norm_hist)
        return;
    // Istogramma normalizzato

    std::unique_ptr<TH1> hist_norm{
        dynamic_cast<TH1*>(hist->Clone(
            (std::string{hist->GetName()} + "_norm").c_str()
        ))
    };

    if (!hist_norm)
        throw std::runtime_error{"Impossibile clonare l'istogramma"};

    hist_norm->SetDirectory(nullptr);
    hist_norm->SetStats(kFALSE);

    const double integral = hist_norm->Integral("width");

    if (integral > 0.0)
        hist_norm->Scale(1.0 / integral);

    TCanvas canvas_norm{
        "canvas_norm",
        hist_norm->GetTitle(),
        HIST_WIDTH,
        HIST_HEIGHT
    };

    hist_norm->Draw("HIST");

    TPaveText stats_norm(
        0.65, 0.72,
        0.88, 0.88,
        "NDC"
    );

    stats_norm.SetFillColor(0);
    stats_norm.SetBorderSize(1);
    stats_norm.SetTextAlign(12);
    stats_norm.SetTextFont(42);
    stats_norm.SetTextSize(0.025);

    stats_norm.AddText(Form("Entries    %.6g", hist_norm->GetEntries()));
    stats_norm.AddText(Form("Mean       %.4g", hist_norm->GetMean()));
    stats_norm.AddText(Form("Std Dev    %.4g", hist_norm->GetStdDev()));
    stats_norm.AddText(Form("Integral   %.6g", hist_norm->Integral("width")));

    stats_norm.Draw();

    // Costruisce filename_norm sostituendo l'estensione
    std::string filename_norm = filename;

    const std::size_t dot = filename_norm.find_last_of('.');

    if (dot != std::string::npos)
        filename_norm.insert(dot, "_norm");
    else
        filename_norm += "_norm";

    canvas_norm.SaveAs(filename_norm.c_str());
}

const char* input_file{ DATA_DIR "/input/parameters.txt" };
const char* output_file{ DATA_DIR "/results/sigma_mass.txt" };


int main(int argc, char** argv)
{
    gStyle->SetOptStat("eimr");
    gROOT->SetBatch(kTRUE);
    ROOT::EnableThreadSafety();

    auto start{ std::chrono::high_resolution_clock::now() };

    TApplication app("app", &argc, argv);

    fs::create_directories(DATA_DIR "/results");


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

    save_histogram(pdf0.get(), DATA_DIR "/results/pdf_Wmass0.png");
    save_histogram(pdf1.get(), DATA_DIR "/results/pdf_Wmass1.png");


    MassSensitivityAnalyzer analyzer{
        pdf0.get(),
        pdf1.get(),
        W_DELTA
    };

    auto ratioHist{ analyzer.releaseRatioHist() };
    double sigma{ analyzer.sigma() };

    save_histogram(ratioHist.get(), DATA_DIR "/results/template_ratio.png", false);

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