#include <TApplication.h>
#include <TCanvas.h>
#include <ROOT/RDataFrame.hxx>
#include <SpectrumBuilder.hpp>
#include <MassSensitivityAnalyzer.hpp>

#include <TCanvas.h>
#include <TH1.h>

#include <TCanvas.h>
#include <TH1.h>

TCanvas* drawHistogram(TH1* hist)
{
    static int counter = 0;

    if (!hist)
        throw std::runtime_error("Istogramma nullo");

    auto canvas = new TCanvas(
        Form("canvas_%d", counter),
        hist->GetTitle(),
        800,
        600
    );

    ++counter;

    hist->Draw("HIST");
    canvas->Update();

    return canvas;
}

const double W_MASS0 = 80.360;
const double W_MASS1 = 80.460;
const double W_DELTA = W_MASS1 - W_MASS0;
const double W_WIDTH = 2.0;

int main(int argc, char** argv) {

    TApplication app("app", &argc, argv);

    ROOT::EnableImplicitMT();

    auto pT_gen{ std::make_unique<PT_Generator>() };

    WDecaySampler sampler0{ W_MASS0, W_WIDTH, pT_gen.get() },
                sampler1{ W_MASS1, W_WIDTH, pT_gen.get() };
    
    auto pdf0{ SpectrumBuilder{sampler0}.releaseHist() };
    auto pdf1{ SpectrumBuilder{sampler1}.releaseHist() };

    drawHistogram(pdf0.get());
    drawHistogram(pdf1.get());
    MassSensitivityAnalyzer analyzer{ pdf0.get(), pdf1.get(), W_DELTA };
    auto ratioHist{ analyzer.releaseRatioHist() };
    auto sigma{ analyzer.sigma() };

    drawHistogram(ratioHist.get());
    std::cout << "sigma: " << sigma << "\n";

    app.Run();

    return 0;
}


/*
const std::size_t EVENT_COUNT = 1000000;

double generate_random_pTW(const TH1* h) {
    // 1. Clona l'istogramma una sola volta per thread (uso di .Clone() generico)
    thread_local auto local_hist = [h]() {
        auto copy = std::unique_ptr<TH1>(static_cast<TH1*>(h->Clone()));
        copy->SetDirectory(0); // Scollega dal registro globale ROOT
        return copy;
    }();

    // 2. Generatore casuale separato per thread (evita la data race su gRandom)
    thread_local TRandom3 rng(0);

    // 3. Passa il generatore a GetRandom
    return local_hist->GetRandom(&rng);
}

const double ETA_MIN = -2.4;
const double ETA_MAX = 2.4;
double generate_random_eta() {
    thread_local std::mt19937 gen(
        std::random_device{}()
    );

    thread_local std::uniform_real_distribution<double> dist(ETA_MIN, ETA_MAX);

    return dist(gen);
}

const double W_MASS = 80.379; // GeV/c^2
const double W_WIDTH = 2.085; // GeV/c^2
double generate_random_invariant_mass() {
    thread_local std::mt19937 gen(
        std::random_device{}()
    );
    
    thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);

    return W_MASS + W_WIDTH * std::tan(M_PI * (dist(gen) - 0.5));
}
const double MUON_MASS = 0.105658; // GeV/c^2
TLorentzVector generate_random_muon_p_rest_frame(double m_mass) {
    thread_local std::mt19937 gen(
        std::random_device{}()
    );

    thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);

    double p = m_mass / 2.0  - MUON_MASS * MUON_MASS / (2.0 * m_mass);
    double theta = dist(gen) * M_PI;
    double phi = dist(gen) * 2.0 * M_PI;

    double px = p * std::sin(theta) * std::cos(phi);
    double py = p * std::sin(theta) * std::sin(phi);
    double pz = p * std::cos(theta);
    double E = std::sqrt(p * p + MUON_MASS * MUON_MASS);
    return TLorentzVector(px, py, pz, E);
}
double generate_random_phi() {
    thread_local std::mt19937 gen(
        std::random_device{}()
    );
    
    thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);

    return dist(gen) * 2.0 * M_PI;
} 
struct Event {
    double pTW;
    double eta;
    double phi;
    double invariant_mass;
    TLorentzVector muon_p_rest_frame;
};

Event generate_random_event(TH1* h) {
    Event event;
    event.pTW = generate_random_pTW(h);
    event.eta = generate_random_eta();
    event.phi = generate_random_phi();
    event.invariant_mass = generate_random_invariant_mass();
    event.muon_p_rest_frame = generate_random_muon_p_rest_frame(event.invariant_mass);
    return event;
}
TLorentzVector calculate_W_p(const Event& event) {
    TLorentzVector W_p;
    W_p.SetPtEtaPhiM(event.pTW, event.eta, event.phi, event.invariant_mass);
    
    return W_p;
}
const double MIN_MUON_PT = 26.0; // GeV/c
const double MAX_MUON_PT = 56.0; // GeV/c
const std::size_t BIN_NUM = 30;
double calculate_muon_pT(const Event& event, const TLorentzVector& W_p) {

    TLorentzVector muon_p = event.muon_p_rest_frame;
    muon_p.Boost(W_p.BoostVector());

    return muon_p.Pt();
}

int main(int argc, char** argv) {

    TApplication app("app", &argc, argv);

    TFile distribution_pTW("distribution_pTW.root");

    if(distribution_pTW.IsZombie()) {
        std::cout << "Impossibile aprire il file\n";
        return -1;
    }

    std::vector<TH1*> histograms;

    TIter next(distribution_pTW.GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next())) {

        TObject *obj = key->ReadObj();

        if (obj->InheritsFrom("TH1")) {
            TH1 *h = dynamic_cast<TH1*>(obj);

            if (h) {
                std::cout << "Istogramma: " << h->GetName() << '\n';
                histograms.push_back(h);
            }
        }
        else if (obj->InheritsFrom("TTree")) {
            TTree *tree = dynamic_cast<TTree*>(obj);
            std::cout << "Tree: " << tree->GetName() << '\n';
        }
    }


    std::vector<TCanvas*> canvases;

    int i = 0;
    for (TH1* h : histograms) {

        TCanvas* c = new TCanvas(
            Form("canvas_%d", i),
            h->GetName(),
            800,
            600
        );

        h->Draw();

        canvases.push_back(c);
        i++;
    }
    TCanvas* c = new TCanvas("PDF", "All Histograms", 800, 600);

    //provo a creare un istogramma direttamente dal File senza iterare
    TH1* h = (TH1*)distribution_pTW.Get("h_pTW");
    if(!h) {
        std::cout << "Impossibile trovare l'istogramma h_pTW\n";
        return -1;
    }
    h->Draw();

    //provo a usare i RDataFrame per generare numeri casuali dall'istogramma precedente
    ROOT::EnableImplicitMT();


    ROOT::RDataFrame df(EVENT_COUNT);

    auto h_muon_pt = df.Define("e", [h](){
        return generate_random_event(h);
    })
    .Define("W_p", calculate_W_p, {"e"})
    .Define("muon_pT", calculate_muon_pT, {"e", "W_p"})
    .Filter([](double muon_pT){ return muon_pT > MIN_MUON_PT && muon_pT < MAX_MUON_PT;}, {"muon_pT"})
    .Histo1D({
        "h_pt",
        "Muon pT;p_{T}^{#mu} [GeV];Events",
        BIN_NUM,      // numero bin
        MIN_MUON_PT,       // minimo
        MAX_MUON_PT      // massimo
    }, "muon_pT" );

    h_muon_pt->Draw();

    app.Run();

    return 0;
}

*/
