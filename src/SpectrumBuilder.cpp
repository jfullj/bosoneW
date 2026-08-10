#include "SpectrumBuilder.hpp"
#include <ROOT/RDataFrame.hxx>

SpectrumBuilder::SpectrumBuilder(const WDecaySampler& sampler, std::size_t event_count)
: event_count{ event_count }
{
    ROOT::RDataFrame df{ event_count };

    auto h_muon_pt = df.Define("muon", [&sampler](){
        return sampler();
    })
    .Define("muon_pT", [](std::pair<double, double> muon){ return muon.first; }, {"muon"})
    .Define("muon_eta", [](std::pair<double, double> muon){ return muon.second; }, {"muon"})
    .Filter([=](double muon_pT, double muon_eta){
         return (muon_pT > MIN_MUON_PT && muon_pT < MAX_MUON_PT) &&
                (muon_eta > MIN_ETA && muon_eta < MAX_ETA);}, {"muon_pT", "muon_eta"})
    .Histo1D({
        "h_pt",
        "Muon pT;p_{T}^{#mu} [GeV];Events",
        BIN_COUNT,
        MIN_MUON_PT,
        MAX_MUON_PT 
    }, "muon_pT" );

    hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(h_muon_pt->Clone("h_pt")));
}

TH1D* SpectrumBuilder::getHist() const
{
    return hist.get();
}

std::unique_ptr<TH1D> SpectrumBuilder::releaseHist()
{
    return std::move(hist);
}

std::size_t SpectrumBuilder::get_event_count() const { return event_count; }