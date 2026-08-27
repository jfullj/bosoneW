#ifndef SPECTRUM_HPP
#define SPECTRUM_HPP

#include <WDecaySampler.hpp>
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>


namespace Acceptance
{
    inline constexpr double MIN_MUON_PT = 26.0; // GeV/c
    inline constexpr double MAX_MUON_PT = 56.0; // GeV/c
    inline constexpr double MIN_ETA = -2.4;
    inline constexpr double MAX_ETA = 2.4;
    inline constexpr double MIN_TRANSVERSE_MASS = 40; //GeV/x

    inline bool standard(WDecaySampler::Event const& e)
    {
        return (e.muon_pT > MIN_MUON_PT && e.muon_pT < MAX_MUON_PT) &&
               (e.muon_eta > MIN_ETA && e.muon_eta < MAX_ETA) &&
                e.mT > MIN_TRANSVERSE_MASS;
    }

    inline bool all(WDecaySampler::Event const& e)
    {
        return true;
    }

    template<typename T>
    concept Type =
        std::invocable<T, const WDecaySampler::Event&> &&
        std::convertible_to<
            std::invoke_result_t<T, const WDecaySampler::Event&>,
            bool
        >;
}

namespace Binning
{
    inline constexpr std::size_t BIN_COUNT = 30;

    struct Parameters
    {
        int bin_count; 
        double min,
            max;
    };

    inline constexpr Parameters standard{
        .bin_count = BIN_COUNT,
        .min = Acceptance::MIN_MUON_PT,
        .max = Acceptance::MAX_MUON_PT
    };
}

class Spectrum
{
public:
    Spectrum() = delete;
    template<Acceptance::Type Acc>
    Spectrum(const WDecaySampler& sampler, std::size_t event_count, Binning::Parameters const& params, Acc&& acceptance)
    {
        ROOT::RDataFrame df{ event_count };
        auto h_muon_pt = df.Define("muon", [&sampler](){
            return sampler();
        })
        .Define("muon_pT", [](WDecaySampler::Event muon){ return muon.muon_pT; }, {"muon"})
        .Filter(acceptance, {"muon"})
        .Histo1D({
            "", "",
            params.bin_count,
            params.min,
            params.max 
        }, "muon_pT" );

        hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(h_muon_pt->Clone("h_pt")));

        double N = hist->Integral();
        for (int i = 1; i <= hist->GetNbinsX(); ++i)
        {
            double content = hist->GetBinContent(i);
            double p = content / N;
            double error   = std::sqrt(N * p * (1 - p));

            hist->SetBinError(i, error);
        }
    }

    Spectrum(const WDecaySampler& sampler, std::size_t event_count)
    : Spectrum(sampler, event_count, Binning::standard) {}
    Spectrum(const WDecaySampler& sampler, std::size_t event_count, Binning::Parameters const& params)
    : Spectrum(sampler, event_count, params, Acceptance::standard) {}

    Spectrum(const Spectrum& other) : hist(std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(other.hist->Clone()))){}

    Spectrum& operator=(const Spectrum& other)
    {
        if (this != &other)
        {
            hist.reset(dynamic_cast<TH1D*>(other.hist->Clone()));
        }

        return *this;
    }

    Spectrum(Spectrum&&) = default;
    Spectrum& operator=(Spectrum&&) = default;    
    
    ~Spectrum() = default;

    TH1D* get_hist() const { return hist.get(); }
    std::unique_ptr<TH1D> release_hist() { return std::move(hist); }
    
    void normalize() { hist->Scale(1. / hist->Integral("width")); }
private:

    std::unique_ptr<TH1D> hist;
};

#endif //SPECTRUM_HPP