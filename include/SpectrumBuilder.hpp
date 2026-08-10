#ifndef SPECTRUMBUILDER_HPP
#define SPECTRUMBUILDER_HPP

#include "WDecaySampler.hpp"
#include <TH1.h>


class SpectrumBuilder
{
public:
    SpectrumBuilder() = delete;
    SpectrumBuilder(const WDecaySampler& sampler, std::size_t event_count);

    SpectrumBuilder(const SpectrumBuilder&) = delete;
    SpectrumBuilder& operator=(const SpectrumBuilder&) = delete;

    SpectrumBuilder(SpectrumBuilder&&) = default;
    SpectrumBuilder& operator=(SpectrumBuilder&&) = default;

    static constexpr double MIN_MUON_PT = 26.0; // GeV/c
    static constexpr double MAX_MUON_PT = 56.0; // GeV/c
    static constexpr double MIN_ETA = -2.4;
    static constexpr double MAX_ETA = 2.4;

    static constexpr std::size_t BIN_COUNT = 30; 
    
    TH1D* getHist() const;
    std::unique_ptr<TH1D> releaseHist();
    std::size_t get_event_count() const;

    ~SpectrumBuilder() = default;
private:
    std::size_t event_count;

    std::unique_ptr<TH1D> hist;
};

#endif //SPECTRUMBUILDER_HPP