#ifndef WDECAYSAMPLER_HPP
#define WDECAYSAMPLER_HPP


#include <Generator.hpp>

//Genera un evento del decadimento
class WDecaySampler
{
public:
    struct Event
    {
        double muon_pT;
        double muon_eta;
        double mT;
    };

    WDecaySampler() = delete;
    explicit WDecaySampler(const Generator<LorentzVector> * const W_generator);

    WDecaySampler(const WDecaySampler&) = default;
    WDecaySampler& operator=(const WDecaySampler&) = default;

    WDecaySampler(WDecaySampler&&) = default;
    WDecaySampler& operator=(WDecaySampler&&) = default;

    Event operator()() const;

    ~WDecaySampler() = default;

    static constexpr double MUON_MASS = 0.105658; // GeV/c^2
    static constexpr double MIN_RAPIDITY = -3.0;
    static constexpr double MAX_RAPIDITY = 3.0;
private:

    std::vector<std::unique_ptr<Generator<LorentzVector>>> w_local_generators;
};


#endif //WDECAYSAMPLER_HPP