#ifndef WDECAYSAMPLER_HPP
#define WDECAYSAMPLER_HPP

#include <Generator.hpp>
#include <thread>

namespace Event
{
    struct Type
    {
        double muon_pT;
        double muon_eta;
        double mT;
    };

    struct DefaultTransformation
    {
        Type operator()(Type const& event, double invariant_mass) const
        {
            return event;
        }
    };

    template<typename T>
    concept Transformation =
        std::invocable<T, Type const&, double> &&
        std::same_as<
            std::invoke_result_t<T, Type const&, double>,
            Type
        >;
}



namespace Impl
{
    ROOT::Math::PxPyPzEVector generate_random_muon_p_rest_frame(double w_mass, double muon_mass);
    Event::Type generate_decay_event(ROOT::Math::PxPyPzEVector const& W_p, double invariant_mass, double muon_mass);
}

//Genera un evento del decadimento
template<Event::Transformation Transformation = Event::DefaultTransformation> 
class WDecaySampler
{
public:
    WDecaySampler() = delete;
    explicit WDecaySampler(
        const Generator<LorentzVector> * const W_generator,
        Transformation t = Event::DefaultTransformation{})
    : transform{std::move(t)}
    {
        for(std::size_t i{}; i < std::thread::hardware_concurrency(); ++i )
            w_local_generators.push_back(W_generator->clone());
    }

    WDecaySampler(const WDecaySampler&) = default;
    WDecaySampler& operator=(const WDecaySampler&) = default;

    WDecaySampler(WDecaySampler&&) = default;
    WDecaySampler& operator=(WDecaySampler&&) = default;

    Event::Type operator()() const
    {
        auto W_p{ (*w_local_generators[ThreadID::get()])() };
        auto invariant_mass{ W_p.M() };

        auto e{ Impl::generate_decay_event(W_p, invariant_mass, MUON_MASS) };
        return transform(e, invariant_mass);
    }

    ~WDecaySampler() = default;

    static constexpr double MUON_MASS = 0.105658; // GeV/c^2
private:
    Transformation transform;
    std::vector<std::unique_ptr<Generator<LorentzVector>>> w_local_generators;
};


#endif //WDECAYSAMPLER_HPP