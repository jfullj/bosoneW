#ifndef WDECAYSAMPLER_HPP
#define WDECAYSAMPLER_HPP

#include <Generator.hpp>
#include <thread>

namespace Event
{
    struct Type
    {
        double w_invariant_mass;
        double muon_pT;
        double muon_eta;
        double mT;
    };

    namespace Transform
    {
        inline Type standard(Type const& event)
        {
            return event;
        }
        template<typename T>
        concept Callable =
            std::invocable<T, Type const&> &&
            std::same_as<
                std::invoke_result_t<T, Type const&>,
                Type
            >;
    }
}

namespace Impl
{
    ROOT::Math::PxPyPzEVector generate_random_muon_p_rest_frame(double w_mass, double muon_mass);
    Event::Type generate_decay_event(ROOT::Math::PxPyPzEVector const& W_p, double invariant_mass, double muon_mass);
}

//Genera un evento del decadimento
class WDecaySampler
{
public:
    WDecaySampler() = delete;
    explicit WDecaySampler(const Generator<LorentzVector> * const W_generator);

    WDecaySampler(const WDecaySampler&) = default;
    WDecaySampler& operator=(const WDecaySampler&) = default;

    WDecaySampler(WDecaySampler&&) = default;
    WDecaySampler& operator=(WDecaySampler&&) = default;

    Event::Type operator()() const;

    ~WDecaySampler() = default;

    static constexpr double MUON_MASS = 0.105658; // GeV/c^2
private:
    std::vector<std::unique_ptr<Generator<LorentzVector>>> w_local_generators;
};


#endif //WDECAYSAMPLER_HPP