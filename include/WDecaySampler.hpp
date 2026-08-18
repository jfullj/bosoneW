#ifndef WDECAYSAMPLER_HPP
#define WDECAYSAMPLER_HPP

#include <functional>
#include <memory>
#include <atomic>
#include <TRandom3.h>
#include <TH1.h>

class ThreadID
{
public:
    static std::size_t get() 
    {
        thread_local std::size_t local_ID{ ID.fetch_add(1) };
        return local_ID;
    }
private:
    static inline std::atomic<std::size_t> ID{ 0 };
};

class Generator
{
public:
    virtual double operator()() = 0;
    virtual std::unique_ptr<Generator> clone() const = 0;
    virtual ~Generator() {};
};

class PT_Generator : public Generator
{
public:
    PT_Generator();

    PT_Generator& operator=(const PT_Generator&) = delete;


    PT_Generator(PT_Generator&&) = default;
    PT_Generator& operator=(PT_Generator&&) = delete;

    virtual double operator()();
    virtual std::unique_ptr<Generator> clone() const;
    virtual ~PT_Generator() = default;

    static const inline char* path{ DATA_DIR "/distribution_pTW.root" };
private:
    PT_Generator(const PT_Generator&);

    std::unique_ptr<TH1D> hist;
    std::unique_ptr<TRandom3> rng;
};

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
    explicit WDecaySampler(double WMass, double WWidth, const Generator * const pT_generator);

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
    double WMass, WWidth;
    std::vector<std::unique_ptr<Generator>> pT_local_generators;
};


#endif //WDECAYSAMPLER_HPP