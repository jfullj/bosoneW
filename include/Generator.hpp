#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <functional>
#include <memory>
#include <random>
#include <atomic>
#include <concepts>

#include <TRandom3.h>
#include <TH1.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

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

//Singleton che genera numeri casuali. Threadsafe
//Così diventa accessibile da tutti i generatori e anche da altre 
//classi che necessitano di generare numeri
class Random
{
public:
    static double get()
    {
        thread_local auto rand = [](){
            thread_local std::mt19937 gen(std::random_device{}());
            thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);

            return dist(gen);
        };
        return rand();
    }
};

//interfaccia per tutte le classi che generano un dato
template<typename T>
class Generator
{
public:
    virtual T operator()() = 0;
    virtual std::unique_ptr<Generator> clone() const = 0;
    virtual ~Generator() {};
};
//genera pT da una distribuzione sperimentale
class PT_Generator : public Generator<double>
{
public:
    PT_Generator();

    PT_Generator& operator=(const PT_Generator&) = delete;


    PT_Generator(PT_Generator&&) = default;
    PT_Generator& operator=(PT_Generator&&) = default;

    virtual double operator()();
    TH1* get_hist();
    virtual std::unique_ptr<Generator> clone() const;
    virtual ~PT_Generator() = default;

    static const inline char* path{ DATA_DIR "/distribution_pTW.root" };
private:
    PT_Generator(const PT_Generator&);

    std::unique_ptr<TH1D> hist;
    std::unique_ptr<TRandom3> rng;
};

//genera pT da una distribuzione a delta di Dirac
class PT_Delta_Generator : public Generator<double>
{
public:
    PT_Delta_Generator() = default;
    PT_Delta_Generator(double pT);

    PT_Delta_Generator& operator=(const PT_Delta_Generator&) = delete;


    PT_Delta_Generator(PT_Delta_Generator&&) = default;
    PT_Delta_Generator& operator=(PT_Delta_Generator&&) = default;

    virtual double operator()();
    virtual std::unique_ptr<Generator> clone() const;
    virtual ~PT_Delta_Generator() = default;
private:
    PT_Delta_Generator(const PT_Delta_Generator&) = default;

    double pT;
};

class BreitWignerGenerator : public Generator<double>
{
public:
    BreitWignerGenerator() = default;
    BreitWignerGenerator(double mass, double width);

    BreitWignerGenerator& operator=(const BreitWignerGenerator&) = delete;


    BreitWignerGenerator(BreitWignerGenerator&&) = default;
    BreitWignerGenerator& operator=(BreitWignerGenerator&&) = default;

    virtual double operator()();
    virtual std::unique_ptr<Generator> clone() const;
    virtual ~BreitWignerGenerator() = default;
private:
    BreitWignerGenerator(const BreitWignerGenerator&) = default;

    double mass, width;
};

using LorentzVector = ROOT::Math::PxPyPzEVector;
//genera la massa del W, richiede in ingresso il generatore del momento trasverso.
class W_Generator : public Generator<LorentzVector>
{
public:
    W_Generator() = default;
    W_Generator(const Generator<double>* const pT, const Generator<double>* const W_mass);

    W_Generator& operator=(const W_Generator&) = delete;


    W_Generator(W_Generator&&) = default;
    W_Generator& operator=(W_Generator&&) = default;

    virtual LorentzVector operator()();
    virtual std::unique_ptr<Generator> clone() const;
    virtual ~W_Generator() = default;

    static constexpr double MIN_RAPIDITY = -3.0;
    static constexpr double MAX_RAPIDITY = 3.0;
protected:
    std::unique_ptr<Generator<double>> pT_gen, W_mass_gen;
private:
    W_Generator(const W_Generator&);
};


#endif //GENERATOR_HPP