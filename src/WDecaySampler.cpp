#include <WDecaySampler.hpp>

#include <TRandom3.h>
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <random>
#include <thread>

#include <TH1.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>


WDecaySampler::WDecaySampler(double WMass, double WWidth, const Generator * const pT_generator)
: WMass{ WMass }, WWidth{ WWidth } 
{
    pT_local_generators.resize(std::thread::hardware_concurrency());
    for(auto& g : pT_local_generators)
        g = pT_generator->clone();
}

template<std::invocable Func>
double generate_random_eta(double eta_min, double eta_max, Func&& rand) {
    return (rand() * (eta_max - eta_min) + eta_min);
}

template<std::invocable Func>
double generate_random_invariant_mass(double w_mass, double w_width, Func&& rand) {
    return w_mass + w_width * std::tan(M_PI * (rand() - 0.5));
}

template<std::invocable Func>
ROOT::Math::PxPyPzEVector generate_random_muon_p_rest_frame(double m_mass, double muon_mass, Func&& rand) {
    double p = m_mass / 2.0  - muon_mass * muon_mass / (2.0 * m_mass);
    double theta = rand() * M_PI;
    double phi = rand() * 2.0 * M_PI;

    double px = p * std::sin(theta) * std::cos(phi);
    double py = p * std::sin(theta) * std::sin(phi);
    double pz = p * std::cos(theta);
    double E = std::sqrt(p * p + muon_mass * muon_mass);
    return ROOT::Math::PxPyPzEVector(px, py, pz, E);
}

template<std::invocable Func>
double generate_random_phi(Func&& rand) {
    return rand() * 2.0 * M_PI;
} 

double WDecaySampler::operator()() const
{
    thread_local auto rand = [](){
        thread_local std::mt19937 gen(
            std::random_device{}()
        );
        thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);

        return dist(gen);
    };

    double pTW{ (*pT_local_generators[ThreadID::get()])() };
    double eta{ generate_random_eta(ETA_MIN, ETA_MAX, rand) };
    double phi{ generate_random_phi(rand)};
    double invariant_mass{ generate_random_invariant_mass(WMass, WWidth, rand) };
    auto muon_p{ generate_random_muon_p_rest_frame(invariant_mass, MUON_MASS, rand) };

    ROOT::Math::PtEtaPhiEVector W_p{ pTW, eta, phi, invariant_mass };

    ROOT::Math::XYZVector beta{
        W_p.Px() / W_p.E(),
        W_p.Py() / W_p.E(),
        W_p.Pz() / W_p.E(),
    };

    ROOT::Math::Boost boost(beta);
    muon_p = boost(muon_p);

    return muon_p.Pt();
}
PT_Generator::PT_Generator()
{
    TFile file(path);

    if(file.IsZombie()) {
        throw std::runtime_error("impossibile aprire il file");
    }

    TH1D* h = dynamic_cast<TH1D*>(file.Get("h_pTW"));
    if(!h) throw std::runtime_error("impossibile creare l'istogramma della pT pdf");

    hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(h->Clone()));
    hist->SetDirectory(nullptr);
    rng = std::make_unique<TRandom3>(0);
}

double PT_Generator::operator()()
{
    return hist->GetRandom(rng.get());
}

PT_Generator::PT_Generator(const PT_Generator& g)
{
    hist = std::unique_ptr<TH1D>(
        dynamic_cast<TH1D*>(g.hist->Clone())
    );
    hist->SetDirectory(nullptr);

    rng = std::make_unique<TRandom3>(*g.rng);
}

std::unique_ptr<Generator> PT_Generator::clone() const
{
    auto copy = std::unique_ptr<PT_Generator>{ new PT_Generator{ *this } };
    return copy;  
}
