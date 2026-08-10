#include <WDecaySampler.hpp>

#include <TRandom3.h>
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
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
double generate_random_gaussian(double mean, double sigma,Func&& rand)
{
    double u1 = rand();
    double u2 = rand();

    while (u1 <= 0.0)
        u1 = rand();

    double z = std::sqrt(-2.0 * std::log(u1))
             * std::cos(2.0 * M_PI * u2);

    return mean + sigma * z;
}

template<std::invocable Func>
double generate_random_eta(Func&& rand) {
    //TODO: definire una distribuzione effettiva

    return generate_random_gaussian(0, 2, rand);
}

template<std::invocable Func>
double generate_random_invariant_mass(double w_mass, double w_width, Func&& rand) {
    return w_mass + w_width * std::tan(M_PI * (rand() - 0.5));
}

template<std::invocable Func>
ROOT::Math::PxPyPzEVector generate_random_muon_p_rest_frame(double m_mass, double muon_mass, Func&& rand) {
    double p = m_mass / 2.0 - muon_mass * muon_mass / (2.0 * m_mass);
    
    double cos_theta = rand() * 2.0 - 1.0;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = rand() * 2.0 * M_PI;

    double px = p * sin_theta * std::cos(phi);
    double py = p * sin_theta * std::sin(phi);
    double pz = p * cos_theta;
    double E = std::sqrt(p * p + muon_mass * muon_mass);
    return ROOT::Math::PxPyPzEVector(px, py, pz, E);
}

template<std::invocable Func>
double generate_random_phi(Func&& rand) {
    return rand() * 2.0 * M_PI;
} 

WDecaySampler::Event WDecaySampler::operator()() const
{
    thread_local auto rand = [](){
        thread_local std::mt19937 gen(
            std::random_device{}()
        );
        thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);

        return dist(gen);
    };

    double pTW{ (*pT_local_generators[ThreadID::get()])() };

    //eta ancora non ben specificato
    double eta{ generate_random_eta(rand) };
    double phi{ generate_random_phi(rand)};
    double invariant_mass{ generate_random_invariant_mass(WMass, WWidth, rand) };
    auto muon_p{ generate_random_muon_p_rest_frame(invariant_mass, MUON_MASS, rand) };

    auto neutrino_p{ - muon_p };

    ROOT::Math::PtEtaPhiMVector W_p{ pTW, eta, phi, invariant_mass };

    ROOT::Math::XYZVector beta{
        W_p.Px() / W_p.E(),
        W_p.Py() / W_p.E(),
        W_p.Pz() / W_p.E(),
    };

    ROOT::Math::Boost boost(beta);
    muon_p = boost(muon_p);
    neutrino_p = boost(neutrino_p);

    double pt_mu = muon_p.Pt();
    double pt_nu = neutrino_p.Pt();

    double dphi = ROOT::Math::VectorUtil::DeltaPhi(muon_p, neutrino_p);

    double mt = std::sqrt(2.0 * pt_mu * pt_nu * (1.0 - std::cos(dphi)));

    return { pt_mu, muon_p.Eta(), mt };

}
PT_Generator::PT_Generator()
{
    TFile file(path);

    if(file.IsZombie()) {
        throw std::runtime_error("impossibile aprire il file");
    }

    TH1D* h = dynamic_cast<TH1D*>(file.Get("h_pTW"));

    if(double I{ h->Integral("width") }; I > 0)
        h->Scale(1.0 / I);
    
        if(!h) throw std::runtime_error("impossibile creare l'istogramma della pT pdf");

    hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(h->Clone()));
    hist->SetDirectory(nullptr);
    rng = std::make_unique<TRandom3>(0);
}

double PT_Generator::operator()()
{
    return hist->GetRandom(rng.get(), "width");
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
