#include <Generator.hpp>

#include <TH1.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>

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

std::unique_ptr<Generator<double>> PT_Generator::clone() const
{
    auto copy = std::unique_ptr<PT_Generator>{ new PT_Generator{ *this } };
    return copy;  
}
TH1* PT_Generator::get_hist()
{
    return hist.get();
}


PT_Delta_Generator::PT_Delta_Generator(double pT) : pT{pT} {}
double PT_Delta_Generator::operator()()
{
    return pT;
}

std::unique_ptr<Generator<double>> PT_Delta_Generator::clone() const
{
    return std::unique_ptr<Generator<double>>{ new PT_Delta_Generator{ pT } };
}

W_Generator::W_Generator(double mass, double width, const Generator<double>* const pT)
: mass{mass}
, width{width}
, pT_gen{ pT->clone() }
{

}

double generate_random_eta() {
    double cos_theta = Random::get() * 2.0 - 1.0;
    return std::atanh(cos_theta);
}

double generate_random_invariant_mass(double w_mass, double w_width) {
    return w_mass + w_width * std::tan(M_PI * (Random::get() - 0.5));
}

double generate_random_phi() {
    return Random::get() * 2.0 * M_PI;
} 

double generate_random_rapidity(double min_rapidity, double max_rapidity) {
    return Random::get() * (max_rapidity - min_rapidity) + min_rapidity;
} 


ROOT::Math::PxPyPzEVector calculate_boson_p(double m_mass,double pTW, double rapidity, double phi)
{
    double mT = std::sqrt(m_mass * m_mass + pTW * pTW);
    return {
        pTW * std::cos(phi),
        pTW * std::sin(phi),
        mT * std::sinh(rapidity),
        mT * std::cosh(rapidity)
    };
}

LorentzVector W_Generator::operator()()
{
    double pTW{ (*pT_gen)() };
    double invariant_mass{ generate_random_invariant_mass(mass, width) };
    double rapidity{ generate_random_rapidity(MIN_RAPIDITY, MAX_RAPIDITY) };
    double phi{ generate_random_phi() };

    return calculate_boson_p(invariant_mass, pTW, rapidity, phi);
}

std::unique_ptr<Generator<LorentzVector>> W_Generator::clone() const
{
    return std::unique_ptr<Generator<LorentzVector>>{ new W_Generator(*this) };
}