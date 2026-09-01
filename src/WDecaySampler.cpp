#include <WDecaySampler.hpp>

#include <TRandom3.h>
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <random>
#include <thread>

ROOT::Math::PxPyPzEVector Impl::generate_random_muon_p_rest_frame(double w_mass, double muon_mass) {
    double p = w_mass / 2.0 - muon_mass * muon_mass / (2.0 * w_mass);
    
    double cos_theta = Random::get() * 2.0 - 1.0;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = Random::get() * 2.0 * M_PI;

    double px = p * sin_theta * std::cos(phi);
    double py = p * sin_theta * std::sin(phi);
    double pz = p * cos_theta;
    double E = std::sqrt(p * p + muon_mass * muon_mass);
    return ROOT::Math::PxPyPzEVector(px, py, pz, E);
}

Event::Type Impl::generate_decay_event(ROOT::Math::PxPyPzEVector const& W_p, double invariant_mass, double muon_mass)
{
    ROOT::Math::XYZVector beta{
        W_p.Px() / W_p.E(),
        W_p.Py() / W_p.E(),
        W_p.Pz() / W_p.E(),
    };

    //momenti del muone e del neutrino nel sistema di riferimento a riposo del bosone W
    auto muon_p{ generate_random_muon_p_rest_frame(invariant_mass, muon_mass) };

    ROOT::Math::PxPyPzEVector neutrino_p{ 
        -muon_p.Px(),
        -muon_p.Py(),
        -muon_p.Pz(),
        muon_p.P() //il neutrino ha massa trascurabile
    };

    //applico il boost alle due particelle per portarle nel sistema del laboratorio
    ROOT::Math::Boost boost(beta);
    muon_p = boost(muon_p);
    neutrino_p = boost(neutrino_p);


    //calcolo della massa trasversa 
    double pt_mu = muon_p.Pt();
    double pt_nu = neutrino_p.Pt();

    double dphi = ROOT::Math::VectorUtil::DeltaPhi(muon_p, neutrino_p);

    double mt = std::sqrt(2.0 * pt_mu * pt_nu * (1.0 - std::cos(dphi)));

    return { pt_mu, muon_p.Eta(), mt };
}