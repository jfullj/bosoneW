#include <FisherInformation.hpp>
#include <cmath>
#include <iostream>

FisherInformation::FisherInformation(TH1D * nominal, TH1D * shifted, double parameter_delta, std::size_t accepted_events)
{
    m_fisher = 0.0;

    const auto norm_nominal{ nominal->Integral("width") };
    const auto norm_shifted{ shifted->Integral("width") };

    for (int i = 1; i <= nominal->GetNbinsX(); ++i)
    {
        const auto p_nominal{ nominal->GetBinContent(i) / norm_nominal };
        const auto p_shifted{ shifted->GetBinContent(i) / norm_shifted };

        const auto d{ (p_shifted - p_nominal) / parameter_delta };

        m_fisher += d * d / p_nominal;
    }
    m_fisher *= accepted_events;

    double fisher_variance{};
    for (int i = 1; i <= nominal->GetNbinsX(); ++i)
    {
        const auto p_nominal{ nominal->GetBinContent(i) / norm_nominal };
        const auto p_shifted{ shifted->GetBinContent(i) / norm_shifted };

        const auto d{ (p_shifted - p_nominal) / parameter_delta };

        const auto sigma_p_nominal{ nominal->GetBinError(i) / norm_nominal };
        const auto sigma_p_shifted{ shifted->GetBinError(i) / norm_shifted };
        const auto sigma_d{ 
            std::sqrt(std::pow(sigma_p_nominal, 2.) + std::pow(sigma_p_shifted, 2.)) / parameter_delta };
        
        fisher_variance += 
            std::pow(2 * d / p_nominal * sigma_d, 2) 
            + std::pow(d * d / (p_nominal * p_nominal) * sigma_p_nominal, 2)
            + 4. * std::pow(d/ p_nominal, 3.) * (sigma_p_nominal * sigma_p_nominal) / parameter_delta;
    }
    fisher_variance *= accepted_events * accepted_events;
    m_fisher_uncertainty = std::sqrt(fisher_variance);

    m_sigma = 1. / std::sqrt(m_fisher);
    m_sigma_uncertainty = m_fisher_uncertainty / (2 * std::pow(m_fisher, 1.5));
}

FisherInformation::FisherInformation(TH1D * nominal, TH1D * derivative, std::size_t accepted_events)
{
    m_fisher = 0.0;
    for (int i = 1; i <= nominal->GetNbinsX(); ++i)
    {
        const auto d{ derivative->GetBinContent(i) };
        const auto p{ nominal->GetBinContent(i) };

        m_fisher += d * d / p;
    }
    m_fisher *= accepted_events;
 
    double fisher_variance{};
    for (int i = 1; i <= nominal->GetNbinsX(); ++i)
    {
        const auto d{ derivative->GetBinContent(i) };
        const auto p{ nominal->GetBinContent(i) };

        const auto sigma_d{ derivative->GetBinError(i) };
        const auto sigma_p{ nominal->GetBinError(i) };

        fisher_variance += std::pow(2 * d / p * sigma_d, 2) + std::pow(d * d / (p * p) * sigma_p, 2);
    }
    fisher_variance *= accepted_events * accepted_events;
    m_fisher_uncertainty = std::sqrt(fisher_variance);

    m_sigma = 1. / std::sqrt(m_fisher);
    m_sigma_uncertainty = m_fisher_uncertainty / (2 * std::pow(m_fisher, 1.5));
}

double FisherInformation::sigma() const
{
    return m_sigma;
}
double FisherInformation::sigma_uncertainty() const
{
    return m_sigma_uncertainty;
}
double FisherInformation::fisher() const
{
    return m_fisher;
}
double FisherInformation::fisher_uncertainty() const
{
    return m_fisher_uncertainty;
}