#include <WidthSensitivity.hpp>


#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>

#include <Math/GSLIntegrator.h>
#include <Math/Functor.h>
#include <TH1.h>

#include <HistUtils.hpp>

double Width::estimate_delta_lower_bound(
    double mass,
    double width,
    const Spectrum& nominal)
{
    auto lambda = [=](double m)
    {
        const double diff = m - mass;
        const double num = diff * diff - width * width;
        const double den = std::pow(diff * diff + width * width, 2);

        return std::abs(num / den) / M_PI;
    };

    ROOT::Math::Functor1D argument(lambda);
    ROOT::Math::GSLIntegrator integrator;

    const double m1 = mass - width;
    const double m2 = mass + width;
    const double mmax = mass + 1000.0 * width;

    const double integral =
        integrator.Integral(argument, 0.0, m1) +
        integrator.Integral(argument, m1, m2) +
        integrator.Integral(argument, m2, mmax);

    auto hist_ptr = nominal.get_hist();

    double max_sigma = 0.0;

    for (int i = 1; i <= hist_ptr->GetNbinsX(); ++i)
        max_sigma = std::max(max_sigma, hist_ptr->GetBinError(i));

    std::cout << "max_sigma" << max_sigma << "\n"
              << "integral" << integral << std::endl;
    
    return max_sigma / integral;
}

std::pair<std::vector<Spectrum>, std::vector<double>> Width::build_spectra(
    double lower_bound,
    double upper_bound, 
    double mass, 
    double width, 
    std::size_t event_count, 
    std::size_t point_count)
{
    std::vector<Spectrum> spectra;
    std::vector<double> deltas;
    const double ratio{ lower_bound / upper_bound };

    auto pT_gen{ std::make_unique<PT_Generator>() };

    for(std::size_t i{}; i < point_count; ++i)
    {
        auto k{ static_cast<double>(i) / (point_count - 1) };

        //la variazione sulla larghezza viene scalata logaritmicamente
        double delta{ upper_bound * std::pow(ratio, k) };

        auto w_gen{ std::make_unique<W_Generator>(
            mass,
            width + delta,
            pT_gen.get()
        )};

        WDecaySampler sampler{ w_gen.get() };

        Spectrum spectrum{sampler, event_count};
        spectrum.normalize();

        spectra.push_back(std::move(spectrum));
        deltas.push_back(delta);
    }

    return {spectra, deltas};
}


std::vector<double> Width::estimate_spectrum_derivative(
    const Spectrum& nominal,
    const Spectrum& shifted,
    double delta_width)
{
    std::vector<double> derivatives;

    auto ptr_nominal{ nominal.get_hist() };
    auto ptr_shifted{ shifted.get_hist() };
    
    for (int i = 1; i <= ptr_nominal->GetNbinsX(); ++i) {
        derivatives.push_back(
            (ptr_shifted->GetBinContent(i) - ptr_nominal->GetBinContent(i)) / delta_width
        );
    }
    return derivatives;
}

std::vector<double> Width::estimate_spectrum_derivative_uncertainties(
    const Spectrum& nominal,
    const Spectrum& shifted,
    double delta_width)
{
    std::vector<double> uncertainties;

    auto ptr_nominal{ nominal.get_hist() };
    auto ptr_shifted{ shifted.get_hist() };

    for (int i = 1; i <= ptr_nominal->GetNbinsX(); ++i) 
    {
        auto sigma_nominal{ ptr_nominal->GetBinError(i) };
        auto sigma_shifted{ ptr_shifted->GetBinError(i) };

        uncertainties.push_back(
            (std::sqrt(sigma_nominal * sigma_nominal + sigma_shifted * sigma_shifted)) / delta_width
        );
    }
    return uncertainties;
}

double Width::estimate_cost_function(
    const std::vector<double>& d0,
    const std::vector<double>& d1,
    const std::vector<double>& sigma0,
    const std::vector<double>& sigma1,
    Spectrum const& nominal,
    double delta0,
    double delta1)
{
    const std::size_t bin_count{ d0.size() };

    double cost{};
    auto ptr_nominal{ nominal.get_hist() };
    for(std::size_t i{}; i < bin_count; ++i)
    {
        auto sigma2{ sigma0[i] * sigma0[i] + sigma1[i] * sigma1[i] - 2 * std::pow(ptr_nominal->GetBinError(i), 2) / (delta0 * delta1) };
        cost += std::pow(d0[i] - d1[i], 2) / sigma2;
    }

    return cost / bin_count;
}

double Width::find_best_delta_width(
    double mass,
    double width,
    std::size_t event_count,
    std::size_t point_count,
    const Spectrum& nominal)
{
    Debug::init();

    auto normalized_nominal{ nominal };
    normalized_nominal.normalize();

    const auto lower_bound{ width / 1000 };//estimate_delta_lower_bound(mass, width, nominal) };
    const auto upper_bound{ width  * 100};

    std::cout << "lower_bound" << lower_bound << "\n"
              << "upper_bound" << upper_bound << std::endl;

    auto [spectra, deltas]{ build_spectra(lower_bound, upper_bound, mass, width, event_count, point_count) };
    
    Debug::save_all_spectra(spectra, deltas);

    std::vector<std::vector<double>> derivative_matrix, uncertainty_matrix;
    
    const double ratio{ upper_bound / lower_bound };
    for(std::size_t i{}; i < spectra.size(); ++i)
    {
        auto& spectrum{ spectra[i] };
        auto delta{ deltas[i] };
        spectrum.normalize();

        derivative_matrix.push_back(estimate_spectrum_derivative(normalized_nominal, spectrum, delta));
        uncertainty_matrix.push_back(estimate_spectrum_derivative_uncertainties(normalized_nominal, spectrum, delta));
    }
    Debug::save_all_derivatives(derivative_matrix, uncertainty_matrix, deltas);

    std::vector<double> costs;
    for(std::size_t i{}; i < derivative_matrix.size() - 1; ++i)
    {
        auto cost{ estimate_cost_function(
            derivative_matrix[i],
            derivative_matrix[i + 1],
            uncertainty_matrix[i],
            uncertainty_matrix[i + 1],
            normalized_nominal,
            deltas[i],
            deltas[i + 1]
        ) };

        Debug::log(std::format("{}: {}, {}, {}", i, cost, deltas[i], deltas[i + 1]));

        costs.push_back(cost);
    }


    for(std::size_t i{}; i < costs.size(); ++i)
        std::cout << deltas[i] << " " << costs[i] << "\n";

    return 0;
}

void Width::Debug::init()
{
    if constexpr(Debug::ENABLE)
    {
        std::filesystem::create_directories(Debug::DIR);
        log_file.open(Debug::LOG_FILE);

        if(!log_file.is_open())
        {
            const auto error_content{ std::format("Width::Debug::init(): impossibile aprire il file di log {}", Debug::LOG_FILE) };
            throw std::runtime_error(error_content);
        }
    }
}

void Width::Debug::save_all_spectra(
    std::vector<Spectrum> const& spectra,
    std::vector<double> const& deltas)
{
    if constexpr(Debug::ENABLE)
    {
        for(std::size_t i{}; auto& spectrum : spectra)
        {
            auto params{ Debug::plot_params };
            params.name = std::format("#Delta #Gamma: {:.2e} GeV", deltas[i]).c_str();
            auto filename{ std::format("{}/spectrum{}.png", Debug::DIR, i) };
            save_plot(spectrum.get_hist(), filename.c_str(), params);
            ++i;
        }
    }
}

void Width::Debug::save_all_derivatives(
    std::vector<std::vector<double>> const& derivatives,
    std::vector<std::vector<double>> const& uncertainties,
    std::vector<double> const& deltas)
{
    if constexpr(Debug::ENABLE)
    {
        const std::size_t bin_count{ derivatives[0].size() };

        auto hist{ std::make_unique<TH1D>("", "", bin_count, 0, static_cast<double>(bin_count)) };
        for(std::size_t i{}; i < derivatives.size(); ++i)
        {
            for (std::size_t j{}; j < bin_count; ++j)
            {
                hist->SetBinContent(j + 1, derivatives[i][j]);
                hist->SetBinError(j + 1, uncertainties[i][j]);
            }
            auto params{ plot_params };
            params.name = std::format("#Delta #Gamma: {:.2e} GeV", deltas[i]).c_str();
            auto filename{ std::format("{}/derivative{}.png", Debug::DIR, i) };
            save_plot(hist.get(), filename.c_str(), params);
        }
    }
}


void Width::Debug::log(std::string_view content)
{
    if constexpr(Debug::ENABLE)
    {
        log_file << content << "\n";
    }
}

