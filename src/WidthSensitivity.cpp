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
#include <FisherInformation.hpp>

void Debug::init()
{
    if constexpr(Debug::ENABLE)
    {
        if(std::filesystem::exists(Debug::DIR))
            std::filesystem::remove_all(Debug::DIR);
        std::filesystem::create_directories(Debug::DIR);
        log_file.open(Debug::LOG_FILE);

        if(!log_file.is_open())
        {
            const auto error_content{ std::format("Width::Debug::init(): impossibile aprire il file di log {}", Debug::LOG_FILE) };
            throw std::runtime_error(error_content);
        }
    }
}

void Debug::save_spectrum(
    std::string const& filename,
    std::string const& title,
    Spectrum const& spectrum)
{
    if constexpr(Debug::ENABLE)
    {
        auto params{ Debug::plot_params };
        save_plot(spectrum.get_hist(), filename, params);
    }
}

void Debug::save_all_spectra(
    std::vector<Spectrum> const& spectra,
    std::vector<double> const& deltas)
{
    if constexpr(Debug::ENABLE)
    {
        for(std::size_t i{}; auto& spectrum : spectra)
        {
            auto filename{ std::format("{}/spectrum{}.png", Debug::DIR, i) };
            auto title{ std::format("#Delta #Gamma: {:.2e} GeV", deltas[i]) };
            Debug::save_spectrum(filename, title, spectrum);
            ++i;
        }
    }
}

void Debug::save_derivative(
    std::string const& filename, 
    std::string const& title, 
    std::vector<double> const& derivative,
    std::vector<double> const& uncertainties)
{
    if constexpr(Debug::ENABLE)
    {
        const std::size_t bin_count{ derivative.size() };

        auto hist{ std::make_unique<TH1D>("", "", bin_count, 0, static_cast<double>(bin_count)) };

        for (std::size_t j{}; j < bin_count; ++j)
        {
            hist->SetBinContent(j + 1, derivative[j]);
            hist->SetBinError(j + 1, uncertainties[j]);
        }
        auto params{ plot_params };
        params.name = title.c_str();
        save_plot(hist.get(), filename, params);
    }
}

void Debug::save_all_derivatives(
    std::vector<std::vector<double>> const& derivatives,
    std::vector<std::vector<double>> const& uncertainties,
    std::vector<double> const& deltas)
{
    if constexpr(Debug::ENABLE)
    {
        for(std::size_t i{}; i < derivatives.size(); ++i)
        {
            auto title{ std::format("#Delta #Gamma: {:.2e} GeV", deltas[i]) };
            auto filename{ std::format("{}/derivative{}.png", Debug::DIR, i) };
            
            Debug::save_derivative(filename, title, derivatives[i], uncertainties[i]);
        }
    }
}


void Debug::log(std::string_view content)
{
    if constexpr(Debug::ENABLE)
    {
        std::cout << content << std::endl;
        log_file << content << std::endl;
    }
}

std::vector<double> Width0::estimate_spectrum_derivative(
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

std::vector<double> Width0::estimate_spectrum_derivative_uncertainties(
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

Spectrum Width0::generate_spectrum(
    double mass,
    double width,
    std::size_t event_count)
{
    auto spectrum{ generate_spectrum_unnormalized(mass, width, event_count) };   
    spectrum.normalize();

    return spectrum;  
}

Spectrum Width0::generate_spectrum_unnormalized(
    double mass,
    double width,
    std::size_t event_count)
{
    auto pT_gen{ std::make_unique<PT_Generator>() };

    auto w_gen{ std::make_unique<W_Generator>(
        mass,
        width,
        pT_gen.get()
    )};

    WDecaySampler sampler{ w_gen.get() };

    Spectrum spectrum{sampler, event_count};

    return spectrum;
}

double Width0::estimate_range_lower_bound(
    double mass,
    double width,
    const Spectrum& nominal)
{
    auto hist_ptr = nominal.get_hist();

    double max_sigma{};

    for (int i = 1; i <= hist_ptr->GetNbinsX(); ++i)
        max_sigma = std::max(max_sigma, hist_ptr->GetBinError(i));

    return max_sigma * width / 2.;
}

double Width0::estimate_derivative_to_noise_ratio(
    std::vector<double> const& derivatives,
    std::vector<double> const& uncertainties)
{
    double dnr{};

    const auto bin_count{ derivatives.size() };
    for(std::size_t i{}; i < bin_count; ++i)
    {
        dnr += std::pow(derivatives[i] / uncertainties[i], 2);
    }
    
    return dnr / (bin_count - 1);
}

double Width0::Impl::estimate_derivative_dnr(
    Spectrum const& nominal,
    double mass,
    double width,
    double delta,
    std::size_t event_count)
{
    auto shifted{ generate_spectrum(mass, width + delta, event_count) };
    
    auto derivatives{ estimate_spectrum_derivative(nominal, shifted, delta) };
    auto uncertainties{ estimate_spectrum_derivative_uncertainties(nominal, shifted, delta) };

    auto dnr{ estimate_derivative_to_noise_ratio(derivatives, uncertainties) };

    {
        static std::size_t counter{};
        ::Debug::save_derivative(
            std::format("{}/derivative_dnr{}.png", Debug::DIR, counter++),
            std::format("#Delta#Gamma = {:.2e}, dnr = {:.2e}", delta, dnr),
            derivatives,
            uncertainties
        );
    }
    

    return dnr;
}

double Width0::Impl::estimate_delta_lower_bound(
    Spectrum const& nominal,
    double mass,
    double width,
    std::size_t event_count,
    double gamma_lower_bound,
    double gamma_upper_bound,
    double lower_dnr,
    double upper_dnr,
    std::size_t max_iteration_count,
    double dnr_lower_bound,
    double epsilon)
{
    Debug::log("estimate_delta_lower_bound()...");
    auto delta{ std::sqrt(gamma_lower_bound * gamma_upper_bound) };
    
    auto dnr{ Impl::estimate_derivative_dnr(nominal, mass, width, delta, event_count) }; 
    
    auto new_lower_bound{ gamma_lower_bound };
    auto new_upper_bound{ gamma_upper_bound };
    auto new_lower_dnr{ lower_dnr };
    auto new_upper_dnr{ upper_dnr };
    if(max_iteration_count == 0)
    {
        Debug::log("concluso per numero massimo di iterazioni...");
        Debug::log(std::format("risultato estimate_delta_lower_bound(): {:.2e}", delta));
        return delta;
    }
    if(std::abs(dnr - dnr_lower_bound) < epsilon)
    {
        Debug::log(std::format("risultato estimate_delta_lower_bound(): {:.2e}", delta));
        return delta;
    }
    else if(dnr > dnr_lower_bound)
    {
        new_upper_bound = delta;
        new_upper_dnr = dnr;
    }
    else
    {
        new_lower_bound = delta;
        new_lower_dnr = dnr;
    }
    
    Debug::log(std::format("chiamata ricorsione con estremi [{:.2e}, {:.2e}]", new_lower_bound, new_upper_bound));

    return estimate_delta_lower_bound(
        nominal,
        mass,
        width,
        event_count,
        new_lower_bound,
        new_upper_bound,
        new_lower_dnr,
        new_upper_dnr,
        max_iteration_count - 1,
        dnr_lower_bound,
        epsilon
    );
}

double Width0::estimate_delta_lower_bound(
    Spectrum const& nominal,
    double mass,
    double width,
    std::size_t event_count,
    double gamma_lower_bound,
    double gamma_upper_bound,
    std::size_t max_iteration_count,
    double dnr_lower_bound,
    double epsilon)
{
    auto lower_dnr{ 
        Impl::estimate_derivative_dnr(
            nominal, 
            mass, 
            width, 
            gamma_lower_bound, 
            event_count
        ) 
    };

    auto upper_dnr{ 
        Impl::estimate_derivative_dnr(
            nominal, 
            mass, 
            width, 
            gamma_upper_bound, 
            event_count
        ) 
    };
    Debug::log(std::format("lower_dnr: {}, upper_dnr: {}", lower_dnr, upper_dnr));

    return Impl::estimate_delta_lower_bound(
        nominal,
        mass,
        width,
        event_count,
        gamma_lower_bound,
        gamma_upper_bound,
        lower_dnr,
        upper_dnr,
        max_iteration_count,
        dnr_lower_bound,
        epsilon
    );
}

double Width0::estimate_cost_function(
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

    return cost / (bin_count - 1);
}

double Width0::estimate_delta_upper_bound(
    double gamma_lower_bound,
    double gamma_upper_bound,
    double mass,
    double width,
    std::size_t event_count,
    std::size_t coarseness,
    std::size_t point_count,
    std::size_t k,
    double m)
{
    Debug::log("estimate_delta_upper_bound()...");

    const auto coarse_event_count{ event_count / coarseness };
    const auto nominal{ generate_spectrum(mass, width, coarse_event_count) };
    const auto dof{ nominal.get_hist()->GetNbinsX() - 1 };
    const double sigma{ std::sqrt(2.0 / dof) };

    Debug::log(std::format("sigma:{}, dof:{}", sigma, dof));

    std::vector<Spectrum> spectra;
    std::vector<std::vector<double>> derivatives, uncertainties;
    std::vector<double> costs;
    std::vector<double> deltas;

    auto add_derivative = [&](std::size_t i)
    {
        spectra.push_back(generate_spectrum(mass, width + deltas[i], coarse_event_count) );
        derivatives.push_back(estimate_spectrum_derivative(nominal, spectra[i], deltas[i]));
        uncertainties.push_back(estimate_spectrum_derivative_uncertainties(nominal, spectra[i], deltas[i]));

        {
            Debug::save_derivative(
                std::format("{}/derivative_cost{}.png", Debug::DIR, i),
                std::format("#Delta#Gamma = {:.2e}, ", deltas[i]),
                derivatives[i],
                uncertainties[i]
            );
        }
    };
    auto add_cost = [&](std::size_t i)
    {
        auto cost{ estimate_cost_function(
                    derivatives[i],
                    derivatives[i + 1],
                    uncertainties[i],
                    uncertainties[i + 1],
                    nominal,
                    deltas[i],
                    deltas[i + 1]
                )
            };
        costs.push_back(cost);
        
        Debug::log(std::format("{}: cost:{:.2e}, delta:{:.2e}", i, costs[i], deltas[i]));
    };

    double ratio{ gamma_lower_bound / gamma_upper_bound };
    for(std::size_t i{}; i  < point_count; ++i)
    {
        auto delta{ gamma_upper_bound * std::pow(ratio, static_cast<double>(i) / (point_count - 1)) };
        deltas.push_back(delta);
    }

    
    const auto confidency_upper_bound{ 1.0 + sigma * k };
    double delta{ deltas[point_count - 2] * deltas[point_count - 1] };

    add_derivative(0);
    for(std::size_t i{}; i < point_count - 1; ++i)
    {
        std::size_t j{};
        bool found{ false };
        while(i + j < point_count - 1)
        {
            add_derivative(i + j + 1);
            add_cost(i + j);

            if(costs[i + j] > confidency_upper_bound)
                break;

            if(j == m - 1)
            {
                delta = std::sqrt(deltas[i] * deltas[i + 1]);
                Debug::log(std::format("trovato... posizione {}", i));
                found = true;
                break;
            }

            j++;
        }
        i += j;
        if(found)
            break;
    }

    Debug::log(std::format("risultato estimate_delta_upper_bound(): {:.2e}", delta / std::sqrt(coarseness)));

    return delta / std::sqrt(coarseness);
}

double Width0::find_best_delta(
    Spectrum const& nominal,
    double mass,
    double width,
    std::size_t event_count)
{
    auto range_lower_bound{ estimate_range_lower_bound(mass, width, nominal) };
    auto range_upper_bound{ width };

    auto delta_upper_bound{
        estimate_delta_upper_bound(
            range_lower_bound,
            range_upper_bound,
            mass,
            width,
            event_count,
            MONTECARLO_COARSENESS,
            POINT_COUNT,
            CONFIDENCE_SIGMA_MULTIPLIER,
            STABILITY_WINDOW_SIZE
        )
    };

    auto delta_lower_bound{ 
        estimate_delta_lower_bound(
            nominal,
            mass,
            width,
            event_count,
            range_lower_bound,
            delta_upper_bound,
            MAX_ITERATION_COUNT,
            DNR_LOWER_BOUND,
            EPSILON
        )
    };

    if(delta_lower_bound > delta_upper_bound)
    {
        throw std::runtime_error{ "Width::find_best_delta(): i limiti trovati non sono compatibili! delta_lower_bound > delta_upper_bound" };
    }

    auto final_delta{ std::sqrt(delta_lower_bound * delta_upper_bound) };
    Debug::log(std::format("-----> delta: {}", final_delta));
    return final_delta;
}

double Width0::estimate_sigma(
    double mass,
    double width,
    std::size_t event_count)
{
    Debug::init();

    auto nominal{ generate_spectrum_unnormalized(mass, width, event_count) };
    auto normalized_nominal{ nominal };
    normalized_nominal.normalize();

    auto delta{ find_best_delta(normalized_nominal, mass, width, event_count) };

    auto shifted{ generate_spectrum_unnormalized(mass, width + delta, event_count) };

    FisherInformation fi{ nominal.get_hist(), shifted.get_hist(), delta };

    auto sigma{ fi.sigma() };
    Debug::log(std::format("incertezza sulla larghezza di decadimento: {}", sigma));

    return sigma;
}