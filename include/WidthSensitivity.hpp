#ifndef WIDTHSENSITIVITY_HPP
#define WIDTHSENSITIVITY_HPP

#include <Spectrum.hpp>
#include <HistUtils.hpp>

#include <vector>
#include <string_view>
#include <fstream>

namespace Width
{
    namespace Debug
    {
        inline const bool ENABLE{ true };
        inline const char* DIR{ DATA_DIR "/debug" };
        inline const char* LOG_FILE{ DATA_DIR "/debug/log.txt" };
        inline std::ofstream log_file;

        inline const SavePlotParams plot_params{
            .width = 1280,
            .height = 720,
            .name = "",
            .x_axis_content = "p_{T}^{#mu} [GeV]",
            .y_axis_content = "pdf",
            .color = kBlue + 1,
            .title_size = 0.055,
            .label_size = 0.045,
            .title_offset_x = 0.80,
            .title_offset_y = 0.60,
            .draw_settings = "HIST E1 SAME"
        };
        void init();
        void log(std::string_view content);
        void save_all_spectra(
            std::vector<Spectrum> const& spectra,
            std::vector<double> const& deltas
        );
        void save_all_derivatives(
            std::vector<std::vector<double>> const& derivatives,
            std::vector<std::vector<double>> const& uncertainties,
            std::vector<double> const& deltas
        );
    }


    //calcola il limite inferiore statistico sul ΔΓ
    double estimate_delta_lower_bound(
        double mass,
        double width,
        const Spectrum& nominal
    );

    //costruisce point_count istogrammi con i valori relativi di ΔΓ
    std::pair<std::vector<Spectrum>, std::vector<double>> build_spectra(
        double lower_bound,
        double upper_bound, 
        double mass, 
        double width, 
        std::size_t event_count, 
        std::size_t point_count
    );

    //calcola le derivate di ogni istogramma rispetto a quelle della pdf nominale
    std::vector<double> estimate_spectrum_derivative(
        const Spectrum& nominal,
        const Spectrum& shifted,
        double delta_width
    );

    //calcola la incertezza su ogni derivata per un certo spettro
    std::vector<double> estimate_spectrum_derivative_uncertainties(
        const Spectrum& nominal,
        const Spectrum& shifted,
        double delta_width
    );

    //calcola quanto le derivate di due pdf sono tra di loro lontane e compatibili con gli errori statistici
    double estimate_cost_function(
        const std::vector<double>& d0,
        const std::vector<double>& d1,
        const std::vector<double>& sigma0,
        const std::vector<double>& sigma1,
        Spectrum const& nominal,
        double delta0,
        double delta1
    );

    double find_best_delta_width(
        double mass,
        double width,
        std::size_t event_count,
        std::size_t point_count,
        const Spectrum& nominal
    );
}


#endif //WIDTHSENSITIVITY_HPP