#ifndef WIDTHSENSITIVITY_HPP
#define WIDTHSENSITIVITY_HPP

#include <Spectrum.hpp>
#include <HistUtils.hpp>

#include <vector>
#include <string_view>
#include <fstream>

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

    void save_histogram(
        std::string const& filename,
        std::string const& title,
        TH1D* histogram
    );

    void save_spectrum(
        std::string const& filename,
        std::string const& title,
        Spectrum const& spectrum
    );

    void save_all_spectra(
        std::vector<Spectrum> const& spectra,
        std::vector<double> const& deltas
    );

    void save_derivative(
        std::string const& filename, 
        std::string const& title,
        std::vector<double> const& derivative,
        std::vector<double> const& uncertainties
    );

    void save_all_derivatives(
        std::vector<std::vector<double>> const& derivatives,
        std::vector<std::vector<double>> const& uncertainties,
        std::vector<double> const& deltas
    );
}

namespace Width0
{
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

    double estimate_derivative_to_noise_ratio(
        std::vector<double> const& derivatives,
        std::vector<double> const& uncertainties
    );

    double estimate_range_lower_bound(
        double mass,
        double width,
        const Spectrum& nominal
    );

    Spectrum generate_spectrum(
        double mass,
        double width,
        std::size_t event_count
    );

    Spectrum generate_spectrum_unnormalized(
        double mass,
        double width,
        std::size_t event_count
    );

    inline const double EPSILON{ 2 };
    inline const double DNR_LOWER_BOUND{ 12 };
    inline const double MAX_ITERATION_COUNT{ 10 };

    namespace Impl
    {
        double estimate_derivative_dnr(
            Spectrum const& nominal,
            double mass,
            double width,
            double delta,
            std::size_t event_count
        );
        double estimate_delta_lower_bound(
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
            double epsilon
        );
    }

    //effettua una binary search per individuare la posizione, con un grado di tolleranza
    //epsilon del punto il cui dnr è prossimo a dnr_lower_bound. La binary search si conclude
    //dopo max_iteration_count iterazioni se non è stato soddisfatto il criterio di arresto
    //in epsilon: |dnr - dnr_lower_bound| < epsilon.
    double estimate_delta_lower_bound(
        Spectrum const& nominal,
        double mass,
        double width,
        std::size_t event_count,
        double gamma_lower_bound,
        double gamma_upper_bound,
        std::size_t max_iteration_count,
        double dnr_lower_bound,
        double epsilon
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

    //suddivide l'intervallo in una griglia di point_count punti distanziati logaritmicamente
    //si cerca il primo punto per cui la funzione di costo è minore di 1 + simga * k e lo sono
    //anche i m - 1 punti successivi, in modo da garantire maggiore stabilità.
    //per ogni punto viene generato un istogramma con event_count / coarseness punti.
    inline const std::size_t MONTECARLO_COARSENESS = 100;
    inline const std::size_t STABILITY_WINDOW_SIZE = 5;
    inline const double CONFIDENCE_SIGMA_MULTIPLIER = 5;
    inline const std::size_t POINT_COUNT = 30;
    double estimate_delta_upper_bound(
        double gamma_lower_bound,
        double gamma_upper_bound,
        double mass,
        double width,
        std::size_t event_count,
        std::size_t coarseness,
        std::size_t point_count,
        std::size_t k,
        double m
    );

    double find_best_delta(
        Spectrum const& nominal,
        double mass,
        double width,
        std::size_t event_count
    );

    double estimate_sigma(
        double mass,
        double width,
        std::size_t event_count
    );
}

namespace Width1
{
    struct Transformation
    {
        double w_mass;
        double w_width;

        explicit Transformation(double w_mass, double w_width)
        : w_mass{w_mass}, w_width{w_width} {}

        Event::Type operator()(Event::Type const& e)
        {
            auto modified_e{ e };

            const double offset{ Acceptance::MUON_PT_RANGE };
            const double lb{ w_mass - w_width };
            const double ub{ w_mass + w_width };

            if(e.w_invariant_mass >= lb && e.w_invariant_mass <= ub)
                modified_e.muon_pT += offset;
            
            return modified_e;
        }
    };

    inline const Binning::Parameters BIN_PARAMS{
        .bin_count = Binning::BIN_COUNT * 2,
        .min = Acceptance::MIN_MUON_PT,
        .max = Acceptance::MIN_MUON_PT + Acceptance::MUON_PT_RANGE * 2
    };

    class BW_WidthDerivativeGenerator : public Generator<double>
    {
    public:
        BW_WidthDerivativeGenerator() = default;
        BW_WidthDerivativeGenerator(double mass, double width);

        BW_WidthDerivativeGenerator& operator=(const BW_WidthDerivativeGenerator&) = delete;


        BW_WidthDerivativeGenerator(BW_WidthDerivativeGenerator&&) = default;
        BW_WidthDerivativeGenerator& operator=(BW_WidthDerivativeGenerator&&) = default;

        virtual double operator()();
        virtual std::unique_ptr<Generator> clone() const;
        virtual ~BW_WidthDerivativeGenerator() = default;
    private:
        BW_WidthDerivativeGenerator(const BW_WidthDerivativeGenerator&) = default;
        double inverse_cumulative_function(double x) const;
        double mass, width;
    };

    Spectrum build_nominal_spectrum(
        double mass,
        double width,
        std::size_t event_count
    );

    std::unique_ptr<TH1D> make_merged_derivative(Spectrum const& spectrum);
    double calculate_normalization_constant(double mass, double width, Spectrum const& spectrum);
    std::unique_ptr<TH1D> build_derivative_spectrum(
        double mass,
        double width,
        std::size_t event_count
    );

    double estimate_sigma(
        double mass,
        double width,
        std::size_t event_count
    );
}


#endif //WIDTHSENSITIVITY_HPP