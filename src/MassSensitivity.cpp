#include <MassSensitivity.hpp>
#include <fstream>

#include <Spectrum.hpp>
#include <FisherInformation.hpp>
#include <HistUtils.hpp>

#include <Spectrum.hpp>
#include <FisherInformation.hpp>
#include <sstream>

void Mass::junk(double wmass0, double wmass1, double wwidth, std::size_t event_count)
{
    const auto wdelta{ wmass1 - wmass0 };

    auto pT_gen{ std::make_unique<PT_Generator>() };
    save_plot(pT_gen->get_hist(), pT_dist_dir, SavePlotParams{
        .width = CANVAS_WIDTH,
        .height = CANVAS_HEIGHT,
        .name = "pTW distribution",
        .x_axis_content = "p^{T}_{W} [GeV]",
        .y_axis_content = "pmf",
        .color = kBlack,
        .title_size = 0.055,
        .label_size = 0.045,
        .title_offset_x = 0.80,
        .title_offset_y = 0.80,
        .draw_settings = ""
    });

    auto W_mass_gen0{ std::make_unique<BreitWignerGenerator>(wmass0, wwidth) };
    auto w_gen0{ std::make_unique<W_Generator>(
        W_mass_gen0.get(), 
        pT_gen.get()
    )};
    auto W_mass_gen1{ std::make_unique<BreitWignerGenerator>(wmass1, wwidth) };
    auto w_gen1{ std::make_unique<W_Generator>(
        W_mass_gen1.get(),
        pT_gen.get()
    )};

    WDecaySampler sampler0{ w_gen0.get() };
    WDecaySampler sampler1{ w_gen1.get() };

    auto pdf0{ Spectrum{sampler0, event_count}.release_hist() };
    auto pdf1{ Spectrum{sampler1, event_count}.release_hist() };

    FisherInformation fi{
        pdf0.get(),
        pdf1.get(),
        wdelta
    };

    double sigma{ fi.sigma() };

    TemplateComparison tc(pdf0.get(), pdf1.get(), sigma, {
        .width = CANVAS_WIDTH,
        .height = CANVAS_HEIGHT,
        .nominal_mass = wmass0,
        .shifted_mass = wmass1,
        .pad_height_ratio = 2,
        .gap = 0.02,
        .nominal_hist_color = kBlue + 1,
        .shifted_hist_color = kRed + 1,
        .upper_x_axis_content = "",
        .lower_x_axis_content = "p_{T}^{#mu} [GeV]",
        .upper_y_axis_content = "events / N",
        .lower_y_axis_content = "#frac{f_{m+#Deltam}}{f_{m}}",
        .title_size = 0.055,
        .label_size = 0.045,
        .title_offset_x = 1.20,
        .title_offset_y = 1.20
    });

    tc.save_as(output_hist_file);

    {
        auto pT_gen{ std::make_unique<PT_Delta_Generator>(8) };
        auto W_gen{ std::make_unique<BreitWignerGenerator>(wmass0, wwidth) };
        auto w_gen{ std::make_unique<W_Generator>(
            W_gen.get(),
            pT_gen.get()
        )}; 

        WDecaySampler sampler{ w_gen.get() };

        Spectrum spectrum{
            sampler,
            event_count,
            Event::Transform::standard,
            Binning::Parameters{
                .bin_count = 100,
                .min = 0,
                .max = 100
            }, Acceptance::all
        };

        auto pdf{ spectrum.release_hist() };
        pdf->Scale(1 / pdf->Integral("width"));

        save_plot(spectrum.get_hist(), DATA_DIR "/results/pTW const.png", SavePlotParams{
            .width = CANVAS_WIDTH,
            .height = CANVAS_HEIGHT,
            .name = "",
            .x_axis_content = "p_{T}^{#mu} [GeV]",
            .y_axis_content = "pdf",
            .color = kBlue + 1,
            .title_size = 0.055,
            .label_size = 0.045,
            .title_offset_x = 0.80,
            .title_offset_y = 0.60,
            .draw_settings = "HIST SAME"
        });
    }

    std::ostringstream content;
    content << "sigma_mass = " << sigma << " GeV\n\n"
            << "generated events = " << event_count << "\n"
            << "selected events = " << fi.get_selected_events_count() << "\n"
            << "acceptance mass 0 = " << (static_cast<double>(pdf0->Integral()) / event_count) << "\n"
            << "acceptance mass 1 = " << (static_cast<double>(pdf1->Integral()) / event_count) << "\n";

    std::ofstream output{ output_file };
    if (!output)
        throw std::runtime_error{ "Impossibile creare sigma_mass.txt" };

    output << content.str();
    output.close();

    std::cout << content.str();
}