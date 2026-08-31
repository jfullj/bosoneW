#ifndef HISTUTILS_HPP
#define HISTUTILS_HPP

#include <TH1.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TCanvas.h>


#include <memory>
#include <string>
#include <concepts>

class TemplateComparison
{
public:
    struct CanvasPropeties{
        std::size_t width;
        std::size_t height;
        double nominal_mass;
        double shifted_mass;
        double pad_height_ratio;
        double gap;
        int nominal_hist_color;
        int shifted_hist_color;
        std::string upper_x_axis_content;
        std::string lower_x_axis_content;
        std::string upper_y_axis_content;
        std::string lower_y_axis_content;
        double title_size;
        double label_size;
        double title_offset_x;
        double title_offset_y;
    };

    TemplateComparison(TH1* nominal, TH1* shifted, double sigmaW, const CanvasPropeties& props);

    TemplateComparison(const TemplateComparison&) = delete;
    TemplateComparison& operator=(const TemplateComparison&) = delete;

    TemplateComparison(TemplateComparison&&) = default;
    TemplateComparison& operator=(TemplateComparison&&) = delete;

    void save_as(const std::string& path);

    ~TemplateComparison() = default;
private:
    std::unique_ptr<TH1> Nominal, Shifted, Ratio, Error;
    CanvasPropeties CanvasProps;
    double SigmaW;
};


struct SavePlotParams
{
    std::size_t width;
    std::size_t height; 
    std::string name;
    std::string x_axis_content;
    std::string y_axis_content;
    int color;
    double title_size;
    double label_size;
    double title_offset_x;
    double title_offset_y;
    std::string draw_settings;
};

namespace save_plot_impl
{
    template<typename Drawable>
    auto init_copy(Drawable* d, const SavePlotParams& params)
    {
        auto copy{ std::unique_ptr<Drawable>(dynamic_cast<Drawable*>(d->Clone())) };
        copy->SetTitle(params.name.c_str());

        copy->GetXaxis()->SetTitle(params.x_axis_content.c_str());
        copy->GetYaxis()->SetTitle(params.y_axis_content.c_str());
    
        copy->GetXaxis()->CenterTitle();
        copy->GetYaxis()->CenterTitle();


        copy->GetXaxis()->SetTitleSize(params.title_size);
        copy->GetXaxis()->SetLabelSize(params.label_size);
        copy->GetXaxis()->SetTitleOffset(params.title_offset_x);
    
        copy->GetYaxis()->SetTitleSize(params.title_size);
        copy->GetYaxis()->SetLabelSize(params.label_size);
        copy->GetYaxis()->SetTitleOffset(params.title_offset_y);

        copy->SetLineColor(params.color);
        copy->SetLineWidth(2);
        copy->SetFillStyle(0);
        
        return copy;
    }
}

template<std::derived_from<TObject> Drawable>
void save_plot(Drawable* d, const std::string& path, const SavePlotParams& params)
{
    auto copy{ save_plot_impl::init_copy(d, params) };
    TCanvas canvas{
        "", "", 
        static_cast<Int_t>(params.width),
        static_cast<Int_t>(params.height)
    };

    canvas.cd();
    
    copy->SetStats(false);
    copy->Draw(params.draw_settings.c_str());

    canvas.SaveAs(path.c_str());
}

template<std::derived_from<TObject> Drawable, std::invocable<Drawable*> F>
void save_plot(Drawable* d, const std::string& path, const SavePlotParams& params, F&& f)
{
    auto copy{ save_plot_impl::init_copy(d, params) };
    TCanvas canvas{
        "", "", 
        static_cast<Int_t>(params.width),
        static_cast<Int_t>(params.height)
    };

    canvas.cd();
    
    copy->Draw(params.draw_settings.c_str());

    f(copy.get());

    canvas.SaveAs(path.c_str());
}

#endif //HISTUTILS_HPP