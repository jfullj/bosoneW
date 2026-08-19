#ifndef HISTUTILS_HPP
#define HISTUTILS_HPP

#include <TH1.h>
#include <TStyle.h>
#include <TPaveText.h>

#include <memory>
#include <string>

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

    TemplateComparison(TH1* nominal, TH1* shifted,const CanvasPropeties& props);

    TemplateComparison(const TemplateComparison&) = delete;
    TemplateComparison& operator=(const TemplateComparison&) = delete;

    TemplateComparison(TemplateComparison&&) = default;
    TemplateComparison& operator=(TemplateComparison&&) = delete;

    void save_as(const std::string& path);

    ~TemplateComparison() = default;
private:
    std::unique_ptr<TH1> Nominal, Shifted, Ratio, Error;
    
    CanvasPropeties CanvasProps;
};

#endif //HISTUTILS_HPP