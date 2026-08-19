#include <HistUtils.hpp>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>


TemplateComparison::TemplateComparison(TH1* nominal, TH1* shifted,const CanvasPropeties& props)
: CanvasProps{ props }
{
    if (!nominal || !shifted)
        throw std::runtime_error("Null histogram");

    Nominal = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(nominal->Clone()));
    Shifted = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(shifted->Clone()));

    Nominal->SetDirectory(nullptr);
    Shifted->SetDirectory(nullptr);

    Nominal->SetTitle("");
    Shifted->SetTitle("");

    const double norm0 = Nominal->Integral("width");
    const double norm1 = Shifted->Integral("width");

    Nominal->Scale(1.0 / norm0);
    Shifted->Scale(1.0 / norm1);
    
    Nominal->GetXaxis()->SetTitle(CanvasProps.upper_x_axis_content.c_str());
    Nominal->GetYaxis()->SetTitle(CanvasProps.upper_y_axis_content.c_str());

    Nominal->GetXaxis()->CenterTitle();
    Nominal->GetYaxis()->CenterTitle();

    Nominal->GetXaxis()->SetTitleSize(CanvasProps.title_size);
    Nominal->GetXaxis()->SetLabelSize(CanvasProps.label_size);
    Nominal->GetXaxis()->SetTitleOffset(CanvasProps.title_offset_x);
    
    Nominal->GetYaxis()->SetTitleSize(CanvasProps.title_size);
    Nominal->GetYaxis()->SetLabelSize(CanvasProps.label_size);
    Nominal->GetYaxis()->SetTitleOffset(CanvasProps.title_offset_y);

    Ratio = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(Shifted->Clone("h_ratio_plot")));
    Ratio->SetDirectory(nullptr);
    Ratio->SetTitle("");
    
    Ratio->Divide(Nominal.get());

    double ratio_min = 1.0;
    double ratio_max = 1.0;

    for (int i = 1; i <= Ratio->GetNbinsX(); ++i)
    {
        const double y = Ratio->GetBinContent(i);
        const double e = Ratio->GetBinError(i);

        ratio_min = std::min(ratio_min, y - e);
        ratio_max = std::max(ratio_max, y + e);
    }

    const double range = ratio_max - ratio_min;

    Ratio->SetMinimum(ratio_min - 0.15 * range);
    Ratio->SetMaximum(ratio_max + 0.15 * range);

    Ratio->GetXaxis()->SetTitle(CanvasProps.lower_x_axis_content.c_str());
    Ratio->GetYaxis()->SetTitle(CanvasProps.lower_y_axis_content.c_str());
    
    Ratio->GetXaxis()->CenterTitle();
    Ratio->GetYaxis()->CenterTitle();


    Ratio->GetXaxis()->SetTitleSize(CanvasProps.title_size * CanvasProps.pad_height_ratio);
    Ratio->GetXaxis()->SetLabelSize(CanvasProps.label_size * CanvasProps.pad_height_ratio);
    Ratio->GetXaxis()->SetTitleOffset(CanvasProps.title_offset_x);
    
    Ratio->GetYaxis()->SetTitleSize(CanvasProps.title_size * CanvasProps.pad_height_ratio);
    Ratio->GetYaxis()->SetLabelSize(CanvasProps.label_size * CanvasProps.pad_height_ratio);
    Ratio->GetYaxis()->SetTitleOffset(CanvasProps.title_offset_y / CanvasProps.pad_height_ratio);

    Ratio->GetYaxis()->SetNdivisions(504, false);

    Error = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(Nominal->Clone()));
    Error->SetTitle("");
    Error->SetDirectory(nullptr);

    Error->SetLineColor(kBlack);
    Error->SetMarkerColor(kBlack);

    Error->SetMarkerStyle(20);
    Error->SetMarkerSize(0.65);

    Nominal->SetStats(false);
    Shifted->SetStats(false);
    Ratio->SetStats(false);
    Error->SetStats(false);
}
void TemplateComparison::save_as(const std::string& path)
{
    TCanvas canvas{
        "c_template_comparison",
        "Template comparison",
        static_cast<Int_t>(CanvasProps.width),
        static_cast<Int_t>(CanvasProps.height)
    };
    canvas.SetFillColor(kWhite);

    double lower_height{ 1 / (1 + CanvasProps.pad_height_ratio) };

    TPad upper{ "upper", "", 0.0, lower_height + CanvasProps.gap / 2., 1.0, 1.0 };
    TPad lower{ "lower", "", 0.0, 0.0, 1.0, lower_height - CanvasProps.gap / 2.};

    upper.SetLeftMargin(0.14);
    upper.SetRightMargin(0.04);
    upper.SetTopMargin(0.07);
    upper.SetBottomMargin(0.08);

    lower.SetLeftMargin(0.14);
    lower.SetRightMargin(0.04);
    lower.SetTopMargin(0.02);
    lower.SetBottomMargin(0.32);

    upper.SetTicks(1, 1);
    lower.SetTicks(1, 1);

    canvas.cd();
    upper.Draw();
    lower.Draw();

    upper.cd();

    Nominal->SetLineColor(CanvasProps.nominal_hist_color);
    Nominal->SetLineWidth(2);
    Nominal->SetFillStyle(0);

    Shifted->SetLineColor(CanvasProps.shifted_hist_color);
    Shifted->SetLineWidth(2);
    Shifted->SetFillStyle(0);

    {
        const double ymax = std::max(
            Nominal->GetMaximum(),
            Shifted->GetMaximum()
        );

        Nominal->SetMaximum(1.25 * ymax);
        Nominal->SetMinimum(0.0);
    }
    Nominal->Draw("HIST");
    Shifted->Draw("HIST SAME");

    Error->Draw("E1 X0 SAME");

    TLegend legend{
        0.58, 0.72,
        0.91, 0.89
    };

    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.037);

    legend.AddEntry(
        Nominal.get(),
        Form("m_{W} = %.3f GeV", CanvasProps.nominal_mass),
        "l"
    );

    legend.AddEntry(
        Shifted.get(),
        Form("m_{W} = %.3f GeV", CanvasProps.nominal_mass),
        "l"
    );

    legend.AddEntry(
        Error.get(),
        "MC statistical uncertainty",
        "lep"
    );

    legend.Draw();

    lower.cd();

    Ratio->SetMarkerStyle(20);
    Ratio->SetMarkerSize(0.65);

    Ratio->Draw("E1 X0");

    TLine unity{
        Nominal->GetXaxis()->GetXmin(),
        1.0,
        Nominal->GetXaxis()->GetXmax(),
        1.0
    };

    unity.SetLineColor(kBlack);
    unity.SetLineStyle(2);
    unity.SetLineWidth(1);

    unity.Draw("SAME");

    canvas.cd();
    canvas.SaveAs(path.c_str());
}
