#include <MassSensitivityAnalyzer.hpp>
#include <cmath>
#include <iostream>

MassSensitivityAnalyzer::MassSensitivityAnalyzer(TH1D* h0, TH1D* h1, double deltaMass)
{
    ratioHist = std::unique_ptr<TH1D>(
        dynamic_cast<TH1D*>(h1->Clone("h_ratio"))
    );
    

    ratioHist->SetTitle("Template ratio; p_{T}^{#mu} [GeV];N(m_{W}+#Delta m)/N(m_{W})");
    ratioHist->Divide(h0);
    
    double norm0 = h0->Integral("width");
    double norm1 = h1->Integral("width");

    double fisher = 0.0;

    for (int i = 1; i <= h0->GetNbinsX(); ++i)
    {
        double p0 = h0->GetBinContent(i) * h0->GetBinWidth(i) / norm0;
        double p1 = h1->GetBinContent(i) * h1->GetBinWidth(i) / norm1;
        double dpdm = (p1 - p0) / deltaMass;

        fisher += dpdm * dpdm / p0;
    }
    selectedEvents = h0->Integral();
    sigmaMass = 1.0 / std::sqrt(selectedEvents * fisher);
}
std::size_t MassSensitivityAnalyzer::get_selected_events_count() const
{
    return selectedEvents;
}

TH1D* MassSensitivityAnalyzer::getRatioHist() const
{
    return ratioHist.get();
}
std::unique_ptr<TH1D> MassSensitivityAnalyzer::releaseRatioHist()
{
    return std::move(ratioHist);
}

double MassSensitivityAnalyzer::sigma() const
{
    return sigmaMass;
}
