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
    
    double norm0 = h0->Integral();
    double norm1 = h1->Integral();

    ratioHist->Scale( norm1/norm0 );

    double fisher = 0.0;

    for (int i = 1; i <= h0->GetNbinsX(); ++i)
    {
        double p0 = h0->GetBinContent(i) / norm0;
        double p1 = h1->GetBinContent(i) / norm1;

        if (p0 <= 0.0)
            continue;

        double dpdm = (p1 - p0) / deltaMass;

        fisher += dpdm * dpdm / p0;
    }

    sigmaMass = 1.0 / std::sqrt(norm0 * fisher);
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
