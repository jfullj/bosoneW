#ifndef MASSSENSITIVITYANALYZER_HPP
#define MASSSENSITIVITYANALYZER_HPP

#include <TH1.h>

class MassSensitivityAnalyzer
{
public:
    MassSensitivityAnalyzer() = delete;
    MassSensitivityAnalyzer(TH1D* h0, TH1D* h1, double deltaMass);

    MassSensitivityAnalyzer(const MassSensitivityAnalyzer&) = delete;
    MassSensitivityAnalyzer& operator=(const MassSensitivityAnalyzer&) = delete;

    MassSensitivityAnalyzer(MassSensitivityAnalyzer&&) = default;
    MassSensitivityAnalyzer& operator=(MassSensitivityAnalyzer&&) = default;

    TH1D* getRatioHist() const;
    std::unique_ptr<TH1D> releaseRatioHist();

    std::size_t get_selected_events_count() const;
    double sigma() const;
    ~MassSensitivityAnalyzer() = default;
private:

    std::unique_ptr<TH1D> ratioHist;
    double sigmaMass;
    std::size_t selectedEvents;
};

#endif //MASSSENSITIVITYANALYZER_HPP