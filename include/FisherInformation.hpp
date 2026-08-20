#ifndef FISHERINFORMATION_HPP
#define FISHERINFORMATION_HPP

#include <TH1.h>

class FisherInformation
{
public:
    FisherInformation() = delete;
    FisherInformation(TH1D* h0, TH1D* h1, double deltaMass);

    FisherInformation(const FisherInformation&) = delete;
    FisherInformation& operator=(const FisherInformation&) = delete;

    FisherInformation(FisherInformation&&) = default;
    FisherInformation& operator=(FisherInformation&&) = default;

    TH1D* getRatioHist() const;
    std::unique_ptr<TH1D> releaseRatioHist();

    std::size_t get_selected_events_count() const;
    double sigma() const;
    ~FisherInformation() = default;
private:

    std::unique_ptr<TH1D> ratioHist;
    double sigmaMass;
    std::size_t selectedEvents;
};

#endif //FISHERINFORMATION_HPP