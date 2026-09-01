#ifndef FISHERINFORMATION_HPP
#define FISHERINFORMATION_HPP

#include <TH1.h>

class FisherInformation
{
public:
    FisherInformation() = delete;
    //richiede istogrammi anche non normalizzati
    FisherInformation(TH1D * nominal, TH1D * shifted, double parameter_delta, std::size_t accepted_events);
    //richiede che nominal sia normalizzato e derivative sia la sua derivata.
    FisherInformation(TH1D * nominal, TH1D * derivative, std::size_t accepted_events);

    FisherInformation(const FisherInformation&) = delete;
    FisherInformation& operator=(const FisherInformation&) = delete;

    FisherInformation(FisherInformation&&) = default;
    FisherInformation& operator=(FisherInformation&&) = default;

    std::size_t get_selected_events_count() const;
    double sigma() const;
    double sigma_uncertainty() const;
    double fisher() const;
    double fisher_uncertainty() const;
    ~FisherInformation() = default;
private:

    double m_sigma,
        m_sigma_uncertainty,
        m_fisher, 
        m_fisher_uncertainty;
};

#endif //FISHERINFORMATION_HPP