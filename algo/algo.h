#ifndef ACG_ALP_LDPC_ALGO_H
#define ACG_ALP_LDPC_ALGO_H

#include "../utils/channel.h"

class Decoder {
public:
    virtual pair<TCodeword, bool> decode(const TMatrix &H, const TFVector &channel_word, double snr) = 0;

    virtual string name() const = 0;
};

std::vector<double> CalculateCoef(const std::vector<double> &y, double snr) {
    std::vector<double> coef(y.size());
    for (int i = 0; i < (int) y.size(); ++i) {
        coef[i] = llr(y[i], snr);
    }

    return coef;
}

#endif //ACG_ALP_LDPC_ALGO_H
