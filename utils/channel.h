#ifndef ACG_ALP_LDPC_CHANNEL_H
#define ACG_ALP_LDPC_CHANNEL_H

#include "codeword.h"

using namespace std;

#define double long double
typedef vector<double> TFVector;

const double EPS = 1e-8;

double llr_variance(double snr) { return pow(10, -(snr / 10)) / 2; }

double llr(double v, double snr) {
    return 2 * v / llr_variance(snr);
}

template<typename Gen>
TFVector transmit(double snr, const TCodeword &c, Gen &rnd) {
    double llrSigma = sqrt(llr_variance(snr));
    TFVector res((int) c.size());
    normal_distribution<double> dst(0, llrSigma);
    for (int i = 0; i < (int) c.size(); i++)
        res[i] = (c[i] ? -1.0 : 1.0) + dst(rnd);
    return res;
}

template<typename Gen>
TCodeword gen_random_codeword(const std::vector <TCodeword> &G, Gen &rnd) {
    assert(!G.empty());
    TCodeword res(G[0].size(), false);
    for (int i = 0; i < (int) G.size(); i++)
        if (rnd() % 2 == 0)
            res = (res ^ G[i]);
    return res;
}

template<typename Gen>
vector <TCodeword> gen_random_codewords(const TMatrix &G, int n, Gen &rnd) {
    vector <TCodeword> tests;
    for (int i = 0; i < n; i++)
        tests.push_back(gen_random_codeword(G, rnd));
    return tests;
}

#endif //ACG_ALP_LDPC_CHANNEL_H
