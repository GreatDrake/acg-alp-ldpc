#include <utility>
#include <memory>

#include "experiment.h"

#include "utils/parse_data.h"
#include "algo/algo.h"
#include "algo/qp_admm.h"

using namespace std;

const int THREADS_NUM = 8;
const int LOG_FREQ = 1000000;
const int TESTS_NUM = 1000;

double estimate_qpadmm(const TMatrix& H, const vector<TCodeword>& codewords, double snr, double alpha, double mu) {
    auto decoder = make_shared<QPADMMDecoder>(alpha, mu, 1000, 1e-5);
    ThreadArgs args(decoder, codewords, H, snr, LOG_FREQ);
    vector <pthread_t> threads(THREADS_NUM);
    for (int q = 0; q < THREADS_NUM; q++)
        assert(!pthread_create(&threads[q], nullptr, exp, (void *) &args));
    ExperimentResult res = ExperimentResult(HammingDistanceTracker());
    for (int q = 0; q < THREADS_NUM; q++) {
        void *data;
        assert(!pthread_join(threads[q], &data));
        merge_exp_results(res, *(ExperimentResult *) data);
        delete (ExperimentResult *) data;
    }
    return res.FER();
}

double linear_function(double L, double R, int cnt, int i) {
    return L + ((R - L) / (cnt - 1)) * i;
}

int main() {
    std::ios::sync_with_stdio(0);
    cout.precision(5);
    cout << fixed;

    //TMatrix H = read_pcm("data/H05.txt");
    //TMatrix G = read_pcm("data/G05.txt");
    TMatrix H = read_pcm("data/optimalH.txt");
    TMatrix G = GetOrtogonal(H).first;

    mt19937 rnd(239);
    vector <TCodeword> codewords = gen_random_codewords(G, TESTS_NUM, rnd);

    cerr << "n=" << H[0].size() << " k=" << H.size() << endl;

    double alpha_l = 0;
    double alpha_r = 3.0;
    int alpha_cnt = 61;
    double mu_l = 0;
    double mu_r = 3.0;
    int mu_cnt = 61;

    double snr = -3.0;

    double best_fer = 2.0;
    double best_alpha = -1;
    double best_mu = -1;

    for (int alpha_i = 0; alpha_i < alpha_cnt; ++alpha_i) {
        for (int mu_i = 0; mu_i < mu_cnt; ++mu_i) {
            double alpha = linear_function(alpha_l, alpha_r, alpha_cnt, alpha_i);
            double mu = linear_function(mu_l, mu_r, mu_cnt, mu_i);
            double fer = estimate_qpadmm(H, codewords, snr, alpha, mu);
            cerr << "alpha=" << alpha << ", mu=" << mu << ": fer=" << fer << endl;
            if (fer < best_fer) {
                best_fer = fer;
                best_alpha = alpha;
                best_mu = mu;
                cout << "new best fer found: " << fer << "| alpha=" << alpha << ", mu=" << mu << endl;
            }
        }
    }

    cout << "Best parameters:" << endl;
    cout << "alpha=" << best_alpha << endl;
    cout << "mu=" << best_mu << endl;
    cout << "fer=" << best_fer << endl;

    return 0;
}
