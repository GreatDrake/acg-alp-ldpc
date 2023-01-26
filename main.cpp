#define USE_GLPK

#include <utility>
#include <memory>
#include <iomanip>

#include "experiment.h"
#include "utils/parse_data.h"
#include "utils/codeword.h"
#include "algo/algo.h"
#include "algo/bp.h"
#include "algo/qp_admm.h"

#ifdef USE_GLPK
#include "algo/full_lp.h"
#include "algo/alp.h"
#include "algo/agc_alp.h"
#endif

using namespace std;

const int THREADS_NUM = 8;
const int LOG_FREQ = 1000000;
const int TESTS_NUM = 1000;

const vector<double> SNRS = {-5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5};
const vector <shared_ptr<Decoder>> decoders{
        make_shared<BeliefPropagationDecoder>(100),
        make_shared<QPADMMDecoder>(1.95, 0.5, 5000, 1e-5),
#ifdef USE_GLPK
        //make_shared<FullLPDecoder>(),
        make_shared<ALPDecoder>(),
        make_shared<AGCALPDecoder>(1000),
#endif
};

int main() {
    std::ios::sync_with_stdio(0);
    cout.precision(5);
    cout << fixed;

    ofstream fdata("report.csv");
    fdata << "Method,SNR,Sigma,FER,Time,AvgHamming,AvgHammingCorrect,AvgHammingWrong" << endl;
    fdata << fixed << setprecision(12);

    for (auto snr: SNRS) {
        cerr << "snr=" << snr << ": var=" << llr_variance(snr) << endl;
    }

    TMatrix H = read_pcm("data/optimalH.txt");
    //    TMatrix G = read_pcm("data/G05.txt");
    TMatrix G = GetOrtogonal(H).first;
//    vector<TCodeword> codewords = read_codewords("data/codewords.txt");
//    codewords.resize(TESTS_NUM);
    mt19937 rnd(239);
    vector <TCodeword> codewords = gen_random_codewords(G, TESTS_NUM, rnd);

    cerr << "n=" << H[0].size() << " k=" << H.size() << "\n";

    for (auto decoder: decoders) {
        cout << "Algo: " << decoder->name() << endl;
        for (double snr: SNRS) {
            ExperimentResult res = multithread_experiment(decoder, codewords, H, snr, THREADS_NUM, LOG_FREQ);

            cout << "\tSNR: " << snr << ", FER: " << res.FER() << ", (time=" << res.avg_time() << "s)" << endl;
            cerr << "\tSNR: " << snr << ", FER: " << res.FER() << ", (time=" << res.avg_time() << "s)" << endl;
            cerr << "\t\tAverage hamming distance: " << res.mean_hamming() << std::endl;
            cerr << "\t\tAverage hamming distance, correctly decoded: " << res.mean_hamming_ok() << endl;
            cerr << "\t\tAverage hamming distance, incorrectly decoded: " << res.mean_hamming_wrong() << endl;

            fdata << decoder->name() << ",";
            fdata << snr << ",";
            fdata << sqrt(llr_variance(snr)) << ",";
            fdata << res.FER() << ",";
            fdata << res.avg_time() << ",";
            fdata << res.mean_hamming() << ",";
            fdata << res.mean_hamming_ok() << ",";
            fdata << res.mean_hamming_wrong() << endl;
        }
        cerr << string(30, '_') << endl;
    }

    return 0;
}
