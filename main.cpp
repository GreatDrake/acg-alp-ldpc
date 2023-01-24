//#define USE_GLPK

#include <utility>
#include <memory>
#include <chrono>

#include "utils/parse_data.h"
#include "algo/algo.h"
#include "algo/bp.h"
#include "algo/qp_admm.h"

#ifdef USE_GLPK
#include "algo/full_lp.h"
#include "algo/alp.h"
#include "algo/agc_alp.h"
#endif

using namespace std;

const int THREADS_NUM = 26;
const int LOG_FREQ = 2000;
const int TESTS_NUM = 1000;
const vector<double> SNRS = {-5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5};
const vector <shared_ptr<Decoder>> decoders{
        make_shared<BeliefPropagationDecoder>(100),
        make_shared<QPADMMDecoder>(0.6, 1.0, 1000, 1e-5),
#ifdef USE_GLPK
        make_shared<FullLPDecoder>(),
        make_shared<ALPDecoder>(),
        make_shared<AGCALPDecoder>(2000),
#endif
};

struct ThreadArgs {
    ThreadArgs(shared_ptr <Decoder> decoder, const vector <TCodeword> &codewords, const TMatrix &H, double snr) :
            decoder(std::move(decoder)), codewords(codewords), H(H), snr(snr) {}

    shared_ptr <Decoder> decoder;
    const vector <TCodeword> codewords;
    const TMatrix &H;
    int last_uuid = 0;
    double snr;
    pthread_mutex_t uuid_id_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;
};

struct HammingDistanceTracker {
    HammingDistanceTracker(int sum_hamming = 0, int sum_hamming_ok = 0, int sum_hamming_wrong = 0) :
            sum_hamming(sum_hamming), sum_hamming_ok(sum_hamming_ok), sum_hamming_wrong(sum_hamming_wrong) {}

    int sum_hamming;
    int sum_hamming_ok;
    int sum_hamming_wrong;

    void new_experiment(const TMatrix &H, const TCodeword &c, const TFVector &y, bool correct) {
        int hamming = 0;
        for (int i = 0; i < (int) H[0].size(); ++i) {
            if (!c[i] && y[i] <= 0)
                hamming++;
            if (c[i] && y[i] > 0)
                hamming++;
        }
        sum_hamming += hamming;
        if (correct)
            sum_hamming_ok += hamming;
        else
            sum_hamming_wrong += hamming;
    }
};

struct ExperimentResult {
    ExperimentResult(HammingDistanceTracker tr, int correct = 0, int pseudo = 0, int total = 0, double time_sec = 0) :
            tr(tr), correct(correct), pseudo(pseudo), total(total), time_sec(time_sec) {}

    HammingDistanceTracker tr;
    int correct;
    int pseudo;
    int total;
    double time_sec;

    double FER() { return (double)(total - correct) / total; }

    double avg_time() { return time_sec / total; }

    double mean_hamming() { return (double)tr.sum_hamming / total; }

    double mean_hamming_ok() { return (double)tr.sum_hamming_ok / max(1, correct); }

    double mean_hamming_wrong() { return (double)tr.sum_hamming_wrong / max(1, total - correct); }
};

void merge_exp_results(ExperimentResult &a, const ExperimentResult &b) {
    a.tr.sum_hamming += b.tr.sum_hamming;
    a.tr.sum_hamming_ok += b.tr.sum_hamming_ok;
    a.tr.sum_hamming_wrong += b.tr.sum_hamming_wrong;
    a.correct += b.correct;
    a.pseudo += b.pseudo;
    a.total += b.total;
    a.time_sec += b.time_sec;
}

void *exp(void *arg) {
    ThreadArgs &args = *(ThreadArgs *) arg;
    int correct = 0, total = 0, pseudo = 0;
    double time = 0;
    HammingDistanceTracker tr;
    while (true) {
        assert(!pthread_mutex_lock(&args.uuid_id_mutex));
        TCodeword codeword;
        if (args.last_uuid < (int) args.codewords.size()) {
            codeword = args.codewords[args.last_uuid];
            args.last_uuid++;
        }
        int uuid = args.last_uuid;
        assert(!pthread_mutex_unlock(&args.uuid_id_mutex));
        if (codeword.empty())
            break;

        mt19937 rnd(args.last_uuid);

        TFVector transmitted = transmit(args.snr, codeword, rnd);
        auto start_time = chrono::steady_clock::now();
        pair<TCodeword, bool> p = args.decoder->decode(args.H, transmitted, args.snr);
        auto duration = chrono::steady_clock::now() - start_time;
        time += (double)duration.count() / 1e9;
        if (uuid % LOG_FREQ == 0) {
            assert(!pthread_mutex_lock(&args.log_mutex));
            cout << uuid << endl;
            assert(!pthread_mutex_unlock(&args.log_mutex));
        }
        bool is_correct = false;
        if (p.second == 1) {
            if (IsCodeword(args.H, p.first)) {
                if (p.first == codeword) {
                    correct++;
                    is_correct = true;
                } else
                    pseudo++;
            }
        }
        total++;
        tr.new_experiment(args.H, codeword, transmitted, is_correct);
    }
    return (void *) new ExperimentResult(tr, correct, pseudo, total, time);
}

int main() {
    std::ios::sync_with_stdio(0);
    cout.precision(5);
    cout << fixed;

//    glp_term_out(GLP_MSG_OFF);

    TMatrix H = read_pcm("data/H05.txt");
//    vector<TCodeword> codewords = read_codewords("data/codewords.txt");
//    codewords.resize(TESTS_NUM);
    TMatrix G = read_pcm("data/G05.txt");
    mt19937 rnd(239);
    vector <TCodeword> codewords = gen_random_codewords(G, TESTS_NUM, rnd);

    cerr << "n=" << H[0].size() << " k=" << H.size() << "\n";

    for (auto decoder: decoders) {
        cout << "Algo: " << decoder->name() << endl;
        for (double snr: SNRS) {
            ThreadArgs args(decoder, codewords, H, snr);
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

            int total = codewords.size();
            cout << "\tSNR: " << snr << ", FER: " << res.FER() << ", (time=" << res.avg_time() << "s)" << endl;
            cerr << "\tSNR: " << snr << ", FER: " << res.FER() << ", (time=" << res.avg_time() << "s)" << endl;
            cerr << "\t\tAverage hamming distance: " << res.mean_hamming() << std::endl;
            cerr << "\t\tAverage hamming distance, correctly decoded: " << res.mean_hamming_ok() << endl;
            cerr << "\t\tAverage hamming distance, incorrectly decoded: " << res.mean_hamming_wrong() << endl;
        }
        cerr << string(30, '_') << endl;
    }

    return 0;
}
