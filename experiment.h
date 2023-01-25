#ifndef ACG_ALP_LDPC_EXPERIMENT_H
#define ACG_ALP_LDPC_EXPERIMENT_H

#include <chrono>
#include "utils/codeword.h"
#include "algo/algo.h"

using namespace std;

struct ThreadArgs {
    ThreadArgs(shared_ptr <Decoder> decoder, const vector <TCodeword> &codewords,
               const TMatrix &H, double snr, int log_freq) :
            decoder(std::move(decoder)), codewords(codewords), H(H), snr(snr), log_freq(log_freq) {}

    shared_ptr <Decoder> decoder;
    const vector <TCodeword> codewords;
    const TMatrix &H;
    int last_uuid = 0;
    double snr;
    int log_freq;
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

    double FER() { return (double) (total - correct) / total; }

    double avg_time() { return time_sec / total; }

    double mean_hamming() { return (double) tr.sum_hamming / total; }

    double mean_hamming_ok() { return (double) tr.sum_hamming_ok / max(1, correct); }

    double mean_hamming_wrong() { return (double) tr.sum_hamming_wrong / max(1, total - correct); }
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
        time += (double) duration.count() / 1e9;
        if (uuid % args.log_freq == 0) {
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

ExperimentResult multithread_experiment(shared_ptr <Decoder> decoder, const vector <TCodeword> &codewords,
                                        const TMatrix &H, double snr, int threads_num, int log_freq = 1e9) {
    ThreadArgs args(decoder, codewords, H, snr, log_freq);
    vector <pthread_t> threads(threads_num);
    for (int q = 0; q < threads_num; q++)
        assert(!pthread_create(&threads[q], nullptr, exp, (void *) &args));
    ExperimentResult res = ExperimentResult(HammingDistanceTracker());
    for (int q = 0; q < threads_num; q++) {
        void *data;
        assert(!pthread_join(threads[q], &data));
        merge_exp_results(res, *(ExperimentResult *) data);
        delete (ExperimentResult *) data;
    }
    return res;
}

#endif //ACG_ALP_LDPC_EXPERIMENT_H
