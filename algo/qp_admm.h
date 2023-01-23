#ifndef ACG_ALP_LDPC_QP_ADMM_H
#define ACG_ALP_LDPC_QP_ADMM_H

#include "algo.h"

struct ADMMProblem {
    std::vector<double> q;
    std::vector<std::vector<std::pair<int, double>>> A;
    std::vector<double> b;
    std::vector<double> e;
};

ADMMProblem ConstructADMMProblem(const std::vector<TCodeword> &H, const std::vector<double> &y, double snr) {
    int n_aux = 0;
    for (int i = 0; i < (int) H.size(); ++i) {
        int sz = 0;
        for (int t = 0; t < (int) H[i].size(); ++t) {
            sz += H[i][t];
        }
        n_aux += std::max(sz - 3, 0);
    }
    int n_var = y.size() + n_aux;

    std::vector<double> q(n_var, 0.0);
    std::vector<std::vector<std::pair<int, double>>> A(n_var);

    auto coef = CalculateCoef(y, snr);
    for (int i = 0; i < (int) y.size(); ++i) {
        q[i] = coef[i];
    }

    std::vector<double> b;

    auto add_three = [&A, &b](int i, int j, int h) {
        b.emplace_back(0.0);
        b.emplace_back(0.0);
        b.emplace_back(0.0);
        b.emplace_back(2.0);
        {
            A[i].emplace_back((int) b.size() - 4, 1.0);
            A[i].emplace_back((int) b.size() - 3, -1.0);
            A[i].emplace_back((int) b.size() - 2, -1.0);
            A[i].emplace_back((int) b.size() - 1, 1.0);
        }
        {
            A[j].emplace_back((int) b.size() - 4, -1.0);
            A[j].emplace_back((int) b.size() - 3, 1.0);
            A[j].emplace_back((int) b.size() - 2, -1.0);
            A[j].emplace_back((int) b.size() - 1, 1.0);
        }
        {
            A[h].emplace_back((int) b.size() - 4, -1.0);
            A[h].emplace_back((int) b.size() - 3, -1.0);
            A[h].emplace_back((int) b.size() - 2, 1.0);
            A[h].emplace_back((int) b.size() - 1, 1.0);
        }
    };

    int pos = y.size();
    for (int i = 0; i < (int) H.size(); ++i) {
        std::vector<int> idx;
        for (int j = 0; j < (int) y.size(); ++j) {
            if (H[i][j]) {
                idx.emplace_back(j);
            }
        }
        if (idx.empty()) {
            continue;
        }
        if (idx.size() == 1) {
            A[idx[0]].emplace_back(b.size(), 1.0);
            b.emplace_back(0.0);
            continue;
        }
        if (idx.size() == 2) {
            b.emplace_back(0.0);
            b.emplace_back(0.0);
            A[idx[0]].emplace_back((int) b.size() - 2, 1.0);
            A[idx[0]].emplace_back((int) b.size() - 1, -1.0);
            A[idx[1]].emplace_back((int) b.size() - 2, -1.0);
            A[idx[1]].emplace_back((int) b.size() - 1, 1.0);
            continue;
        }
        int idx_last = idx[0];
        for (int j = 1; j < (int) idx.size() - 2; ++j) {
            int idx_mid = idx[j];
            int idx_aux = pos++;
            add_three(idx_last, idx_mid, idx_aux);
            idx_last = idx_aux;
        }
        add_three(idx_last, idx[(int) idx.size() - 2], idx.back());
    }

    std::vector<double> e(n_var, 0.0);
    for (int i = 0; i < (int) A.size(); ++i) {
        for (auto f: A[i]) {
            e[i] += f.second * f.second;
        }
    }

    return ADMMProblem{q, A, b, e};
}

pair<TCodeword, bool> DecodeQPADMM(const std::vector<TCodeword> &H, const std::vector<double> &y,
                                   double snr, double alpha, double mu, int max_iter, double eps_stop) {
    auto prob = ConstructADMMProblem(H, y, snr);

    double e_min = 1e9;
    for (double e: prob.e) {
        e_min = std::min(e_min, e);
    }
    assert(e_min * mu > alpha);

    std::vector<double> v(prob.q.size(), 0.0);
    for (int i = 0; i < (int) prob.q.size(); ++i) {
        v[i] = prob.q[i] > 0.0 ? 1.0 : 0.0;
    }
    std::vector<double> z(prob.b.size(), 0.0);
    std::vector<double> yl(prob.b.size(), 0.0);

    std::vector<double> r(prob.b.size());
    for (int iter = 0; iter < max_iter; ++iter) {
        // update v
        for (int i = 0; i < (int) v.size(); ++i) {
            double A = (mu * prob.e[i] - alpha) / 2;
            double B = prob.q[i] + (alpha / 2);
            for (auto p: prob.A[i]) {
                auto j = p.first;
                auto cf = p.second;
                B += cf * (yl[j] + mu * (z[j] - prob.b[j]));
            }
            v[i] = -B / (2 * A);
            v[i] = std::max(v[i], 0.0);
            v[i] = std::min(v[i], 1.0);
        }

        for (int i = 0; i < (int) r.size(); ++i) {
            r[i] = prob.b[i];
        }
        for (int i = 0; i < (int) v.size(); ++i) {
            for (auto p: prob.A[i]) {
                r[p.first] -= p.second * v[i];
            }
        }

        // update z, yl
        double sum2 = 0;
        for (int i = 0; i < (int) z.size(); ++i) {
            z[i] = std::max(0.0, r[i] - yl[i]);
            yl[i] = std::max(0.0, yl[i] - r[i]);
            sum2 += (z[i] - r[i]) * (z[i] - r[i]);
        }

        if (sum2 < eps_stop) {
            break;
        }
    }

    bool answer = true;
    TCodeword result((int) y.size());
    for (int i = 0; i < (int) y.size(); ++i) {
        auto val = v[i];
        if (val <= 0.5) {
            result[i] = false;
        } else {
            result[i] = true;
        }
        if (val >= EPS && val <= 1.0 - EPS) {
            answer = false;
        }
    }

    return {result, answer};
}

class QPADMMDecoder : public Decoder {
public:
    explicit QPADMMDecoder(double alpha, double mu, int max_iter = 2000, double eps_stop = 1e-5) :
            _alpha(alpha), _mu(mu), _max_iter(max_iter), _eps_stop(eps_stop) {}

    pair<TCodeword, bool> decode(const TMatrix &H, const TFVector &channel_word, double snr) override {
        return DecodeQPADMM(H, channel_word, snr, _alpha, _mu, _max_iter, _eps_stop);
    }

    string name() const override { return "QP-ADMM"; }

private:
    double _alpha, _mu, _eps_stop;
    int _max_iter;
};

#endif //ACG_ALP_LDPC_QP_ADMM_H
