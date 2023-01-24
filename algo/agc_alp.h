#ifndef ACG_ALP_LDPC_AGC_ALP_H
#define ACG_ALP_LDPC_AGC_ALP_H

#include "algo/alp.h"
#include "algo/full_lp.h"
#include "glpk.h"


std::vector <TCodeword> ShuffleColumns(const std::vector <TCodeword> &H, const std::vector<int> &p) {
    std::vector <TCodeword> res(H.size(), TCodeword(p.size(), false));
    for (int i = 0; i < (int) H.size(); ++i) {
        for (int j = 0; j < (int) p.size(); ++j) {
            res[i][j] = H[i][p[j]];
        }
    }
    return res;
}

std::vector <TCodeword> CalculateGauss(const std::vector <TCodeword> &H0, const std::vector<double> &u) {
    std::vector<int> non_int, zeros, ones;
    for (int i = 0; i < (int) u.size(); ++i) {
        if (u[i] < EPS) {
            zeros.emplace_back(i);
        } else if (u[i] > 1 - EPS) {
            ones.emplace_back(i);
        } else {
            non_int.emplace_back(i);
        }
    }
    std::stable_sort(non_int.begin(), non_int.end(), [&](int i, int j) {
        return abs(u[i] - 0.5) < abs(u[j] - 0.5);
    });
    std::vector<int> p = non_int;
    for (int i: zeros) {
        p.emplace_back(i);
    }
    for (int i: ones) {
        p.emplace_back(i);
    }
    std::vector<int> p_inv((int) u.size());
    for (int i = 0; i < (int) u.size(); ++i) {
        p_inv[p[i]] = i;
    }
    int col = 0;
    auto H = ShuffleColumns(H0, p);
    std::vector<int> pos(H.size());
    for (int i = 0; i < (int) H.size(); ++i) {
        while (col < (int) u.size()) {
            bool found = false;
            for (int t = i; t < (int) H.size(); ++t) {
                if (H[t][col]) {
                    std::swap(H[i], H[t]);
                    found = true;
                    break;
                }
            }
            if (found) {
                break;
            }
            ++col;
        }
        assert(col < (int) u.size());
        pos[i] = col;
        ++col;
        for (int k = 0; k < (int) H.size(); ++k) {
            if (k != i && H[k][pos[i]]) {
                for (int t = 0; t < (int) H[k].size(); ++t) {
                    H[k][t] = H[k][t] ^ H[i][t];
                }
            }
        }
    }
    return ShuffleColumns(H, p_inv);
}

class AGCALPDecoder : public Decoder {
public:
    explicit AGCALPDecoder(int max_rows) : _max_rows(max_rows) {}

    pair<TCodeword, bool> decode(const TMatrix &H, const TFVector &y, double snr) override {
        glp_term_out(GLP_MSG_OFF);

        glp_prob *lp = glp_create_prob();

        glp_set_obj_dir(lp, GLP_MIN);

        auto coef = CalculateCoef(y, snr);
        glp_add_cols(lp, y.size());
        for (int i = 0; i < y.size(); ++i) {
            glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
            glp_set_obj_coef(lp, i + 1, coef[i]);
        }

        glp_smcp parm;
        glp_init_smcp(&parm);
        parm.meth = GLP_DUALP;
        glp_simplex(lp, &parm);

        while (glp_get_num_rows(lp) < _max_rows &&
               (AddRowsALP(H, lp) != 0 || AddRowsALP(CalculateGauss(H, GetSolution(lp, y.size())), lp) != 0)) {
            glp_simplex(lp, &parm);
        }

        pair<TCodeword, bool> res = DecodeFromLp(lp, y);

        glp_delete_prob(lp);

        if (res.second) {
            assert(IsCodeword(H, res.first));
        }

        return res;
    }

    string name() const override { return "AGC-ALP"; }

private:
    int _max_rows;
};

#endif //ACG_ALP_LDPC_AGC_ALP_H
