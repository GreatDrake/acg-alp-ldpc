#ifndef ACG_ALP_LDPC_ALP_H
#define ACG_ALP_LDPC_ALP_H

#include "glpk.h"

#include "../utils/channel.h"
#include "../utils/codeword.h"
#include "../utils/parse_data.h"
#include "algo/algo.h"
#include "algo/full_lp.h"


std::vector<double> GetSolution(glp_prob *lp, int N) {
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = glp_get_col_prim(lp, i + 1);
    }
    return res;
}

int AddRowsALP(const std::vector <TCodeword> &H, glp_prob *lp) {
    auto u = GetSolution(lp, H[0].size());
    int added_rows = 0;
    for (int i = 0; i < (int) H.size(); ++i) {
        int n_size = 0;
        int v_size = 0;
        int j_best = 0;
        double val_best = 10;
        for (int j = 0; j < (int) H[i].size(); ++j) {
            if (H[i][j]) {
                ++n_size;
                double current = abs(u[j] - 0.5);
                if (current < val_best) {
                    j_best = j;
                    val_best = current;
                }
                if (u[j] > 0.5) {
                    v_size++;
                }
            }
        }
        if (n_size == 0) {
            continue;
        }
        int v_final_size = 0;
        TCodeword is_v(H[i].size(), false);
        for (int j = 0; j < (int) H[i].size(); ++j) {
            if (H[i][j]) {
                if (j == j_best && v_size % 2 == 0) {
                    if (u[j] <= 0.5) {
                        ++v_final_size;
                        is_v[j] = true;
                    }
                } else {
                    if (u[j] > 0.5) {
                        ++v_final_size;
                        is_v[j] = true;
                    }
                }
            }
        }
        assert(v_final_size % 2 == 1);
        double sum = 0;
        for (int j = 0; j < (int) H[i].size(); ++j) {
            if (H[i][j]) {
                if (is_v[j]) {
                    sum += 1 - u[j];
                } else {
                    sum += u[j];
                }
            }
        }
        if (sum < 1.0 - EPS) {
            glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_UP, 0, v_final_size - 1);

            ++added_rows;
            std::vector<int> idx(n_size + 1, 0);
            std::vector<double> coef(n_size + 1, 0.0);
            int pos = 1;
            for (int j = 0; j < (int) H[i].size(); ++j) {
                if (H[i][j]) {
                    idx[pos] = j + 1;
                    if (is_v[j]) {
                        coef[pos] = 1.0;
                    } else {
                        coef[pos] = -1.0;
                    }
                    ++pos;
                }
            }

            glp_set_mat_row(lp, glp_get_num_rows(lp), n_size, idx.data(), coef.data());
        }
    }
    return added_rows;
}

class ALPDecoder : public Decoder {
public:
    ALPDecoder() {}

    pair<TCodeword, bool> decode(const TMatrix &H, const TFVector &y, double snr) override {
        glp_term_out(GLP_MSG_OFF);

        glp_prob *lp = glp_create_prob();

        glp_set_obj_dir(lp, GLP_MIN);

        auto coef = CalculateCoef(y, snr);
        glp_add_cols(lp, y.size());
        for (int i = 0; i < (int) y.size(); ++i) {
            glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
            glp_set_obj_coef(lp, i + 1, coef[i]);
        }

        glp_smcp parm;
        glp_init_smcp(&parm);
        parm.meth = GLP_DUALP;
        glp_simplex(lp, &parm);

        while (AddRowsALP(H, lp) != 0) {
            glp_simplex(lp, &parm);
        }

        pair<TCodeword, bool> res = DecodeFromLp(lp, y);

        glp_delete_prob(lp);

        if (res.second) {
            assert(IsCodeword(H, res.first));
        }

        return res;
    }

    string name() const override { return "ALP"; }
};

#endif //ACG_ALP_LDPC_ALP_H
