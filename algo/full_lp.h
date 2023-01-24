#ifndef ACG_ALP_LDPC_FULL_LP_H
#define ACG_ALP_LDPC_FULL_LP_H

#include "algo/algo.h"
#include "glpk.h"

void AddThreeVariableInnequalities(glp_prob *lp, int i, int j, int h) {
    glp_add_rows(lp, 4);
    glp_set_row_bnds(lp, glp_get_num_rows(lp) - 3, GLP_UP, 0, 0);
    glp_set_row_bnds(lp, glp_get_num_rows(lp) - 2, GLP_UP, 0, 0);
    glp_set_row_bnds(lp, glp_get_num_rows(lp) - 1, GLP_UP, 0, 0);
    glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_UP, 0, 2.0);
    std::vector<int> row_idx(4, 0);
    row_idx[1] = i;
    row_idx[2] = j;
    row_idx[3] = h;
    std::vector<double> row_coef(4, 0.0);
    {
        row_coef[1] = 1.0;
        row_coef[2] = -1.0;
        row_coef[3] = -1.0;
        glp_set_mat_row(lp, glp_get_num_rows(lp) - 3, 3, row_idx.data(), row_coef.data());
    }
    {
        row_coef[1] = -1.0;
        row_coef[2] = 1.0;
        row_coef[3] = -1.0;
        glp_set_mat_row(lp, glp_get_num_rows(lp) - 2, 3, row_idx.data(), row_coef.data());
    }
    {
        row_coef[1] = -1.0;
        row_coef[2] = -1.0;
        row_coef[3] = 1.0;
        glp_set_mat_row(lp, glp_get_num_rows(lp) - 1, 3, row_idx.data(), row_coef.data());
    }
    {
        row_coef[1] = 1.0;
        row_coef[2] = 1.0;
        row_coef[3] = 1.0;
        glp_set_mat_row(lp, glp_get_num_rows(lp), 3, row_idx.data(), row_coef.data());
    }
}

pair<TCodeword, bool> DecodeFromLp(glp_prob *lp, const TFVector &y) {
    bool answer = true;
    TCodeword result((int) y.size());
    for (int i = 0; i < (int) y.size(); ++i) {
        auto val = glp_get_col_prim(lp, i + 1);
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

pair<TCodeword, bool> DecodeFullLP(const std::vector<TCodeword> &H, const std::vector<double> &y, double snr) {
    glp_term_out(GLP_MSG_OFF);

    glp_prob *lp = glp_create_prob();

    glp_set_obj_dir(lp, GLP_MIN);

    int n_aux = 0;
    for (int i = 0; i < (int) H.size(); ++i) {
        int sz = 0;
        for (int t = 0; t < (int) H[i].size(); ++t) {
            sz += H[i][t];
        }
        n_aux += std::max(sz - 3, 0);
    }

    auto coef = CalculateCoef(y, snr);
    glp_add_cols(lp, y.size() + n_aux);
    for (int i = 0; i < y.size(); ++i) {
        glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(lp, i + 1, coef[i]);
    }
    for (int i = 0; i < n_aux; ++i) {
        glp_set_col_bnds(lp, y.size() + i + 1, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(lp, y.size() + i + 1, 0.0);
    }

    int pos = y.size();
    for (const auto &i: H) {
        std::vector<int> idx;
        for (int j = 0; j < (int) y.size(); ++j) {
            if (i[j]) {
                idx.emplace_back(j);
            }
        }
        if (idx.empty()) {
            continue;
        }
        if (idx.size() == 1) {
            glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_UP, 0, 0);
            {
                std::vector<int> row_idx(2, 0);
                std::vector<double> row_coef(2, 0.0);
                row_idx[1] = idx[0] + 1;
                row_coef[1] = 1.0;
                glp_set_mat_row(lp, glp_get_num_rows(lp), 1, row_idx.data(), row_coef.data());
            }
            continue;
        }
        if (idx.size() == 2) {
            glp_add_rows(lp, 2);
            glp_set_row_bnds(lp, glp_get_num_rows(lp) - 1, GLP_UP, 0, 0);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_UP, 0, 0);
            std::vector<int> row_idx(3, 0);
            row_idx[1] = idx[0] + 1;
            row_idx[2] = idx[1] + 1;
            {
                std::vector<double> row_coef(3, 0.0);
                row_coef[1] = 1.0;
                row_coef[2] = -1.0;
                glp_set_mat_row(lp, glp_get_num_rows(lp) - 1, 2, row_idx.data(), row_coef.data());
            }
            {
                std::vector<double> row_coef(3, 0.0);
                row_coef[1] = -1.0;
                row_coef[2] = 1.0;
                glp_set_mat_row(lp, glp_get_num_rows(lp), 2, row_idx.data(), row_coef.data());
            }
            continue;
        }
        int idx_last = idx[0] + 1;
        for (int j = 1; j < (int) idx.size() - 2; ++j) {
            int idx_mid = idx[j] + 1;
            int idx_aux = ++pos;
            AddThreeVariableInnequalities(lp, idx_last, idx_mid, idx_aux);
            idx_last = idx_aux;
        }
        AddThreeVariableInnequalities(lp, idx_last, idx[(int) idx.size() - 2] + 1, idx.back() + 1);
    }

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.meth = GLP_DUALP;
    glp_simplex(lp, &parm);

    pair<TCodeword, bool> res = DecodeFromLp(lp, y);

    glp_delete_prob(lp);

    if (res.second) {
        assert(IsCodeword(H, res.first));
    }

    return res;
}

class FullLPDecoder : public Decoder {
public:
    FullLPDecoder() {}

    pair<TCodeword, bool> decode(const TMatrix &H, const TFVector &channel_word, double snr) override {
        return DecodeFullLP(H, channel_word, snr);
    }

    string name() const override { return "FullLP"; }
};

#endif //ACG_ALP_LDPC_FULL_LP_H
