#include "glpk.h"

#include <bits/stdc++.h>

//constexpr int N = 128;
constexpr double EPS = 1e-6;

using TCodeword = std::vector<bool>;
using TDecoderFunction = std::function<bool(const std::vector<TCodeword> &, const std::vector<double> &, TCodeword &,
                                            double)>;
#define TIME (clock() * 1.0 / CLOCKS_PER_SEC)

std::istream &operator>>(std::istream &in, TCodeword &res) {
    std::string s;
    in >> s;
    res.resize(s.size());
    for (int i = 0; i < (int) res.size(); ++i) {
        res[i] = (bool) (s[i] - '0');
    }
    return in;
}

std::ostream &operator<<(std::ostream &out, const TCodeword &c) {
    for (int i = 0; i < (int) c.size(); i++) {
        out << c[i];
    }
    return out;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &vc) {
    for (const auto &t: vc) {
        out << t << "\n";
    }
    return out;
}

double GetllrSigma(double snr) {
    return sqrt(8 * 0.5 * pow(10, (snr / 10)));
}

double GetllrMean(double snr) {
    return 4 * 0.5 * pow(10, (snr / 10));
}

std::vector<TCodeword> ReadCodewords(const std::string &filename) {
    std::ifstream fin(filename);
    int cnt;
    fin >> cnt;
    std::vector<TCodeword> res(cnt);
    for (int i = 0; i < cnt; ++i) {
        fin >> res[i];
    }
    return res;
}

std::vector<TCodeword> ReadParityCheckMatrix(const std::string &filename, bool is_raw = false) {
    std::ifstream fin(filename);
    if (is_raw) {
        int m = -1;
        std::string line;
        std::vector<TCodeword> H;
        int max_cnt = 0;
        while (getline(fin, line)) {
            if (line.empty()) {
                continue;
            }
            if (m == -1) {
                m = ((int) line.size() + 1) / 2;
            } else {
                assert(m == ((int) line.size() + 1) / 2);
            }
            H.emplace_back(m, false);
            int cnt = 0;
            for (int i = 0; i < line.size(); i += 2) {
                if (line[i] == '1') {
                    H.back()[i / 2] = true;
                    cnt++;
                }
            }
            max_cnt = std::max(max_cnt, cnt);
            assert(cnt != 0);
        }
        std::cerr << "n=" << H.size() << " " << "m=" << m << " " << "r=" << max_cnt << std::endl;
        return H;
    } else {
        int n, m;
        fin >> n >> m;
        //assert(n == N);
        int dc, dr;
        fin >> dc >> dr;
        std::vector<int> d_col(n);
        for (int i = 0; i < n; ++i) {
            fin >> d_col[i];
        }
        std::vector<int> d_row(m);
        for (int i = 0; i < m; ++i) {
            fin >> d_row[i];
        }
        std::vector<TCodeword> result(m, TCodeword(n, false));
        for (int i = 0; i < n; ++i) {
            std::vector<int> id(dc);
            for (int j = 0; j < dc; ++j) {
                fin >> id[j];
                --id[j];
            }
            for (int j = 0; j < dc; ++j) {
                if (id[j] >= 0) {
                    result[id[j]][i] = true;
                }
            }
        }
        return result;
    }
}

bool IsCodeword(const std::vector<TCodeword> &H, const TCodeword &c) {
    for (const auto &vec: H) {
        bool parity = false;
        for (int i = 0; i < (int) vec.size(); ++i) {
            parity = parity ^ (vec[i] && c[i]);
        }
        if (parity) {
            return false;
        }
    }
    return true;
}

std::vector<TCodeword> GetOrtogonal(std::vector<TCodeword> H) {
    assert(!H.empty());
    std::vector<int> pos(H.size());
    TCodeword is_main(H[0].size(), false);
    for (int i = 0; i < (int) H.size(); ++i) {
        pos[i] = -1;
        for (int j = 0; j < (int) H[i].size(); ++j) {
            if (H[i][j]) {
                pos[i] = j;
                break;
            }
        }
        assert(pos[i] != -1);
        for (int k = 0; k < (int) H.size(); ++k) {
            if (k != i && H[k][pos[i]]) {
                for (int t = 0; t < (int) H[i].size(); ++t) {
                    H[k][t] = H[k][t] ^ H[i][t];
                }
            }
        }
        is_main[pos[i]] = true;
    }
    std::vector<TCodeword> res(H[0].size() - H.size(), TCodeword(H[0].size(), false));
    int idx = 0;
    for (int j = 0; j < (int) H[0].size(); ++j) {
        if (!is_main[j]) {
            res[idx][j] = true;
            for (int i = 0; i < (int) H.size(); ++i) {
                if (H[i][j]) {
                    res[idx][pos[i]] = true;
                }
            }
            ++idx;
        }
    }
    return res;
}

template<typename Gen>
TCodeword GenerateRandomCodeword(const std::vector<TCodeword> &G, Gen &rnd) {
    assert(!G.empty());
    TCodeword res(G[0].size(), false);
    for (const auto &vec: G) {
        if (rnd() % 2 == 0) {
            for (int t = 0; t < (int) vec.size(); ++t) {
                res[t] = res[t] ^ vec[t];
            }
        }
    }
    return res;
}


std::vector<double> GetSolution(glp_prob *lp, int N) {
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = glp_get_col_prim(lp, i + 1);
    }
    return res;
}

int AddRowsALP(const std::vector<TCodeword> &H, glp_prob *lp) {
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

double CalculateLogDensity(double mean, double std, double value) {
    return -(value - mean) * (value - mean);
}

std::vector<double> CalculateCoef(const std::vector<double> &y, double snr) {
    double llrMean = GetllrMean(snr);
    double llrSigma = GetllrSigma(snr);

    std::vector<double> coef(y.size());
    for (int i = 0; i < (int) y.size(); ++i) {
        coef[i] = CalculateLogDensity(llrMean, llrSigma, -y[i]) - CalculateLogDensity(llrMean, llrSigma, y[i]);
    }

    return coef;
}

bool DecodeALP(const std::vector<TCodeword> &H, const std::vector<double> &y, TCodeword &result, double snr) {
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

    bool answer = true;
    result.assign(y.size(), false);
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

    glp_delete_prob(lp);

    if (answer) {
        assert(IsCodeword(H, result));
    }

    return answer;
}

std::vector<TCodeword> ShuffleColumns(const std::vector<TCodeword> &H, const std::vector<int> &p) {
    std::vector<TCodeword> res(H.size(), TCodeword(p.size(), false));
    for (int i = 0; i < (int) H.size(); ++i) {
        for (int j = 0; j < (int) p.size(); ++j) {
            res[i][j] = H[i][p[j]];
        }
    }
    return res;
}

std::vector<TCodeword> CalculateGauss(const std::vector<TCodeword> &H0, const std::vector<double> &u) {
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

bool DecodeAGCALP(const std::vector<TCodeword> &H, const std::vector<double> &y, TCodeword &result, double snr,
                  int max_rows = 2000) {
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

    while (glp_get_num_rows(lp) < max_rows &&
           (AddRowsALP(H, lp) != 0 || AddRowsALP(CalculateGauss(H, GetSolution(lp, y.size())), lp) != 0)) {
        glp_simplex(lp, &parm);
    }

    bool answer = true;
    result.assign(y.size(), false);
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

    glp_delete_prob(lp);

    if (answer) {
        assert(IsCodeword(H, result));
    }
    return answer;
}

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

bool DecodeFullLP(const std::vector<TCodeword> &H, const std::vector<double> &y, TCodeword &result, double snr) {
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

    bool answer = true;
    result.assign(y.size(), false);
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
    glp_delete_prob(lp);

    if (answer) {
        assert(IsCodeword(H, result));
    }
    return answer;
}

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
        for (auto [_, cf]: A[i]) {
            e[i] += cf * cf;
        }
    }

    return ADMMProblem{q, A, b, e};
}

bool DecodeQPADMM(const std::vector<TCodeword> &H, const std::vector<double> &y, TCodeword &result, double snr,
                  double alpha, double mu, int max_iter = 2000, double eps_stop = 1e-5) {
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
            for (auto [j, cf]: prob.A[i]) {
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
            for (auto [j, cf]: prob.A[i]) {
                r[j] -= cf * v[i];
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
    result.assign(y.size(), false);
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

    return answer;
}

template<typename Gen>
std::vector<double> Transmit(double snr, const TCodeword &c, Gen &rnd) {
    double llrMean = GetllrMean(snr);
    double llrSigma = GetllrSigma(snr);
    std::vector<double> res(c.size());
    std::normal_distribution<double> dst(llrMean, llrSigma);
    for (int i = 0; i < (int) c.size(); ++i) {
        res[i] = (c[i] ? 1.0 : -1.0) * dst(rnd);
    }
    return res;
}

struct ExperimentResult {
    int total;
    int correct;
    int wrong;
    int pseudo;

    void Print() const {
        std::cout << correct << "/" << total << " correct, ";
        std::cout << wrong << "/" << total << " wrong, ";
        std::cout << pseudo << "/" << total << " pseudocodewords" << std::endl;
        std::cout << "Success percent " << ((double) correct / total) * 100 << "%" << std::endl;
        std::cout << "FER: " << ((double) (wrong + pseudo) / total) << std::endl;
    }

    double GetFER() const {
        return ((double) (wrong + pseudo) / total);
    }
};

ExperimentResult MakeExperiment(
        const TDecoderFunction &decoding_func,
        double snr,
        int seed,
        const std::vector<TCodeword> &H,
        const std::vector<TCodeword> &tests,
        int n_iter = -1,
        bool full_verbose = true
) {
    assert(!H.empty());
    ExperimentResult result{0, 0, 0, 0};
    std::mt19937 rnd(seed);
    int sum_hamming = 0;
    int sum_hamming_ok = 0;
    int sum_hamming_wrong = 0;
    for (const auto &c: tests) {
        ++result.total;
        auto y = Transmit(snr, c, rnd);
        int hamming = 0;
        for (int i = 0; i < (int) H[0].size(); ++i) {
            if (c[i] && y[i] <= 0) {
                hamming++;
            }
            if (!c[i] && y[i] > 0) {
                hamming++;
            }
        }
        sum_hamming += hamming;
        TCodeword res(H[0].size(), false);
        int cd = 0;
        if (decoding_func(H, y, res, snr)) {
            if (res == c) {
                cd = 1;
                ++result.correct;
                sum_hamming_ok += hamming;
            } else {
                ++result.wrong;
                sum_hamming_wrong += hamming;
            }
        } else {
            ++result.pseudo;
            sum_hamming_wrong += hamming;
        }
        if (full_verbose || result.total % 1000 == 0) {
            std::cout << result.total << ": " << result.correct << ", hamming: " << hamming << std::endl;
        }
        if (result.pseudo + result.wrong == n_iter) {
            break;
        }
    }
    std::cout << "Average hamming distance after transmittion: " << (double) sum_hamming / result.total << std::endl;
    std::cout << "Average hamming distance of correctly decoded codes: "
              << (result.correct == 0 ? -1 : (double) sum_hamming_ok / result.correct) << std::endl;
    std::cout << "Average hamming distance of incorrectly decoded codes: "
              << (result.total == result.correct ? -1 : (double) sum_hamming_wrong / (result.total - result.correct))
              << std::endl;
    return result;
}

struct ExperimentReport {
    std::string method;
    double snr, fer, time;
};

int main() {
    std::ios::sync_with_stdio(0);

    glp_term_out(GLP_MSG_OFF);

    auto H = ReadParityCheckMatrix("data/H05.txt", true);
    auto G = ReadParityCheckMatrix("data/G05.txt", true);
    /*auto G = GetOrtogonal(H);*/
    for (const auto &vec: G) {
        assert(IsCodeword(H, vec));
    }

    //auto tests = ReadCodewords("codewords.txt");
    int tests_count = 10000;
    std::mt19937 rnd(239);
    std::vector<TCodeword> tests;
    tests.reserve(tests_count);
    for (int i = 0; i < tests_count; ++i) {
        tests.emplace_back(GenerateRandomCodeword(G, rnd));
    }
    std::cerr << "There are " << tests.size() << " codewords" << std::endl;

    for (const auto &vec: tests) {
        assert(IsCodeword(H, vec));
    }

    std::map<std::string, std::vector<ExperimentReport>> info;
    std::vector<double> snr_list = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<std::string> method_list = {"AGC-ALP", "ALP", "Full-LP", "QP-ADMM"};
    for (const std::string &name: method_list) {
        std::cout << "Decoding method: " << name << std::endl;
        TDecoderFunction decoder;
        if (name == "AGC-ALP") {
            decoder = [](const std::vector<TCodeword> &H, const std::vector<double> &y, TCodeword &result, double snr) {
                return DecodeAGCALP(H, y, result, snr);
            };
        } else if (name == "ALP") {
            decoder = DecodeALP;
        } else if (name == "Full-LP") {
            decoder = DecodeFullLP;
        } else if (name == "QP-ADMM") {
            decoder = [](const std::vector<TCodeword> &H, const std::vector<double> &y, TCodeword &result, double snr) {
                return DecodeQPADMM(H, y, result, snr, 0.6, 1.0, 1000);
            };
        } else {
            assert(false);
        }
        for (double snr: snr_list) {
            std::cout << "SNR=" << snr << std::endl;
            double start_time = TIME;
            auto result = MakeExperiment(decoder, snr, 239, H, tests, 100, false);
            double working_time = (TIME - start_time) / result.total;
            result.Print();
            info[name].emplace_back(ExperimentReport{name, snr, result.GetFER(), working_time});
        }
        std::cout << "______________________________________________________________________________" << std::endl;
    }

    std::cout << std::fixed << std::setprecision(6);
    for (const auto &base: info) {
        std::cout << std::string(10, '_') << " " << base.first << " " << std::string(10, '_') << std::endl;
        for (const auto &exp: base.second) {
            std::cout << "snr=" << exp.snr << ": fer=" << exp.fer << " (time=" << exp.time << "s)" << std::endl;
        }
    }

    return 0;
}
