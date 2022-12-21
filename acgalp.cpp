#include "glpk.h"

#include <bits/stdc++.h>

constexpr int N = 128;
constexpr double EPS = 1e-6;

using TCodeword = std::bitset<N>;
#define TIME (clock() * 1.0 / CLOCKS_PER_SEC)

std::istream& operator>>(std::istream& in, TCodeword& res) {
    std::string s;
    in >> s;
    for (int i = 0; i < N; ++i) {
        res[i] = (bool)(s[i] - '0');
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, const TCodeword& c) {
    for (int i = 0; i < N; i++) {
        out << c[i];
    }
    return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vc) {
    for (const auto& t : vc) {
        out << t << "\n";
    }
    return out;
}

std::vector<TCodeword> ReadCodewords(const std::string& filename) {
    std::ifstream fin(filename);
    int cnt;
    fin >> cnt;
    std::vector<TCodeword> res(cnt);
    for (int i = 0; i < cnt; ++i) {
        fin >> res[i];
    }
    return res;
}

std::vector<TCodeword> ReadParityCheckMatrix(const std::string& filename) {
    std::ifstream fin(filename);
    int n, m;
    fin >> n >> m;
    assert(n == N);
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
    std::vector<TCodeword> result(m, 0);
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

bool IsCodeword(const std::vector<TCodeword>& H, const TCodeword& c) {
    for (const auto& vec : H) {
        if ((vec & c).count() % 2 == 1) {
            return false;
        }
    }
    return true;
}

std::vector<TCodeword> GetOrtogonal(std::vector<TCodeword> H) {
    std::vector<int> pos(H.size());
    TCodeword is_main;
    for (int i = 0; i < (int)H.size(); ++i) {
        pos[i] = -1;
        for (int j = 0; j < N; ++j) {
            if (H[i][j]) {
                pos[i] = j;
                break;
            }
        }
        assert(pos[i] != -1);
        for (int k = 0; k < (int)H.size(); ++k) {
            if (k != i && H[k][pos[i]]) {
                H[k] ^= H[i];
            }
        }
        is_main[pos[i]] = true;
    }
    std::vector<TCodeword> res(N - (int)H.size());
    int idx = 0;
    for (int j = 0; j < N; ++j) {
        if (!is_main[j]) {
            res[idx][j] = true;
            for (int i = 0; i < (int)H.size(); ++i) {
                if (H[i][j]) {
                    res[idx][pos[i]] = true;
                }
            }
            ++idx;
        }
    }
    return res;
}

template <typename Gen>
TCodeword GenerateRandomCodeword(const std::vector<TCodeword>& G, Gen& rnd) {
    TCodeword res;
    for (const auto& vec : G) {
        if (rnd() % 2 == 0) {
            res ^= vec;
        }
    }
    return res;
}


std::vector<double> GetSolution(glp_prob* lp) {
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = glp_get_col_prim(lp, i + 1);
    }
    return res;
}

int AddRowsALP(const std::vector<TCodeword>& H, glp_prob* lp) {
    auto u = GetSolution(lp);
    int added_rows = 0;
    for (int i = 0; i < (int)H.size(); ++i) {
        int n_size = 0;
        int v_size = 0;
        int j_best = 0;
        double val_best = 10;
        for (int j = 0; j < N; ++j) {
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
        TCodeword is_v;
        for (int j = 0; j < N; ++j) {
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
        for (int j = 0; j < N; ++j) {
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
            for (int j = 0; j < N; ++j) {
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

            glp_set_mat_row(lp, glp_get_num_rows(lp), n_size, idx.data(),coef.data());
        }
    }
    return added_rows;
}

bool DecodeALP(const std::vector<TCodeword>& H, TCodeword y, TCodeword& result) {
    glp_prob *lp = glp_create_prob();

    glp_set_obj_dir(lp, GLP_MIN);

    glp_add_cols(lp, N);
    for (int i = 0; i < N; ++i) {
        glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(lp, i + 1, y[i] ? -1.0 : 1.0);
    }

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.meth = GLP_DUALP;
    glp_simplex(lp, &parm);

    while (AddRowsALP(H, lp) != 0) {
        glp_simplex(lp, &parm);
    }

    bool answer = true;
    result = 0;
    for (int i = 0; i < N; ++i) {
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

std::vector<TCodeword> ShuffleColumns(const std::vector<TCodeword>& H, const std::vector<int>& p) {
    std::vector<TCodeword> res(H.size(), 0);
    for (int i = 0; i < (int)H.size(); ++i) {
        for (int j = 0; j < N; ++j) {
            res[i][j] = H[i][p[j]];
        }
    }
    return res;
}

std::vector<TCodeword> CalculateGauss(const std::vector<TCodeword>& H0, const std::vector<double>& u) {
    std::vector<int> non_int, zeros, ones;
    for (int i = 0; i < N; ++i) {
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
    for (int i : zeros) {
        p.emplace_back(i);
    }
    for (int i : ones) {
        p.emplace_back(i);
    }
    std::vector<int> p_inv(N);
    for (int i = 0; i < N; ++i) {
        p_inv[p[i]] = i;
    }
    int col = 0;
    auto H = ShuffleColumns(H0, p);
    std::vector<int> pos(H.size());
    for (int i = 0; i < (int)H.size(); ++i) {
        while (col < N) {
            bool found = false;
            for (int t = i; t < (int)H.size(); ++t) {
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
        assert(col < N);
        pos[i] = col;
        ++col;
        for (int k = 0; k < (int)H.size(); ++k) {
            if (k != i && H[k][pos[i]]) {
                H[k] ^= H[i];
            }
        }
    }
    return ShuffleColumns(H, p_inv);
}

bool DecodeAGCALP(const std::vector<TCodeword>& H, TCodeword y, TCodeword& result) {
    glp_prob *lp = glp_create_prob();

    glp_set_obj_dir(lp, GLP_MIN);

    glp_add_cols(lp, N);
    for (int i = 0; i < N; ++i) {
        glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(lp, i + 1, y[i] ? -1.0 : 1.0);
    }

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.meth = GLP_DUALP;
    glp_simplex(lp, &parm);

    /*std::vector<std::vector<TCodeword>> mats;
    mats.emplace_back(H);
    while (mats.size() < 2) {
        bool found = false;
        for (const auto& mat : mats) {
            if (AddRowsALP(mat, lp) != 0) {
                found = true;
                break;
            }
        }
        //std::cerr << found << "\n";
        if (!found) {
            auto H_cur = CalculateGauss(H, GetSolution(lp));
            if (AddRowsALP(H_cur, lp) == 0) {
                break;
            }
            mats.emplace_back(H_cur);
        }
    }*/
    while (AddRowsALP(H, lp) != 0 || AddRowsALP(CalculateGauss(H, GetSolution(lp)), lp) != 0) {
        glp_simplex(lp, &parm);
    }

    bool answer = true;
    result = 0;
    for (int i = 0; i < N; ++i) {
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

    return IsCodeword(H, result);

    /*if (answer) {
        assert(IsCodeword(H, result));
    }

    return answer;*/
}

template <typename Gen>
TCodeword Transmit(double std, const TCodeword& c, Gen& rnd) {
    std::normal_distribution<double> dst(0.0, std);
    TCodeword res;
    for (int i = 0; i < N; ++i) {
        double x = c[i] ? 1.0 : -1.0;
        x += dst(rnd);
        if (x > 0) {
            res[i] = true;
        }
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
        std::cout << "Success percent " << ((double)correct / total) * 100 << "%" << std::endl;
    }
};

template <typename Func>
ExperimentResult MakeExperiment(
    Func decoding_func,
    double std,
    int seed,
    const std::vector<TCodeword>& H,
    const std::vector<TCodeword>& tests,
    int n_iter = -1
) {
    ExperimentResult result{0, 0, 0, 0};
    std::mt19937 rnd(seed);
    for (const auto& c : tests) {
        ++result.total;
        auto y = Transmit(std, c, rnd);
        TCodeword res;
        int cd = 0;
        if (decoding_func(H, y, res)) {
            if (res == c) {
                cd = 1;
                ++result.correct;
            } else {
                ++result.wrong;
            }
        } else {
            ++result.pseudo;
        }
        std::cout << result.total << ": " << cd << ", hamming: " << (c ^ y).count() << std::endl;
        if (result.pseudo + result.wrong == n_iter) {
            break;
        }
    }
    return result;
}

int main() {
#ifdef ONPC
    freopen("input", "r", stdin);
#endif
    std::ios::sync_with_stdio(0); std::cin.tie(0); std::cout.tie(0);

    glp_term_out(GLP_MSG_OFF);

    auto H = ReadParityCheckMatrix("H.txt");
    auto G = GetOrtogonal(H);
    for (const auto& vec : G) {
        assert(IsCodeword(H, vec));
    }

    auto tests = ReadCodewords("codewords.txt");
    std::cerr << "There are " << tests.size() << " codewords" << std::endl;

    for (const auto& vec : tests) {
        assert(IsCodeword(H, vec));
    }

    auto result = MakeExperiment(DecodeAGCALP, 0.6, 239, H, tests, 50);
    result.Print();
    std::cerr << TIME << "\n";

    /*std::mt19937 rnd(239);
    std::ofstream fout("codewords.txt");
    for (int i = 0; i < 10000; ++i) {
        auto c = GenerateRandomCodeword(G, rnd);
        fout << c << "\n";
    }
    return 0;*/

    /*for (int it = 0; it < 100; ++it) {
        assert(IsCodeword(H, GenerateRandomCodeword(G, rnd)));
    }*/

    /*std::vector<int> vals{1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1,
                          1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0,
                          0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    TCodeword some;
    for (int i = 0; i < N; i++) {
        some[i] = (bool)vals[i];
    }
    std::cout << IsCodeword(H, some) << "\n";
    auto modif = some;
    modif[N - 10] = 1 - modif[N - 10];
    modif[N - 20] = 1 - modif[N - 20];
    modif[N - 30] = 1 - modif[N - 30];
    modif[N - 40] = 1 - modif[N - 40];
    modif[N - 50] = 1 - modif[N - 50];
    modif[N - 1] = 1 - modif[N - 1];
    modif[N - 2] = 1 - modif[N - 2];
    modif[N - 3] = 1 - modif[N - 3];
    modif[N - 4] = 1 - modif[N - 4];
    modif[N - 5] = 1 - modif[N - 5];
    modif[N - 6] = 1 - modif[N - 6];
    modif[N - 7] = 1 - modif[N - 7];
    modif[N - 8] = 1 - modif[N - 8];
    modif[N - 9] = 1 - modif[N - 9];
    modif[N - 10] = 1 - modif[N - 10];
    std::cout << modif << "\n";
    std::cout << IsCodeword(H, modif) << "!\n";

    TCodeword res;
    std::cout << DecodeAGCALP(H, modif, res) << "\n";
    std::cout << modif << "\n";
    std::cout << res << "\n";
    std::cout << some << "\n";
    std::cout << IsCodeword(H, res) << "??\n";
    assert(res == some);*/

    return 0;
}