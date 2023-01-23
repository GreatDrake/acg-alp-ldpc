#ifndef ACG_ALP_LDPC_CODEWORD_H
#define ACG_ALP_LDPC_CODEWORD_H

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cassert>
#include <cassert>
#include <unordered_map>
#include <utility>
#include <vector>
#include <fstream>
#include <random>

using namespace std;

typedef vector<bool> TCodeword;
typedef vector<TCodeword> TMatrix;

std::istream &operator>>(std::istream &in, TCodeword &res) {
    std::string s;
    in >> s;
    res.resize((int) s.length());
    for (int i = 0; i < (int) s.length(); ++i) {
        res[i] = (bool) (s[i] - '0');
    }
    return in;
}

std::ostream &operator<<(std::ostream &out, const TCodeword &c) {
    for (auto i: c) {
        out << i;
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

TCodeword operator^(const TCodeword &a, const TCodeword &b) {
    assert(a.size() == b.size()), "Incorrect sizes in ^";
    TCodeword c(a.size());
    for (int i = 0; i < (int) a.size(); i++)
        c[i] = a[i] ^ b[i];
    return c;
}

TCodeword operator&(const TCodeword &a, const TCodeword &b) {
    assert(a.size() == b.size()), "Incorrect sizes in &";
    TCodeword c(a.size());
    for (int i = 0; i < (int) a.size(); i++)
        c[i] = a[i] & b[i];
    return c;
}

TMatrix operator*(const TMatrix &a, const TMatrix &b) {
    assert(a[0].size() == b.size()), "Incorrect shapes for matmul";
    TMatrix c(a.size());
    for (int i = 0; i < (int) a.size(); i++)
        c[i].resize((int)b[0].size(), false);
    for (int k = 0; k < (int) b.size(); k++)
        for (int i = 0; i < (int) a.size(); i++)
            for (int j = 0; j < (int)b[0].size(); j++)
                c[i][j] = (c[i][j] ^ (a[i][k] & b[k][j]));
    return c;
}

TCodeword operator*(const TMatrix &H, const TCodeword &v) {
    TMatrix c((int)v.size());
    for (int i = 0; i < (int)v.size(); i++)
        c[i].push_back(v[i]);
    TMatrix p = H * c;
    TCodeword res(p.size());
    for (int i = 0; i < (int)res.size(); i++)
        res[i] = p[i][0];
    return res;
}

TCodeword operator*(const TCodeword &v, const TMatrix &H) {
    TMatrix c{v};
    TMatrix p = c * H;
    return p[0];
}

bool IsCodeword(const std::vector<TCodeword> &H, const TCodeword &c) {
    for (bool x: H * c)
        if (x)
            return false;
    return true;
}

std::vector<TCodeword> GetOrtogonal(std::vector<TCodeword> H) {
    std::vector<int> pos(H.size());
    TCodeword is_main;
    int n = H.size(), m = H[0].size();
    for (int i = 0; i < n; ++i) {
        pos[i] = -1;
        for (int j = 0; j < m; ++j) {
            if (H[i][j]) {
                pos[i] = j;
                break;
            }
        }
        assert(pos[i] != -1);
        for (int k = 0; k < n; ++k) {
            if (k != i && H[k][pos[i]]) {
                H[k] = (H[k] ^ H[i]);
            }
        }
        is_main[pos[i]] = true;
    }
    std::vector<TCodeword> res(m - n);
    int idx = 0;
    for (int j = 0; j < m; ++j) {
        if (!is_main[j]) {
            res[idx][j] = true;
            for (int i = 0; i < n; ++i) {
                if (H[i][j]) {
                    res[idx][pos[i]] = true;
                }
            }
            ++idx;
        }
    }
    return res;
}

#endif //ACG_ALP_LDPC_CODEWORD_H
