#ifndef ACG_ALP_LDPC_PARSE_DATA_H
#define ACG_ALP_LDPC_PARSE_DATA_H

#include "codeword.h"

TMatrix read_pcm(const string &filename) {
    ifstream fin(filename);
    TMatrix result;
    string s;
    bool t;
    while (fin >> s) {
        if (s.empty())
            continue;
        result.push_back({});
        if (s.back() != ',')
            s += ',';
        for (char i : s) {
            if (i == ',') {
                result.back().push_back(t);
            } else
                t = (i == '1');
        }
    }
    return result;
}


vector<TCodeword> read_codewords(const string &filename) {
    ifstream fin(filename);
    int n;
    fin >> n;
    vector<TCodeword> codewords;
    for (int i = 0; i < n; i++) {
        string s;
        fin >> s;
        TCodeword codeword;
        for (char c: s)
            codeword.push_back(c == '0');
        codewords.push_back(codeword);
    }
    return codewords;
}

void save_matrix(const TMatrix &H, const string &filepath) {
    ofstream fout(filepath);
    for (auto v : H) {
        for (int i = 0; i < (int)v.size(); i++) {
            fout << v[i];
            if (i != (int)v.size() - 1)
                fout << ',';
        }
        fout << endl;
    }
}

#endif //ACG_ALP_LDPC_PARSE_DATA_H
