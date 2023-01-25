#include <utility>
#include <memory>

#include "experiment.h"
#include "utils/parse_data.h"
#include "utils/codeword.h"
#include "algo/algo.h"
#include "algo/qp_admm.h"

using namespace std;

const int THREADS_NUM = 200;
double SNR = -3.0;
shared_ptr<QPADMMDecoder> decoder = make_shared<QPADMMDecoder>(1.95, 0.5, 1000, 1e-5);

double FER(const TMatrix &H, int tests_num = 1000) {
    TMatrix G = GetOrtogonal(H);
    mt19937 rnd(239);
    vector <TCodeword> codewords = gen_random_codewords(G, tests_num, rnd);
    ExperimentResult res = multithread_experiment(decoder, codewords, H, SNR, THREADS_NUM);
    return res.FER();
}

struct PermutationsMatrix {
public:
    PermutationsMatrix(int block_size, const TMatrix &blocks, const vector <vector<int>> diagonals) :
            _block_size(block_size), _blocks(blocks), _diagonals(diagonals) {}

    PermutationsMatrix(int block_size, const TMatrix &H) : _block_size(block_size) {
        assert((int) H.size() % block_size == 0 && (int) H[0].size() % block_size == 0);
        for (int i = 0; i < (int) H.size(); i += block_size) {
            _diagonals.push_back({});
            _blocks.push_back({});
            for (int j = 0; j < (int) H[0].size(); j += block_size) {
                int s = -block_size;
                for (int k = 0; k < block_size; k++)
                    for (int l = 0; l < block_size; l++)
                        if (H[i + k][j + l]) {
                            int ns = (l - k + block_size) % block_size;
                            assert(s == -block_size || s == ns);
                            s = ns;
                        }
                _diagonals.back().push_back(s);
                _blocks.back().push_back(s != -block_size);
            }
        }
        for (int i = 0; i < (int)H.size(); i++) {
            if (to_tmatrix()[i] != H[i]) {
                cout << i << endl;
                cout << to_tmatrix()[i] << endl;
                cout << H[i] << endl;
                cout << endl;
            }
        }
        assert(to_tmatrix() == H);
    }

    TMatrix to_tmatrix() const {
        TMatrix H(_block_size * _blocks.size());
        for (int i = 0; i < (int) H.size(); i++)
            H[i].resize(_block_size * _blocks[0].size(), 0);
        for (int i = 0; i < (int) _blocks.size(); i++)
            for (int j = 0; j < (int) _blocks[i].size(); j++)
                if (_blocks[i][j]) {
                    int s = _diagonals[i][j];
                    assert(0 <= s && s < _block_size);
                    for (int k = 0; k < _block_size; k++) {
                        int l = (s + k + _block_size) % _block_size;
                        H[i * _block_size + k][j * _block_size + l] = 1;
                    }
                }
        return H;
    }

    template<typename Gen>
    PermutationsMatrix random_permute(Gen &rnd) const {
        int i = rnd() % (int)_blocks.size();
        int j = rnd() % (int)_blocks[0].size();
        TMatrix blocks = _blocks;
        vector <vector<int>> diagonals = _diagonals;
        if (!blocks[i][j] or rnd() % 2 == 0)
            blocks[i][j] = !blocks[i][j];
        diagonals[i][j] = rnd() % _block_size;
        return PermutationsMatrix(_block_size, blocks, diagonals);
    }

private:
    int _block_size;
    TMatrix _blocks;
    vector <vector<int>> _diagonals;
};

template<typename Gen>
PermutationsMatrix optimize(PermutationsMatrix H, Gen &rnd, int iters, const string &save_filepath) {
    double error = FER(H.to_tmatrix());
    cout << "initial FER=" << error << endl;
    for (int i = 0; i < iters; i++) {
        PermutationsMatrix newH = H.random_permute(rnd);
        double new_error = FER(newH.to_tmatrix());
        cout << "\tproposal: FER=" << new_error << endl;
        if (new_error < error) {
            H = newH;
            error = new_error;
            cout << "accept, FER=" << error << endl;
            save_matrix(H.to_tmatrix(), save_filepath);
        }
    }
    return H;
}

int main() {
    std::ios::sync_with_stdio(0);
    cout.precision(5);
    cout << fixed;

    PermutationsMatrix H0(20, read_pcm("data/H05.txt"));

    mt19937 rnd(239);
    TMatrix H = optimize(H0, rnd, 10000, "data/optimalH.txt").to_tmatrix();

    cout << FER(H, 10000) << endl;

    return 0;
}
