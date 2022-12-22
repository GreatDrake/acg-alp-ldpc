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

const float eps = 1e-8;
const float SNR = 1.0;
const int MAX_ITER = 100;
const int THREADS_NUM = 26;
const int LOG_FREQ = 100;

typedef vector<vector<int> > TMatrix;
typedef vector<int> TVector;
typedef vector<float> TFVector;

float llr_mean(float snr) { return 4 * 0.5 * pow(10, (snr / 10)); }

float llr_variance(float snr) { return 8 * 0.5 * pow(10, (snr / 10)); }

TVector operator*(const TMatrix &H, const TVector &v) {
    assert(H[0].size() == v.size()), "Incorrect shapes for matmul";
    TVector w(H.size(), 0);
    for (int i = 0; i < (int) H.size(); i++)
        for (int j = 0; j < (int) H[i].size(); j++)
            w[i] = (w[i] + H[i][j] * v[j]) % 2;
    return w;
}

TMatrix read_pcm(const string &filename) {
    ifstream fin(filename);
    int n, m;
    fin >> n >> m;
    int dc, dr;
    fin >> dc >> dr;
    vector<int> d_col(n);
    for (int i = 0; i < n; ++i)
        fin >> d_col[i];
    vector<int> d_row(m);
    for (int i = 0; i < m; ++i)
        fin >> d_row[i];
    TMatrix result(m, vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        std::vector<int> id(dc);
        for (int j = 0; j < dc; ++j) {
            fin >> id[j];
            --id[j];
        }
        for (int j = 0; j < dc; ++j)
            if (id[j] >= 0)
                result[id[j]][i] = 1;
    }
    return result;
}

vector<TVector> read_codewords(const string &filename) {
    ifstream fin(filename);
    int n;
    fin >> n;
    vector<TVector> codewords;
    for (int i = 0; i < n; i++) {
        string s;
        fin >> s;
        TVector codeword;
        for (char c: s)
            codeword.push_back(c - '0');
        codewords.push_back(codeword);
    }
    return codewords;
}

class Node {
public:
    static int counter;

    Node() : _uuid(counter++) {};

    void register_neighbor(int neighbor_uuid) { _neighbors.push_back(neighbor_uuid); }

    void receive_message(int n, float message) { _received_messages[n] = message; }

    virtual float message(const Node &n) = 0;

    int uuid() const { return _uuid; }

    vector<int> neighbors() const { return _neighbors; }

    vector<int> _neighbors;
    int _uuid;
    unordered_map<int, float> _received_messages;
};

int Node::counter = 0;

float phi(float x) { return -log(tanh(x / 2)); }

class CNode : public Node {
public:
    explicit CNode() {
        _uuid = counter++;
    }

    void init() {
        for (int n: _neighbors)
            _received_messages[n] = 0;
    }

    float message(const Node &n) override {
        float sum = 0, sgn = 1;
//        cout << "\t message from C:" << uuid() << " to V:" << n.uuid() << " : (";
        for (pair<int, float> p: _received_messages)
            if (p.first != n.uuid()) {
                sum += phi(abs(p.second));
//                cout << "phi(|" << p.second << "|)(" << phi(abs(p.second)) << ") +";
                if (p.second < -eps)
                    sgn *= -1;
                if (abs(p.second) <= eps)
                    sgn = 0;
            }
//        cout << ") * " << sgn << " : " << ((abs(sum) < eps) ? 0 : sgn * phi(sum)) << endl;
        if (abs(sum) < eps)
            return 0;
        return sgn * phi(sum);
    }
};

float log_density(float v, float snr) {
    float mean = llr_mean(snr);
    float var = llr_variance(snr);
    return -(v - mean) * (v - mean) / (2 * var);
}

class VNode : public Node {
public:
    explicit VNode(float snr, int channel_symbol) : _channel_symbol(channel_symbol) {
        _uuid = counter++;
        _channel_symbol = channel_symbol;
//        _channel_llr = log((1 - BSC_P) / BSC_P);
//        if (channel_symbol)
//            _channel_llr = -_channel_llr;
        _channel_llr = log_density(-channel_symbol, snr) - log_density(channel_symbol, snr);
//        cout << "VNode " << uuid() << " : " << _channel_llr << endl;
    }

    void init() {
        for (int n: _neighbors)
            _received_messages[n] = 0;
    }

    float message(const Node &n) override {
        float sum = 0;
//        cout << "\t message from V:" << uuid() << " to C:" << n.uuid() << " : ";
        for (pair<int, float> p: _received_messages)
            if (p.first != n.uuid()) {
                sum += p.second;
//                cout << p.second << " + ";
            }
//        cout << "llr(" << _channel_llr << ")" << endl;
        return _channel_llr + sum;
    }

    float estimate() const {
        float sum = 0;
        for (pair<int, float> p: _received_messages)
            sum += p.second;
//        cout << "estimation in " << uuid() << " : " << _channel_llr << " + " << sum << endl;
        return _channel_llr + sum;
    }

private:
    int _channel_symbol;
    float _channel_llr;
};

class TannerGraph {
public:
    int add_v_node(float snr, float channel_value) {
        VNode v(snr, channel_value);
        _uuid_to_v_nodes[v.uuid()] = (int) _v_nodes.size();
        _v_nodes.push_back(v);
        _v_nodes_to_uuid[_uuid_to_v_nodes[v.uuid()]] = v.uuid();
        return v.uuid();
    }

    int add_c_node() {
        CNode c;
        _uuid_to_c_nodes[c.uuid()] = (int) _c_nodes.size();
        _c_nodes.push_back(c);
        _c_nodes_to_uuid[_uuid_to_c_nodes[c.uuid()]] = c.uuid();
        return c.uuid();
    }

    size_t v_size() const { return _v_nodes.size(); }

    size_t c_size() const { return _c_nodes.size(); }

    void add_edge(int vnode_uuid, int cnode_uuid) {
        _v_nodes[_uuid_to_v_nodes[vnode_uuid]].register_neighbor(cnode_uuid);
        _c_nodes[_uuid_to_c_nodes[cnode_uuid]].register_neighbor(vnode_uuid);
    }

    CNode &c_node(int i) { return _c_nodes[i]; }

    VNode &v_node(int i) { return _v_nodes[i]; }

    CNode &c_node_by_uuid(int uuid) { return c_node(_uuid_to_c_nodes[uuid]); }

    VNode &v_node_by_uuid(int uuid) { return v_node(_uuid_to_v_nodes[uuid]); }

protected:
    vector<VNode> _v_nodes;
    vector<CNode> _c_nodes;
    unordered_map<int, int> _uuid_to_v_nodes, _uuid_to_c_nodes;
    unordered_map<int, int> _v_nodes_to_uuid, _c_nodes_to_uuid;
};

TannerGraph from_biadjacency_matrix(TMatrix h, float snr, const TFVector &channel_word) {
    TannerGraph g;
    int m = (int) h.size(), n = (int) h[0].size();
    vector<int> v_uuids, c_uuids;
    for (int i = 0; i < n; i++)
        v_uuids.push_back(g.add_v_node(snr, channel_word[i]));
    for (int i = 0; i < m; i++)
        c_uuids.push_back(g.add_c_node());
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            if (abs(h[i][j]) > eps)
                g.add_edge(v_uuids[j], c_uuids[i]);
    for (int i = 0; i < n; i++)
        g.v_node(i).init();
    for (int j = 0; j < m; j++)
        g.c_node(j).init();
    return g;
}

class BeliefPropagation {
public:
    BeliefPropagation(const TannerGraph &g, TMatrix h, int max_iter) :
            _h(std::move(h)), _graph(g), _max_iter(max_iter) {}

    void c_receive_messages() {
        for (int i = 0; i < _graph.c_size(); i++) {
            CNode &c = _graph.c_node(i);
            for (int n: c.neighbors()) {
                VNode &v = _graph.v_node_by_uuid(n);
                c.receive_message(n, v.message(c));
            }
        }
    }

    void v_receive_messages() {
        for (int i = 0; i < _graph.v_size(); i++) {
            VNode &v = _graph.v_node(i);
            for (int n: v.neighbors()) {
                CNode &c = _graph.c_node_by_uuid(n);
                auto msg = c.message(v);
//                cout << "\t" << v.uuid() << " receives " << msg << " from " << c.uuid() << endl;
                v.receive_message(n, msg);
            }
        }
    }

    pair<TVector, bool> decode() {
        c_receive_messages();
        for (int _ = 0; _ < _max_iter; _++) {
            v_receive_messages();
            c_receive_messages();

            TVector estimation;
            for (int i = 0; i < _graph.v_size(); i++)
                estimation.push_back((_graph.v_node(i).estimate() <= 0) ? 1 : 0);

//            cout << "\testimation: ";
//            for (int i = 0; i < _graph.v_size(); i++)
//                cout << estimation[i];
//            cout << endl;

            TVector syndrome = _h * estimation;
            bool fits = true;
            for (float c: syndrome)
                if (abs(c) > eps)
                    fits = false;
            if (fits)
                return {estimation, true};
        }
        return {TVector(), false};
    }

private:
    TannerGraph _graph;
    TMatrix _h;
    int _max_iter;
};

pair<TVector, bool> test(const TMatrix &H, const TFVector &channel_word, float snr, int max_iter) {
    TannerGraph graph = from_biadjacency_matrix(H, snr, channel_word);
    BeliefPropagation bp(graph, H, max_iter);
    return bp.decode();
}

template<typename Gen>
TFVector transmit(double snr, const TVector &c, Gen &rnd) {
    double llrVariance = llr_variance(snr);
    double llrMean = llr_mean(snr);
    double llrSigma = sqrt(llrVariance);
    TFVector res((int) c.size());
    normal_distribution<double> dst(llrMean, llrSigma);
    for (int i = 0; i < (int) c.size(); i++)
        res[i] = (c[i] ? 1.0 : -1.0) * dst(rnd);
    return res;
}

struct ThreadArgs {
    ThreadArgs(const vector<TVector> &codewords, const TMatrix &H) : codewords(codewords), H(H) {}

    const vector<TVector> codewords;
    const TMatrix &H;
    int last_uuid = 0;
    pthread_mutex_t uuid_id_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;
};

void *exp(void *arg) {
    ThreadArgs &args = *(ThreadArgs *) arg;
    int result = 0;
    while (true) {
        assert(!pthread_mutex_lock(&args.uuid_id_mutex));
        TVector codeword;
        if (args.last_uuid < (int) args.codewords.size()) {
            codeword = args.codewords[args.last_uuid];
            args.last_uuid++;
        }
        int uuid = args.last_uuid;
        assert(!pthread_mutex_unlock(&args.uuid_id_mutex));
        if (codeword.empty())
            break;

        std::mt19937 rnd(args.last_uuid);
        TFVector transmitted = transmit(SNR, codeword, rnd);
        pair<TVector, bool> p = test(args.H, transmitted, SNR, MAX_ITER);
        if (uuid % LOG_FREQ == 0) {
            assert(!pthread_mutex_lock(&args.log_mutex));
            cout << "Checking " << uuid << endl;
            assert(!pthread_mutex_unlock(&args.log_mutex));
        }
        result += (p.second ? 1 : 0);
    }
    return (void *) new int(result);
}

int main() {
//    float bsc_p = 0.1;
//    float std = 0.3;
//    int max_iter = 1000;
//
//    TMatrix H{{1, 1, 1, 1, 0, 0, 0, 0, 0, 0},
//              {1, 0, 0, 0, 1, 1, 1, 0, 0, 0},
//              {0, 1, 0, 0, 1, 0, 0, 1, 1, 0},
//              {0, 0, 1, 0, 0, 1, 0, 1, 0, 1},
//              {0, 0, 0, 1, 0, 0, 1, 0, 1, 1}};
//
//    mt19937 rnd(239);
//    TVector v = {1, 1, 0, 0, 1, 0, 0, 0, 0, 0};
//    TFVector tr = transmit(SNR, v, rnd);
//    pair<TVector, bool> p = test(H, tr, SNR, MAX_ITER);
//    cout << p.second << endl;

    TMatrix H = read_pcm("H.txt");
    vector<TVector> codewords = read_codewords("codewords.txt");
    vector<int> uuids;
    ThreadArgs args(codewords, H);
    vector<pthread_t> threads(THREADS_NUM);
    for (int q = 0; q < THREADS_NUM; q++)
        assert(!pthread_create(&threads[q], nullptr, exp, (void *) &args));
    int correct = 0;
    for (int q = 0; q < THREADS_NUM; q++) {
        void *data;
        assert(!pthread_join(threads[q], &data));
        correct += *(int *) data;
        delete (double *) data;
    }

    cout << (float) correct / (float) codewords.size() << endl;

    return 0;
}