#ifndef ACG_ALP_LDPC_BP_H
#define ACG_ALP_LDPC_BP_H

#include <utility>

#include "algo/algo.h"
#include "utils/codeword.h"

typedef pair <double, double> Message;

class Node {
public:
    static int counter;

    Node() : _uuid(counter++) {};

    void register_neighbor(int neighbor_uuid) { _neighbors.push_back(neighbor_uuid); }

    void receive_message(int n, Message message) { _received_messages[n] = message; }

    virtual Message message(const Node &n) = 0;

    int uuid() const { return _uuid; }

    vector<int> neighbors() const { return _neighbors; }

    vector<int> _neighbors;
    int _uuid;
    unordered_map<int, Message> _received_messages;
};

int Node::counter = 0;

double phi(double x) { return -log(tanh(x / 2)); }

class CNode : public Node {
public:
    CNode(string name) : _name(name) {
        _uuid = counter++;
    }

    void init() {
        for (int n: _neighbors)
            _received_messages[n] = {0.0, 1.0};
    }

    string name() const { return _name; }

    Message message(const Node &n) override {
        double sum = 0, sgn = 1;
        for (pair<int, Message> p: _received_messages)
            if (p.first != n.uuid()) {
                sum += p.second.first;
                sgn *= p.second.second;
            }
        return {sgn * phi(sum), 1};
    }

private:
    string _name;
};

class VNode : public Node {
public:
    VNode(double snr, double channel_symbol, string name) :
            _channel_llr(llr(channel_symbol, snr)), _name(name) {
        _uuid = counter++;
    }

    void init() {
        for (int n: _neighbors)
            _received_messages[n] = {0.0, 1.0};
    }

    string name() const { return _name; }

    Message message(const Node &n) override {
        double sum = 0;
        for (pair<int, Message> p: _received_messages)
            if (p.first != n.uuid())
                sum += p.second.first;
        return {phi(fabs(_channel_llr + sum)), (_channel_llr + sum <= 0) ? -1 : 1};
    }

    double estimate() const {
        double sum = 0;
        for (pair<int, Message> p: _received_messages)
            sum += p.second.first;
        return _channel_llr + sum;
    }

private:
    double _channel_llr;
    string _name;
};

class TannerGraph {
public:
    int add_v_node(double snr, double channel_value, int i) {
        VNode v(snr, channel_value, "V" + to_string(i));
        _uuid_to_v_nodes[v.uuid()] = (int) _v_nodes.size();
        _v_nodes.push_back(v);
        return v.uuid();
    }

    int add_c_node(int i) {
        CNode c("C" + to_string(i));
        _uuid_to_c_nodes[c.uuid()] = (int) _c_nodes.size();
        _c_nodes.push_back(c);
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
    vector <VNode> _v_nodes;
    vector <CNode> _c_nodes;
    unordered_map<int, int> _uuid_to_v_nodes, _uuid_to_c_nodes;
};

TannerGraph from_biadjacency_matrix(TMatrix h, double snr, const TFVector &channel_word) {
    TannerGraph g;
    int m = (int) h.size(), n = (int) h[0].size();
    vector<int> v_uuids, c_uuids;
    for (int i = 0; i < n; i++)
        v_uuids.push_back(g.add_v_node(snr, channel_word[i], i));
    for (int i = 0; i < m; i++)
        c_uuids.push_back(g.add_c_node(i));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            if (h[i][j])
                g.add_edge(v_uuids[j], c_uuids[i]);
    for (int i = 0; i < n; i++)
        g.v_node(i).init();
    for (int j = 0; j < m; j++)
        g.c_node(j).init();
    return g;
}

class BeliefPropagation {
public:
    BeliefPropagation(TannerGraph g, TMatrix h, int max_iter) :
            _h(std::move(h)), _graph(std::move(g)), _max_iter(max_iter) {}

    void c_receive_messages() {
        for (int i = 0; i < _graph.c_size(); i++) {
            CNode &c = _graph.c_node(i);
            for (int n: c.neighbors()) {
                VNode &v = _graph.v_node_by_uuid(n);
//                cout << "\t" << v.name() << " -> " << c.name() << " : " << v.message(c) << endl;
                c.receive_message(n, v.message(c));
            }
        }
    }

    void v_receive_messages() {
        for (int i = 0; i < _graph.v_size(); i++) {
            VNode &v = _graph.v_node(i);
            for (int n: v.neighbors()) {
                CNode &c = _graph.c_node_by_uuid(n);
//                cout << "\t" << c.name() << " -> " << v.name() << " : " << c.message(v) << endl;
                auto msg = c.message(v);
                v.receive_message(n, msg);
            }
        }
    }

    pair<TCodeword, bool> decode() {
        c_receive_messages();
//        cout << "----------------------------------------" << endl;
        for (int _ = 0; _ < _max_iter; _++) {
            v_receive_messages();
            c_receive_messages();
//            cout << "----------------------------------------" << endl;

            TCodeword estimation;
            for (int i = 0; i < _graph.v_size(); i++)
                estimation.push_back(_graph.v_node(i).estimate() <= 0);

            if (IsCodeword(_h, estimation))
                return {estimation, true};
        }
        return {TCodeword(), false};
    }

private:
    TannerGraph _graph;
    TMatrix _h;
    int _max_iter;
};


class BeliefPropagationDecoder : public Decoder {
public:
    explicit BeliefPropagationDecoder(int max_iter) : _max_iter(max_iter) {}

    pair<TCodeword, bool> decode(const TMatrix &H, const TFVector &channel_word, double snr) override {
        TannerGraph graph = from_biadjacency_matrix(H, snr, channel_word);
        BeliefPropagation bp(graph, H, _max_iter);
        return bp.decode();
    }

    string name() const override { return "BP"; }

private:
    int _max_iter;
};

#endif //ACG_ALP_LDPC_BP_H