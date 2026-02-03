#include "WernerAlgo2_time.h"
#include <fstream> 
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

WernerAlgo2_time::WernerAlgo2_time(Graph graph, vector<pair<int,int>> requests, map<SDpair, vector<Path>> paths, double eps)
    : AlgorithmBase(graph, requests, paths) {
    algorithm_name = "ZFA2_time";
    epsilon = eps;
}

void WernerAlgo2_time::variable_initialize() {
    int m = (int)requests.size() + graph.get_num_nodes() * graph.get_time_limit();
    Float delta = (1 + (Float)epsilon) * (1.0L / powl((1 + (Float)epsilon) * m, 1.0L / (Float)epsilon));
    obj = m * delta;
    alpha.assign(requests.size(), delta);
    x.clear();
    x.resize(requests.size());
    
    int V = graph.get_num_nodes();
    int T = graph.get_time_limit();
    dpp.eps_bucket = graph.get_bucket_eps();
    Float F_th = graph.get_fidelity_threshold();
    Float w_th = (4.0L * F_th - 1.0L) / 3.0L;
    dpp.Zhat = sqrtl(-logl(w_th)) + 1e-9L;
    dpp.Zmin = graph.get_Zmin();
    dpp.T = time_limit - 1;
    dpp.tau_max = min(time_limit - 1, 5);
    dpp.eta = (Float)graph.get_tao() / (Float)graph.get_time_limit();
    
    beta.assign(V, vector<Float>(T, (Float)INF));
    for (int v = 0; v < V; ++v) {
        for (int t = 0; t < T; ++t) {
            int cap = graph.get_node_memory_at(v, t);
            beta[v][t] = (cap == 0) ? (Float)INF : (delta / (Float)cap);
        }
    }
}

Shape_vector WernerAlgo2_time::separation_oracle() {
    Float most_violate = 1e18L;
    Shape_vector todo_shape;
    total_labels_generated = 0;
    total_labels_after_bucket = 0;

    for (size_t i = 0; i < requests.size(); i++) { // 修復 Sign-compare
        int src = requests[i].first, dst = requests[i].second;
        vector<Path> paths = get_paths(src, dst);
        for (size_t p = 0; p < paths.size(); p++) {
            int T_size = (int)dpp.T + 5;
            int n_size = (int)paths[p].size() + 5;
            DP_table.assign(T_size, vector<vector<vector<ZLabel>>>(n_size, vector<vector<ZLabel>>(n_size)));

            for (int t = 1; t <= dpp.T; t++) {
                pair<long long, long long> s = run_dp_with_stats(paths[p], dpp, t);
                total_labels_generated += s.first;
                total_labels_after_bucket += s.second;

                auto cur_val = eval_best_J(0, (int)paths[p].size() - 1, t, alpha[i]);
                if (cur_val.first / (Float)graph.path_Pr(paths[p]) < most_violate) {
                    most_violate = cur_val.first / (Float)graph.path_Pr(paths[p]);
                    todo_shape = backtrack_shape(cur_val.second, paths[p]);
                }
            }
        }
    }

    if (total_labels_generated > 0) {
        cerr << fixed << setprecision(2) << "[DP Label Stats] Total Gen: " << total_labels_generated 
             << " | After Bucket: " << total_labels_after_bucket 
             << " | Compression: " << (1.0 - (double)total_labels_after_bucket / total_labels_generated) * 100.0 << "%" << endl;
    }
    return todo_shape;
}

pair<long long, long long> WernerAlgo2_time::run_dp_with_stats(const Path& path, const DPParam& dpp, int t) {
    long long gen_t = 0;
    long long kept_t = 0;
    const int n = (int)path.size();

    for (int a = 0; a < n - 1; a++) {
        for (int b = a + 1; b < n; b++) {
            int s = path[a], e = path[b];
            vector<ZLabel> cand;

            if (a + 1 == b) {
                for (int tlen = 1; tlen <= dpp.tau_max; tlen++) {
                    int st = t - tlen;
                    if (st < 1) continue;
                    Float Bleaf = 0.0;
                    for (int tt = st; tt <= t; tt++) Bleaf += beta[s][tt] + beta[e][tt];
                    Float Zcur = 1.0L - (1.0L - (Float)graph.get_link_werner(s, e)) / (Float)tlen;
                    if (Zcur <= 0) continue;
                    Float Zleaf = sqrtl(-logl(Zcur));
                    if (Zleaf <= dpp.Zhat) {
                        ZLabel L(Bleaf, Zleaf, Op::LEAF, a, b, t, -1);
                        L.ent_time = {st, t};
                        cand.push_back(L);
                    }
                }
            }

            if (t > 1) {
                auto& pre = DP_table[t - 1][a][b];
                for (size_t p_id = 0; p_id < pre.size(); p_id++) {
                    Float Zp = pre[p_id].Z + dpp.eta;
                    if (Zp <= dpp.Zhat) {
                        Float Bp = pre[p_id].B + beta[s][t] + beta[e][t];
                        cand.push_back(ZLabel(Bp, Zp, Op::CONT, a, b, t, -1, (int)p_id));
                    }
                }
            }

            if (t > 1) {
                for (int k = a + 1; k < b; k++) {
                    auto &L1 = DP_table[t - 1][a][k], &L2 = DP_table[t - 1][k][b];
                    for (size_t lid = 0; lid < L1.size(); lid++) {
                        for (size_t rid = 0; rid < L2.size(); rid++) {
                            Float Zp = sqrtl(powl(L1[lid].Z + dpp.eta, 2) + powl(L2[rid].Z + dpp.eta, 2));
                            if (Zp <= dpp.Zhat) {
                                Float Bp = L1[lid].B + L2[rid].B + beta[s][t] + beta[e][t];
                                cand.push_back(ZLabel(Bp, Zp, Op::MERGE, a, b, t, k, -1, (int)lid, (int)rid));
                            }
                        }
                    }
                }
            }

            gen_t += (long long)cand.size();
            bucket_by_Z(cand);
            kept_t += (long long)cand.size();
            DP_table[t][a][b] = cand;
        }
    }
    return {gen_t, kept_t};
}

void WernerAlgo2_time::run_dp_in_t(const Path& path, const DPParam& dpp, int t) {
    run_dp_with_stats(path, dpp, t);
}

void WernerAlgo2_time::pareto_prune_byZ(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    sort(cand.begin(), cand.end(), [](const ZLabel& x, const ZLabel& y) {
        if (x.Z != y.Z) return x.Z < y.Z;
        return x.B < y.B;
    });
    vector<ZLabel> kept;
    Float bestB = (Float)INF;
    for (auto& L : cand) {
        if (L.B + 1e-15L < bestB) {
            kept.push_back(L);
            bestB = L.B;
        }
    }
    cand.swap(kept);
}

void WernerAlgo2_time::bucket_by_Z(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    Float q = 1.0L + dpp.eps_bucket;
    Float invLogQ = 1.0L / logl(q);
    map<long long, ZLabel> buckets; // 使用 long long 為索引避開 double 問題
    for (auto& L : cand) {
        long long k = (L.Z <= dpp.Zmin) ? 0 : (long long)floorl(logl(L.Z / dpp.Zmin) * invLogQ + 1e-12L);
        if (k < 0) k = 0;
        if (buckets.find(k) == buckets.end() || L.B + 1e-15L < buckets[k].B)
            buckets[k] = L;
    }
    vector<ZLabel> bucketed;
    for (auto& it : buckets) bucketed.push_back(it.second);
    pareto_prune_byZ(bucketed);
    sort(bucketed.begin(), bucketed.end(), [](const ZLabel& x, const ZLabel& y) { return x.Z < y.Z; });
    cand.swap(bucketed);
}

Shape_vector WernerAlgo2_time::backtrack_shape(ZLabel leaf, const vector<int>& path) {
    int left_id = path[leaf.a], right_id = path[leaf.b];
    if (leaf.op == Op::LEAF) {
        Shape_vector result;
        result.push_back({left_id, {{leaf.ent_time[0], leaf.ent_time[1]}}});
        result.push_back({right_id, {{leaf.ent_time[0], leaf.ent_time[1]}}});
        return result;
    }
    if (leaf.op == Op::CONT) {
        ZLabel pre_label = DP_table[leaf.t - 1][leaf.a][leaf.b][leaf.parent_id];
        Shape_vector last = backtrack_shape(pre_label, path);
        last.front().second[0].second++;
        last.back().second[0].second++;
        return last;
    }
    if (leaf.op == Op::MERGE) {
        ZLabel left_leaf = DP_table[leaf.t - 1][leaf.a][leaf.k][leaf.left_id];
        ZLabel right_leaf = DP_table[leaf.t - 1][leaf.k][leaf.b][leaf.right_id];
        Shape_vector left_res = backtrack_shape(left_leaf, path);
        Shape_vector right_res = backtrack_shape(right_leaf, path);
        Shape_vector result = left_res;
        result.back().second.push_back(right_res.front().second.front());
        for (size_t i = 1; i < right_res.size(); i++) result.push_back(right_res[i]);
        result.front().second[0].second++;
        result.back().second[0].second++;
        return result;
    }
    return {};
}

long long WernerAlgo2_time::split_dis(int s, int d, WernerAlgo2_time::ZLabel& L) {
    if (L.op != WernerAlgo2_time::Op::MERGE || L.k < 0) return 200000000000000000LL; // 修復溢位
    return (long long)abs(((s + d) / 2) - L.k);
}

pair<WernerAlgo2_time::Float, WernerAlgo2_time::ZLabel> WernerAlgo2_time::eval_best_J(int s, int d, int t, Float alp) {
    Float bestJ = 1e18L;
    long long bestdis = 200000000000000000LL; // 修復溢位
    bool flag = false;
    ZLabel tmp;
    for (auto& L : DP_table[t][s][d]) {
        Float J = (alp + L.B) * expl(L.Z * L.Z);
        long long dis = split_dis(s, d, L);
        if (J + 1e-15L < bestJ || (fabsl(J - bestJ) <= 1e-15L && dis < bestdis)) {
            bestJ = J; tmp = L; bestdis = dis; flag = true;
        }
    }
    return flag ? make_pair(bestJ, tmp) : make_pair((Float)INF, tmp);
}

void WernerAlgo2_time::run() {
    int round = 1;
    while (round-- && !requests.empty()) {
        variable_initialize();
        Float eps_v = 1e-4L;
        while (obj + eps_v < 1.0L) {
            Shape_vector shape = separation_oracle();
            if (shape.empty()) break;
            Float q = 1.0L;
            for (size_t i = 0; i < shape.size(); i++) {
                map<int, int> need;
                for (auto& range : shape[i].second)
                    for (int t = range.first; t <= range.second; t++) need[t]++;
                for (auto& p : need) q = min(q, (Float)graph.get_node_memory_at(shape[i].first, p.first) / (Float)p.second);
            }
            if (q <= 1e-10L) break;
            int req_idx = -1;
            for (size_t i = 0; i < requests.size(); i++) {
                if (requests[i] == make_pair(shape.front().first, shape.back().first))
                    if (req_idx == -1 || alpha[req_idx] > alpha[i]) req_idx = (int)i;
            }
            if (req_idx == -1) break;
            x[req_idx][shape] += (double)q;
            Float ori = alpha[req_idx];
            alpha[req_idx] *= (1.0L + (Float)epsilon * q);
            obj += (alpha[req_idx] - ori);
            for (size_t i = 0; i < shape.size(); i++) {
                map<int, int> need;
                for (auto& range : shape[i].second)
                    for (int t = range.first; t <= range.second; t++) need[t]++;
                for (auto& p : need) {
                    int node_id = shape[i].first, t = p.first;
                    Float original = beta[node_id][t];
                    beta[node_id][t] *= (1.0L + (Float)epsilon * (q / ((Float)graph.get_node_memory_at(node_id, t) / (Float)p.second)));
                    obj += (beta[node_id][t] - original) * (Float)p.second;
                }
            }
        }
        
        vector<pair<double, Shape_vector>> sorted_shapes;
        for(size_t i = 0; i < requests.size(); i++)
            for(auto& P : x[i]) sorted_shapes.push_back({P.second, P.first});
        sort(sorted_shapes.begin(), sorted_shapes.end(), [](auto& a, auto& b){ return a.first > b.first; });

        vector<bool> used(requests.size(), false);
        vector<int> finished;
        for(auto& P : sorted_shapes) {
            Shape shape_obj = Shape(P.second);
            int idx = -1;
            for(size_t i=0; i<requests.size(); i++) {
                if(!used[i] && requests[i] == make_pair(shape_obj.get_node_mem_range().front().first, shape_obj.get_node_mem_range().back().first)) {
                    idx = (int)i; break;
                }
            }
            if(idx != -1 && graph.check_resource(shape_obj, true, true)) {
                used[idx] = true; graph.reserve_shape(shape_obj, true); finished.push_back(idx);
            }
        }
        sort(finished.rbegin(), finished.rend());
        for(int fin : finished) requests.erase(requests.begin() + fin);
    }
    update_res();
    cerr << "[" << algorithm_name << "] end" << endl;
}