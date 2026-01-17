#ifndef __WERNER_ALGOUB_H
#define __WERNER_ALGOUB_H

#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"

#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cassert>
#include <cmath>

using namespace std;

class WernerAlgo_UB : public AlgorithmBase {
public:
    // 為了保持精度，特定計算使用 long double，但介面仍維持 double
    #define ld long double
    
    WernerAlgo_UB(Graph graph,
               vector<pair<int,int>> requests,
               map<SDpair, vector<Path>> paths);

    void run();

private:
    // ===== Werner DP Label =====
    enum class Op : unsigned char { LEAF = 0, CONT = 1, MERGE = 2 };
    class ZLabel {
    public:
        ld B = 1000.0;
        ld Z = 1000.0;   // 成本 / 目標
        int a = -1, b = -1, t = -1, k = -1; 
        int left_id = -1, right_id = -1, parent_id = -1;
        Op op = Op::LEAF;
        vector<int> ent_time; // 用於儲存 entanglement 產生的時間區段
        
        ZLabel(){}
        ZLabel(ld _B, ld _Z, Op _op, int _a, int _b, int _t, int _k, int pid = -1, int lid = -1, int rid = -1)
        : B(_B), Z(_Z), op(_op), a(_a), b(_b), t(_t), k(_k), parent_id(pid), left_id(lid), right_id(rid) {}
    };

    struct DPParam{
        ld eps_bucket, Zhat, Zmin, eta, T;
        int tau_max;
    } dpp;

    // ===== 參數 / 對偶變數 =====
    double epsilon = 0.35;
    double obj = 0.0;
    
    // alpha[i]: 第 i 個 request 的 dual
    vector<double> alpha;                 
    // beta[v][t]: 節點 v 在時間 t 的 dual
    vector<vector<double>> beta;          
    
    // Primal 變數 x[request_id][shape] = flow amount
    vector<map<Shape_vector, double>> x;

    // ===== 全時間 DP 表：DP_table[time][a][b] =====
    // a, b 是 path 上的 index
    vector<vector<vector<vector<ZLabel>>>> DP_table;

    // ===== 內部函式 =====
    void variable_initialize();
    Shape_vector separation_oracle();

    // DP 核心
    void run_dp_in_t(const Path& path, const DPParam& dpp, int t);

    // 輔助函式
    void pareto_prune_byZ(vector<ZLabel>& cand);
    void bucket_by_Z(vector<ZLabel>& cand);
    Shape_vector backtrack_shape(ZLabel leaf, const vector<int>& path);
    
    // 評估函式
    int split_dis(int s, int d, const ZLabel& L);
    pair<double, ZLabel> eval_best_J(int s, int d, int t, double alp);
};

#endif // __WERNER_ALGOUB_H