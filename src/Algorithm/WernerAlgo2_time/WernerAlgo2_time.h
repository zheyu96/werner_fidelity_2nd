#ifndef __WERNER_ALGO2_TIME_H
#define __WERNER_ALGO2_TIME_H

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

using namespace std;

class WernerAlgo2_time : public AlgorithmBase {
public:
    // 使用專用別名代替全域 #define，避免污染標準庫
    using Float = long double;

    WernerAlgo2_time(Graph graph,
               vector<pair<int,int>> requests,
               map<SDpair, vector<Path>> paths, double epsilon);

    void run() override;

private:
    enum class Op : unsigned char { LEAF = 0, CONT = 1, MERGE = 2 };
    
    class ZLabel {
    public:
        Float B = 1000.0;
        Float Z = 1000.0;   
        int a = -1, b = -1, t = -1, k = -1; 
        int left_id = -1, right_id = -1, parent_id = -1;
        Op op = Op::LEAF;
        vector<int> ent_time;
        
        ZLabel(){}
        ZLabel(Float _B, Float _Z, Op _op, int _a, int _b, int _t, int _k, int pid = -1, int lid = -1, int rid = -1)
        : B(_B), Z(_Z), op(_op), a(_a), b(_b), t(_t), k(_k), parent_id(pid), left_id(lid), right_id(rid) {}
    };

    struct DPParam {
        Float eps_bucket, Zhat, Zmin, eta, T;
        int tau_max;
    } dpp;

    double epsilon = 0.35;
    Float obj = 0.0;
    vector<Float> alpha;                 
    vector<vector<Float>> beta;          
    vector<map<Shape_vector, double>> x;
    vector<vector<vector<vector<ZLabel>>>> DP_table;

    // 統計數據
    long long total_labels_generated = 0;
    long long total_labels_after_bucket = 0;

    void variable_initialize();
    Shape_vector separation_oracle();

    pair<long long, long long> run_dp_with_stats(const Path& path, const DPParam& dpp, int t);
    void run_dp_in_t(const Path& path, const DPParam& dpp, int t);

    void pareto_prune_byZ(vector<ZLabel>& cand);
    void bucket_by_Z(vector<ZLabel>& cand);

    Shape_vector backtrack_shape(ZLabel leaf, const vector<int>& path);
    long long split_dis(int s, int d, WernerAlgo2_time::ZLabel& L); // 修改為 long long
    pair<Float, WernerAlgo2_time::ZLabel> eval_best_J(int s, int d, int t, Float alp);
};

#endif