#ifndef __MYALGO1_H
#define __MYALGO1_H

#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"
#include <cfloat>

using namespace std;
class state{
public:
    double B,W;
    state(double B=DBL_MAX,double W=DBL_MAX):B(b),W(w){};
    state operator +(const state& other) const{
        return state(B+other.B,W+other.W);
    }
    state& operator +=(const state& other) {
        B+=other.B;
        W+=other.W;
        return *this;
    }
    double get_eval(double lambda){
        return B+lambda*W;
    }
};
class strategy{
public:
    vector<vector<vector<state>>> dp;
    vector<vector<vector<bool>>> caled;
    vector<vector<vector<int>>> par;
    state min_shape_in_path
    strategy(int path_len,int time_len){
        dp.clear();
        dp.resize(path_len);
        par.clear();
        par.resize(path_len);
        caled.clear();
        caled.resize(path_len);
        for(int i = 0; i < path_len; i++) {
            dp[i].resize(path_len);
            par[i].resize(path_len);
            caled[i].resize(path_len);
            for(int j = 0; j < path_len; j++) {
                dp[i][j].resize(time_len,{INF,INF});
                par[i][j].resize(time_len, -2);
                caled[i][j].resize(time_len, false);
            }
        }
    }
};
class MyAlgo7 : public AlgorithmBase {
    vector<double> alpha;
    vector<vector<double>> beta;
    vector<map<Shape_vector, double>> x;
    map<double,strategy> dp_table;
    double epsilon, obj;
    void variable_initialize();
    Shape_vector separation_oracle();
    pair<Shape_vector, double> find_min_shape(int src, int dst, double alp);
    state recursion_calculate_min_shape(int left, int right, int t, vector<int> &path);
    Shape_vector recursion_find_shape(int left, int right, int t, vector<int> &path);
public:
    MyAlgo1(Graph graph, vector<pair<int, int>> requests, map<SDpair, vector<Path>> paths);
    void run();
};

#endif