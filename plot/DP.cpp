#include <bits/stdc++.h>
using namespace std;
bool DEBUG = false;
const double EPS = 1e-9;
const double INF = 1e9;
const double mT = 40e-3;
const double A = 0.25, B = 0.75, n = 2;
double bar(double F) {
    return (1.0 - F);
}
double Fswap(double Fa, double Fb, bool maximize = true) {
    if(Fa <= A + EPS || Fb <= A + EPS) return (maximize ? 0 : INF);
    if(Fa >= INF / 2 || Fb >= INF / 2) return (maximize ? 0 : INF);
    return Fa * Fb + (1.0 / 3.0) * bar(Fa) * bar(Fb);
}
double t2F(double t, bool maximize = true) {
    if(t >= INF / 1e2) return (maximize ? 0 : INF);
    return A + B * exp(-pow(t / mT, n));
}
double F2t(double F, bool maximize = true) {
    if(F <= A + EPS) return (maximize ? INF : INF);
    return mT * pow(-log((F - A) / B), 1.0 / n);
}
double pass_tao(double F, double tao, bool maximize = true) {
    if(F >= INF / 1e2) return INF;
    return t2F(F2t(F, maximize) + tao, maximize);
}

class DP {
    vector<vector<vector<vector<double>>>> dp;
    vector<vector<vector<vector<bool>>>> caled;
    double tao;
    double solve_fidelity_max(int left, int right, int t, int state);
    double solve_fidelity_min(int left, int right, int t, int state);
public:
    double max_fidelity(int T, double _tao);
    double min_fidelity(int T, double _tao);
};

double DP::max_fidelity(int T, double _tao) {
    if(T == 1) return 0.97;
    tao = _tao;
    double res = 0;
    int cnt = 0;

    for(int l = min(T, 7); l <= min(T, 7); l++) {
        dp.clear();
        dp.resize(l);
        caled.clear();
        caled.resize(l);
        for(int i = 0; i < l; i++) {
            dp[i].resize(l);
            caled[i].resize(l);
            for(int j = 0; j < l; j++) {
                dp[i][j].resize(T);
                caled[i][j].resize(T);
                for(int t = 0; t < T; t++) {
                    dp[i][j][t].resize(4, 0);
                    caled[i][j][t].resize(4, false);
                }
            }
        }
        double tmp = solve_fidelity_max(0, l - 1, T - 1, 0);
        res += (tmp <= 0.25 ? 0.25 : tmp);
        cnt++;
    }
    return res / cnt;
}
double DP::min_fidelity(int T, double _tao) {
    if(T == 1) return 0.97;
    tao = _tao;
    double res = 0;
    int cnt = 0;
    for(int l = min(T, 7); l <= min(T, 7); l++) {
        dp.clear();
        dp.resize(l);
        caled.clear();
        caled.resize(l);
        for(int i = 0; i < l; i++) {
            dp[i].resize(l);
            caled[i].resize(l);
            for(int j = 0; j < l; j++) {
                dp[i][j].resize(T);
                caled[i][j].resize(T);
                for(int t = 0; t < T; t++) {
                    dp[i][j][t].resize(4, 0);
                    caled[i][j][t].resize(4, false);
                }
            }
        }
        double tmp = solve_fidelity_min(0, l - 1, T - 1, 0);
        res += (tmp > 10 ? 0.25 : tmp);
        cnt++;
    }
    return res / cnt;
}

// state = 0, left right no limit
// state = 1, left limit
// state = 2, right limit
// state = 3, left and right limit
double DP::solve_fidelity_max(int left, int right, int t, int state) {
    int left_remain = INF, right_remain = INF;
    if(state == 0 && (left_remain <= 0 || right_remain <= 0)) return 0;
    if(state == 1 && (left_remain <= 1 || right_remain <= 0)) return 0;
    if(state == 2 && (left_remain <= 0 || right_remain <= 1)) return 0;
    if(state == 3 && (left_remain <= 1 || right_remain <= 1)) return 0;
    
    if(t <= 0) return 0;
    if(left == right - 1) {
        int left_last_remain = INF;
        int right_last_remain = INF;
        if(state == 0 && (left_last_remain <= 0 || right_last_remain <= 0)) return 0;
        if(state == 1 && (left_last_remain <= 1 || right_last_remain <= 0)) return 0;
        if(state == 2 && (left_last_remain <= 0 || right_last_remain <= 1)) return 0;
        if(state == 3 && (left_last_remain <= 1 || right_last_remain <= 1)) return 0;
        return pass_tao(0.97, tao);
    }
    if(caled[left][right][t][state]) return dp[left][right][t][state];
    double best = pass_tao(solve_fidelity_max(left, right, t - 1, state), tao);
    for(int k = left + 1; k < right; k++) {
        for(int s = 1; s <= 2; s++) {
            int l_state = (state & 1) | ((s & 1) << 1);
            int r_state = (state & 2) | ((s & 2) >> 1);

            double left_result = solve_fidelity_max(left, k, t - 1, l_state);
            double right_result = solve_fidelity_max(k, right, t - 1, r_state);
            double result = Fswap(pass_tao(left_result, tao), pass_tao(right_result, tao));
            best = max(best, result);
        }
    }
    caled[left][right][t][state] = true;
    return dp[left][right][t][state] = best;
    return best;
}

double DP::solve_fidelity_min(int left, int right, int t, int state) {
    int left_remain = INF, right_remain = INF;
    if(state == 0 && (left_remain <= 0 || right_remain <= 0)) return INF;
    if(state == 1 && (left_remain <= 1 || right_remain <= 0)) return INF;
    if(state == 2 && (left_remain <= 0 || right_remain <= 1)) return INF;
    if(state == 3 && (left_remain <= 1 || right_remain <= 1)) return INF;
    
    if(t <= 0) return INF;
    if(left == right - 1) {
        int left_last_remain = INF;
        int right_last_remain = INF;
        if(state == 0 && (left_last_remain <= 0 || right_last_remain <= 0)) return INF;
        if(state == 1 && (left_last_remain <= 1 || right_last_remain <= 0)) return INF;
        if(state == 2 && (left_last_remain <= 0 || right_last_remain <= 1)) return INF;
        if(state == 3 && (left_last_remain <= 1 || right_last_remain <= 1)) return INF;
        return pass_tao(0.97, tao, false);
    }
    if(caled[left][right][t][state]) return dp[left][right][t][state];
    double best = pass_tao(solve_fidelity_min(left, right, t - 1, state), tao, false);
    for(int k = left + 1; k < right; k++) {
        for(int s = 1; s <= 2; s++) {
            int l_state = (state & 1) | ((s & 1) << 1);
            int r_state = (state & 2) | ((s & 2) >> 1);

            double left_result = solve_fidelity_min(left, k, t - 1, l_state);
            double right_result = solve_fidelity_min(k, right, t - 1, r_state);
            double result = Fswap(pass_tao(left_result, tao, false), pass_tao(right_result, tao , false), false);
            best = min(best, result);
        }
    }
    caled[left][right][t][state] = true;
    return dp[left][right][t][state] = best;
}

int main() {
    // double F = 0.97;
    // for(int i = 0; i < 4; i++) {
    //     F = pass_tao(F, 1e-2);
    // }
    // cerr << pass_tao(F, 1e-2) << "\n";
    // exit(0);

    DP solver;
    for(int T = 1; T <= 15; T++) {
        for(double tao = 1e-3; tao <= 2.5e-3 + EPS; tao += 2e-4) {
            double F_min = solver.min_fidelity(T, tao), F_max = solver.max_fidelity(T, tao);
            cerr << fixed << setprecision(5) << T << " " << tao << " ";
            cerr << fixed << setprecision(10) << F_min << " " << F_max << "\n";
            cout << fixed << setprecision(5) << T << " " << tao << " " << F_max / F_min << "\n";
        }
    }
}