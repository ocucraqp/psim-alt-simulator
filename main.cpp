#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cstring>
#include <queue>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <bitset>


using namespace std;

static const int PORT = 3;
static const double L[] = {37 * pow(10, -6), 37 * pow(10, -6), 37 * pow(10, -6)};
static const double F_SW = 20000;
static const double T_PERIOD = 1 / F_SW;

vector<pair<double, int>> calc_timing(vector<double> phi, vector<double> delta) {
    // 電圧変動時間（周期中のタイミング）の計算
    vector<pair<double, int>> T;
    for (int i = 0; i < PORT; i++) {
        vector<double> t(4, 0);
        t[0] = -(M_PI / 2) - (delta[i] / 2) + phi[i];
        t[1] = -(M_PI / 2) + (delta[i] / 2) + phi[i];
        t[2] = (M_PI / 2) - (delta[i] / 2) + phi[i];
        t[3] = (M_PI / 2) + (delta[i] / 2) + phi[i];
        for (int j = 0; j < 4; ++j) {
            // -pi~piの範囲外のものを補正
            if (t[j] >= M_PI) {
                t[j] -= 2 * M_PI;
            } else if (t[j] < -M_PI) {
                t[j] += 2 * M_PI;
            }
            T.emplace_back(t[j], i);
        }
    }

    T.emplace_back(M_PI, -1);
    T.emplace_back(-M_PI, -1);

    sort(T.begin(), T.end());

    return T;
}

vector<vector<pair<double, double>>>
calc_ui(vector<double> v, vector<double> phi, vector<double> delta, vector<pair<double, int>> T) {
    vector<vector<pair<double, double>>> ui(PORT);

    for (int i = 0; i < PORT; i++) {
        // 各ポートの電圧変動タイミングの取得
        vector<double> t;
        for (auto &j : T) {
            if (j.second == i) {
                t.emplace_back(j.first);
            }
        }
        sort(t.begin(), t.end());

        // 電圧計算
        if (phi[i] <= -(M_PI / 2 - delta[i] / 2)) {
            // 左へシフト
            ui[i].emplace_back(0, t[0]);
            ui[i].emplace_back(v[i], t[1]);
            ui[i].emplace_back(0, t[2]);
            ui[i].emplace_back(-v[i], t[3]);
            ui[i].emplace_back(0, M_PI);
        } else if (phi[i] >= M_PI / 2 - delta[i] / 2) {
            // 右へシフト
            ui[i].emplace_back(0, t[0]);
            ui[i].emplace_back(-v[i], t[1]);
            ui[i].emplace_back(0, t[2]);
            ui[i].emplace_back(v[i], t[3]);
            ui[i].emplace_back(0, M_PI);
        } else {
            ui[i].emplace_back(-v[i], t[0]);
            ui[i].emplace_back(0, t[1]);
            ui[i].emplace_back(v[i], t[2]);
            ui[i].emplace_back(0, t[3]);
            ui[i].emplace_back(-v[i], M_PI);
        }
    }

    return ui;
}

vector<double>
calc_ut(vector<vector<pair<double, double>>> ui, vector<pair<double, int>> T) {
    vector<double> ut;

    for (auto &t : T) {
        // 時間ごとの各ポートの電圧を確認
        vector<double> u(3, 0);
        for (int i = 0; i < PORT; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (t.first <= ui[i][j].second) {
                    u[i] = ui[i][j].first;
                    break;
                }
            }
        }

        double u_t = (L[1] * L[2] * u[0] + L[0] * L[2] * u[1] + L[0] * L[1] * u[2]) /
                     (L[1] * L[2] + L[0] * L[2] + L[0] * L[1]);
        ut.emplace_back(u_t);
    }

    return ut;
}

vector<vector<double>>
calc_il(vector<pair<double, int>> T, vector<vector<pair<double, double>>> ui, vector<double> ut) {
    vector<vector<double>> il(PORT);

    // t=0のときの電圧
    vector<double> u_0(3, 0);
    for (int i = 0; i < PORT; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (0 <= ui[i][j].second) {
                u_0[i] = ui[i][j].first;
                break;
            }
        }
    }
    double ut_0 = (L[1] * L[2] * u_0[0] + L[0] * L[2] * u_0[1] + L[0] * L[1] * u_0[2]) /
                  (L[1] * L[2] + L[0] * L[2] + L[0] * L[1]);

    // t=0のときの電流
    vector<double> i_0(3, 0);
    for (int i = 0; i < PORT; ++i) {
        for (int j = 1; j < T.size(); ++j) {
            // T[j]時のuiを求める
            double ui_tj = 0;
            for (auto &ui_tmp : ui[i]) {
                if (T[j].first <= ui_tmp.second) {
                    ui_tj = ui_tmp.first;
                    break;
                }
            }

            // t=0で計算終了
            if (T[j].first >= 0) {
                i_0[i] += (u_0[i] - ut_0) * T_PERIOD * ((0 - T[j - 1].first) / (2 * M_PI));
                i_0[i] *= -(1 / (2 * L[i]));
                break;
            }

            // t<0のとき
            i_0[i] += (ui_tj - ut[j]) * T_PERIOD * ((T[j].first - T[j - 1].first) / (2 * M_PI));
        }
    }

    for (int i = 0; i < PORT; ++i) {
        il[i].emplace_back(i_0[i]);
        for (int j = 1; j < T.size(); ++j) {
            // T[j]時のuiを求める
            double ui_tj = 0;
            for (auto &ui_tmp : ui[i]) {
                if (T[j].first <= ui_tmp.second) {
                    ui_tj = ui_tmp.first;
                    break;
                }
            }

            double i_tmp = il[i][j - 1] +
                           (1 / L[i]) * (ui_tj - ut[j]) * T_PERIOD * ((T[j].first - T[j - 1].first) / (2 * M_PI));
            il[i].emplace_back(i_tmp);
        }
    }

    return il;
}

vector<double> calc_il_peak(vector<vector<double>> il) {
    vector<double> il_peak(PORT, 0);

    for (int i = 0; i < PORT; ++i) {
        for (int j = 0; j < il[i].size(); ++j) {
            il_peak[i] = max(il_peak[i], il[i][j]);
        }
    }

    return il_peak;
}

void calc() {
    // input
    vector<double> v(3, 0), phi(3, 0), delta(3, 0);
//    for (int i = 0; i < PORT; ++i) {
//        cin >> v[i];
//    }
//    for (int i = 0; i < PORT; ++i) {
//        cin >> phi[i];
//    }
//    for (int i = 0; i < PORT; ++i) {
//        cin >> delta[i];
//    }
    v[0] = 400, v[1] = 400, v[2] = 200;
    phi[0] = 0, phi[1] = 30, phi[2] = 40;
    delta[2] = 70, delta[0] = delta[2] * v[2] / v[0], delta[1] = delta[2] * v[2] / v[1];
    // phiとdeltaはラジアンに変換
    for (int i = 0; i < PORT; ++i) {
        phi[i] *= M_PI / 180;
        delta[i] *= M_PI / 180;
    }

    // T: vector<pair<時間，ポート番号>>
    vector<pair<double, int>> T;
    T = calc_timing(phi, delta);

    // ui: vector<pair<電圧，時間>>
    vector<vector<pair<double, double>>> ui;
    ui = calc_ui(v, phi, delta, T);

    // ut: vector<電圧>
    vector<double> ut;
    ut = calc_ut(ui, T);

    // il: vector<vector<電流>>
    vector<vector<double>> il;
    il = calc_il(T, ui, ut);

    // il_peak: vector<電流>　
    vector<double> il_peak;
    il_peak = calc_il_peak(il);
};

int main() {
    calc();

    return 0;
}
