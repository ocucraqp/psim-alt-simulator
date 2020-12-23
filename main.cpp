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

vector<pair<double, double>>
calc_ut(vector<vector<pair<double, double>>> ui, vector<pair<double, int>> T) {
    vector<pair<double, double>> ut;

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
        ut.emplace_back(u_t, t.first);
    }

    return ut;
}

double get_u_t(vector<pair<double, double>> u, double t) {
    // T[j]時のuiを求める
    double u_t = 0;
    for (auto &u_tmp : u) {
        if (t <= u_tmp.second) {
            u_t = u_tmp.first;
            break;
        }
    }
    return u_t;
}

vector<double>
calc_il_t(vector<pair<double, int>> T, vector<vector<pair<double, double>>> ui, vector<pair<double, double>> ut,
          double t) {
    vector<double> il_t(3, 0);
    for (int i = 0; i < PORT; ++i) {
        for (int j = 1; j < T.size(); ++j) {
            // T[j]時のuiを求める
            double ui_tj = get_u_t(ui[i], T[j].first);

            // t=0で計算終了
            if (T[j].first >= t) {
                il_t[i] += (ui_tj - ut[j].first) * T_PERIOD * ((0 - T[j - 1].first) / (2 * M_PI));
                il_t[i] *= -(1 / (2 * L[i]));
                break;
            }

            // t<0のとき
            il_t[i] += (ui_tj - ut[j].first) * T_PERIOD * ((T[j].first - T[j - 1].first) / (2 * M_PI));
        }
    }

    return il_t;
}

vector<vector<double>>
calc_il(vector<pair<double, int>> T, vector<vector<pair<double, double>>> ui, vector<pair<double, double>> ut,
        vector<double> il_0) {
    vector<vector<double>> il(PORT);

    for (int i = 0; i < PORT; ++i) {
        il[i].emplace_back(il_0[i]);
        for (int j = 1; j < T.size(); ++j) {
            // T[j]時のuiを求める
            double ui_tj = get_u_t(ui[i], T[j].first);

            double i_tmp = il[i][j - 1] +
                           (1 / L[i]) * (ui_tj - ut[j].first) * T_PERIOD * ((T[j].first - T[j - 1].first) / (2 * M_PI));
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

vector<double>
calc_il_rms(vector<vector<pair<double, double>>> ui, vector<pair<double, double>> ut, vector<double> il_0,
            int partition_n) {
    vector<double> il_rms(PORT, 0);

    for (int i = 0; i < PORT; ++i) {
        double il = il_0[i];
        for (int j = 1; j < partition_n; ++j) {
            double t = -M_PI + ((2 * M_PI / partition_n) * j);
            double t_1 = -M_PI + ((2 * M_PI / partition_n) * (j - 1));
            // t時のuiを求める
            double ui_t = get_u_t(ui[i], t);
            double ut_t = get_u_t(ut, t);
            // ilの計算
            il += (1 / L[i]) * (ui_t - ut_t) * T_PERIOD * ((t - t_1) / (2 * M_PI));
            il_rms[i] += pow(il, 2);
        }
        il_rms[i] = sqrt(il_rms[i] / partition_n);
    }

    return il_rms;
}

void calc(vector<double> v, vector<double> phi, vector<double> delta) {
    // T: vector<pair<時間，ポート番号>>
    vector<pair<double, int>> T;
    T = calc_timing(phi, delta);

    // ui: vector<pair<電圧，時間>>
    vector<vector<pair<double, double>>> ui;
    ui = calc_ui(v, phi, delta, T);

    // ut: vector<電圧>
    vector<pair<double, double>> ut;
    ut = calc_ut(ui, T);

    // il_0: vector<電流>
    vector<double> il_0;
    il_0 = calc_il_t(T, ui, ut, 0);

    // il: vector<vector<電流>>
    vector<vector<double>> il;
    il = calc_il(T, ui, ut, il_0);

    // il_peak: vector<ピーク電流>　
    vector<double> il_peak;
    il_peak = calc_il_peak(il);

    // il_rms: vector<実効電流>　
    vector<double> il_rms;
    il_rms = calc_il_rms(ui, ut, il_0, 100);

    // output
    cout << "Input" << endl;
    cout << "V1:" << v[0] << ", V2:" << v[1] << ", V3:" << v[2] << endl;
    cout << "phi1:" << phi[0] * 180 / M_PI << ", phi2:" << phi[1] * 180 / M_PI << ", phi3:" << phi[2] * 180 / M_PI
         << endl;
    cout << "delta1:" << delta[0] * 180 / M_PI << ", delta2:" << delta[1] * 180 / M_PI << ", delta3:"
         << delta[2] * 180 / M_PI << endl;
    cout << endl << "Output" << endl;
    cout << "il1_peak:" << il_peak[0] << ", il2_peak:" << il_peak[1] << ", il3_peak:" << il_peak[2] << endl;
    cout << "il1_rms:" << il_rms[0] << ", il2_rms:" << il_rms[1] << ", il3_rms:" << il_rms[2] << endl;
    cout << endl;

}

void calc_by_input_std() {
    vector<double> v(3, 0), phi(3, 0), delta(3, 0);

    for (int i = 0; i < PORT; ++i) {
        cin >> v[i];
    }
    for (int i = 0; i < PORT; ++i) {
        cin >> phi[i];
    }
    for (int i = 0; i < PORT; ++i) {
        cin >> delta[i];
    }

    // deltaとphiはラジアンに変換
    for (int i = 0; i < PORT; ++i) {
        phi[i] *= M_PI / 180;
        delta[i] *= M_PI / 180;
    }

    calc(v, phi, delta);
}

void calc_by_input_std_only_delta3() {
    vector<double> v(3, 0), phi(3, 0), delta(3, 0);

    // 電圧
    for (int i = 0; i < PORT; ++i) {
        cin >> v[i];
    }

    // 位相差
    for (int i = 0; i < PORT; ++i) {
        cin >> phi[i];
    }
    // ラジアンに変換
    for (int i = 0; i < PORT; ++i) {
        phi[i] *= M_PI / 180;
    }

    // ゼロ電圧動作区間
    cin >> delta[2];
    delta[2] *= M_PI / 180;
    delta[0] = M_PI - (M_PI - delta[2]) * v[2] / v[0];
    delta[1] = M_PI - (M_PI - delta[2]) * v[2] / v[1];

    calc(v, phi, delta);
}

int main() {
//    calc_by_input_std();
    calc_by_input_std_only_delta3();
    return 0;
}
