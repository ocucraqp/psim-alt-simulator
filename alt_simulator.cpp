#include "alt_simulator.h"

void AltSimulator::set_condition(vector<double> input_v, vector<double> input_phi, vector<double> input_delta) {
    v = input_v;
    phi = input_phi;
    delta = input_delta;
}

void AltSimulator::calc_timing() {
    T.clear();

    // 電圧変動時間（周期中のタイミング）の計算
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
}

void AltSimulator::calc_ui() {
    ui.clear();
    ui.assign(PORT, vector<pair<double, double>>{});

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
}

void AltSimulator::calc_ut() {
    ut.clear();

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

        // 変換器の電圧を算出
        double u_t = (L[1] * L[2] * u[0] + L[0] * L[2] * u[1] + L[0] * L[1] * u[2]) /
                     (L[1] * L[2] + L[0] * L[2] + L[0] * L[1]);
        ut.emplace_back(u_t, t.first);
    }
}

double AltSimulator::get_u_t(vector<pair<double, double>> u, double t) {
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

void AltSimulator::calc_il_0() {
    il_0.clear();
    il_0.assign(PORT, 0);

    for (int i = 0; i < PORT; ++i) {
        for (int j = 1; j < T.size(); ++j) {
            // T[j]時のuiを求める
            double ui_tj = get_u_t(ui[i], T[j].first);

            // t=0で計算終了
            if (T[j].first >= 0) {
                il_0[i] += (ui_tj - ut[j].first) * T_PERIOD * ((0 - T[j - 1].first) / (2 * M_PI));
                il_0[i] *= -(1 / (2 * L[i]));
                break;
            }

            // t<0のとき
            il_0[i] += (ui_tj - ut[j].first) * T_PERIOD * ((T[j].first - T[j - 1].first) / (2 * M_PI));
        }
    }
}

void AltSimulator::calc_il() {
    il.clear();
    il.assign(PORT, vector<double>{});

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
}

void AltSimulator::calc_il_peak() {
    il_peak.clear();
    il_peak.assign(PORT, 0);

    for (int i = 0; i < PORT; ++i) {
        for (int j = 0; j < il[i].size(); ++j) {
            il_peak[i] = max(il_peak[i], il[i][j]);
        }
    }
}

void AltSimulator::calc_il_rms() {
    il_rms.clear();
    il_rms.assign(PORT, 0);

    for (int i = 0; i < PORT; ++i) {
        double il_detail = il_0[i];
        for (int j = 1; j < PARTITION_N; ++j) {
            double t = -M_PI + ((2 * M_PI / PARTITION_N) * j);
            double t_1 = -M_PI + ((2 * M_PI / PARTITION_N) * (j - 1));
            // t時のuiを求める
            double ui_t = get_u_t(ui[i], t);
            double ut_t = get_u_t(ut, t);
            // ilの計算
            il_detail += (1 / L[i]) * (ui_t - ut_t) * T_PERIOD * ((t - t_1) / (2 * M_PI));
            il_rms[i] += pow(il_detail, 2);
        }
        il_rms[i] = sqrt(il_rms[i] / PARTITION_N);
    }
}

void AltSimulator::calc() {
    calc_timing();
    calc_ui();
    calc_ut();
    calc_il_0();
    calc_il();
    calc_il_peak();
    calc_il_rms();

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