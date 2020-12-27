#ifndef ALT_SIMULATOR_ALT_SIMULATOR_H
#define ALT_SIMULATOR_ALT_SIMULATOR_H

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
static const int PARTITION_N = 1000;

class AltSimulator {
    vector<double> v, phi, delta;
    // T: vector<pair<時間，ポート番号>>
    vector<pair<double, int>> T;
    // ui: vector<vector<pair<電圧，時間>>>
    vector<vector<pair<double, double>>> ui;
    // ut: vector<pair<電圧，時間>>
    vector<pair<double, double>> ut;
    // il_0: vector<電流>
    vector<double> il_0;
    // il: vector<vector<電流>>
    vector<vector<double>> il;
    // il_peak: vector<ピーク電流>　
    vector<double> il_peak;
    // il_rms: vector<実効電流>　
    vector<double> il_rms;

public:
    void set_condition(vector<double> input_v, vector<double> input_phi, vector<double> input_delta);

    void calc_timing();

    void calc_ui();

    void calc_ut();

    double get_u_t(vector<pair<double, double>> u, double t);

    void calc_il_0();

    void calc_il();

    void calc_il_peak();

    void calc_il_rms();

    void calc();
};

#endif //ALT_SIMULATOR_ALT_SIMULATOR_H
