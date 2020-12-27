#include <iostream>
#include <vector>
#include <cstdlib>
#include "alt_simulator.h"
#include "cmdline.h"

using namespace std;


void calc_by_input_std(AltSimulator alt_simulator) {
    vector<double> v(PORT), phi(PORT), delta(PORT);

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

    alt_simulator.set_condition(v, phi, delta);
    alt_simulator.calc();
}

void calc_by_input_std_only_delta3(AltSimulator alt_simulator) {
    vector<double> v(PORT), phi(PORT), delta(PORT);

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

    alt_simulator.set_condition(v, phi, delta);
    alt_simulator.calc();
}

void calc_by_input_csv(AltSimulator alt_simulator) {
    vector<double> v(PORT), phi(PORT), delta(PORT);

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

    alt_simulator.set_condition(v, phi, delta);
    alt_simulator.calc();
}


int main(int argc, char *argv[]) {
    // コマンドライン引数の処理
    cmdline::parser parser;
    parser.add<string>("file", 'f', "file name", false, "");
    parser.add("power", 'p', "power");
    parser.add("only", 'o', "only port3");
    parser.parse_check(argc, argv);

    AltSimulator alt_simulator;

    if (parser.exist("only")) {
        calc_by_input_std_only_delta3(alt_simulator);
    } else {
        calc_by_input_std(alt_simulator);
    }

    return 0;
}
