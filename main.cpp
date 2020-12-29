#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include "alt_simulator.h"
#include "cmdline.h"

using namespace std;

vector<double> phi_to_power(vector<double> v, vector<double> phi, vector<double> delta) {
    double L_s = L[0] * L[1] + L[1] * L[2] + L[2] * L[0];
    double A = (4 * T_PERIOD) / (pow(M_PI, 3) * L_s);
    double A_0 = A * L[0] * v[0] * v[1] * cos(delta[0] / 2) * cos(delta[1] / 2);
    double A_1 = A * L[1] * v[0] * v[2] * cos(delta[0] / 2) * cos(delta[2] / 2);
    double A_2 = A * L[2] * v[1] * v[2] * cos(delta[1] / 2) * cos(delta[2] / 2);
    vector<double> p(PORT);
    p[0] = A_0 * sin(phi[1]) + A_1 * sin(phi[2]);
    p[1] = A_0 * sin(phi[1]) - A_2 * sin(phi[2] - phi[1]);
    p[2] = A_1 * sin(phi[2]) + A_2 * sin(phi[2] - phi[1]);

    cout << "POWER" << endl;
    cout << "PORT1:" << p[0] << ", PORT2: " << p[1] << ", PORT3: " << p[2] << endl;

    return p;
}

vector<double> power_to_phi(vector<double> v, vector<double> p, vector<double> delta) {
    vector<double> phi(PORT, 0);
    vector<double> v_delta(PORT);
    for (int i = 0; i < PORT; ++i) {
        v_delta[i] = v[i] * cos(delta[i] / 2);
    }

    double A = (4 * T_PERIOD) / (pow(M_PI, 3) * (L[0] * L[1] + L[1] * L[2] + L[2] * L[0]));
    vector<double> G(4), H(4);
    G[0] = A * cos(delta[1] / 2) * (v_delta[0] * L[2] + v_delta[2] * L[0]);
    G[1] = -A * L[0] * cos(delta[1] / 2) * v_delta[2];
    G[2] = -A * L[0] * v_delta[1] * cos(delta[2] / 2);
    G[3] = A * cos(delta[2] / 2) * (v_delta[0] * L[1] + v_delta[1] * L[0]);
    double G_c = 1 / (G[0] * G[3] - G[1] * G[2]);
    H[0] = G_c * G[3];
    H[1] = -G_c * G[1];
    H[2] = -G_c * G[2];
    H[3] = G_c * G[0];
    phi[1] = H[0] * p[1] / v[1] + H[1] * p[2] / v[2];
    phi[2] = H[2] * p[1] / v[1] + H[3] * p[2] / v[2];

    cout << "phi" << endl;
    cout << "phi1:" << phi[0] * 180 / M_PI << ", phi2: " << phi[1] * 180 / M_PI << ", phi3: " << phi[2] * 180 / M_PI
         << endl;

    return phi;
}

void calc_by_input_std(AltSimulator alt_simulator, bool power = false, bool delta3_only = false) {
    vector<double> v(PORT), p(PORT), phi(PORT), delta(PORT);

    // 電圧
    for (int i = 0; i < PORT; ++i) {
        cin >> v[i];
    }

    // 電圧or位相差
    if (power) {
        for (int i = 0; i < PORT; ++i) {
            cin >> p[i];
        }
        // ラジアンに変換
        for (int i = 0; i < PORT; ++i) {
            p[i] *= M_PI / 180;
        }
        phi = power_to_phi(v, p, delta);
    } else {
        for (int i = 0; i < PORT; ++i) {
            cin >> phi[i];
        }
        // ラジアンに変換
        for (int i = 0; i < PORT; ++i) {
            phi[i] *= M_PI / 180;
        }
    }

    // ゼロ電圧動作区間
    if (delta3_only) {
        cin >> delta[2];
        delta[2] *= M_PI / 180;
        delta[0] = M_PI - (M_PI - delta[2]) * v[2] / v[0];
        delta[1] = M_PI - (M_PI - delta[2]) * v[2] / v[1];
    } else {
        for (int i = 0; i < PORT; ++i) {
            cin >> delta[i];
        }
        // ラジアンに変換
        for (int i = 0; i < PORT; ++i) {
            delta[i] *= M_PI / 180;
        }
    }

    alt_simulator.set_condition(v, phi, delta);
    alt_simulator.calc();
}

void calc_by_input_csv(AltSimulator alt_simulator, const string &input_filename, const string &output_filename,
                       bool power = false,
                       bool delta3_only = false) {
    string str_buf;
    string str_conma_buf;
    ifstream ifs_csv_file(input_filename);
    ofstream ofs_csv_file(output_filename);

    while (getline(ifs_csv_file, str_buf)) {
        vector<double> v(PORT), p(PORT), phi(PORT), delta(PORT);

        istringstream i_stream(str_buf);

        // 電圧
        for (int i = 0; i < PORT; ++i) {
            getline(i_stream, str_conma_buf, ',');
            v[i] = stod(str_conma_buf);
        }

        // 電圧or位相差
        if (power) {
            for (int i = 0; i < PORT; ++i) {
                getline(i_stream, str_conma_buf, ',');
                p[i] = stod(str_conma_buf);
            }
            // ラジアンに変換
            for (int i = 0; i < PORT; ++i) {
                p[i] *= M_PI / 180;
            }
            phi = power_to_phi(v, p, delta);
        } else {
            for (int i = 0; i < PORT; ++i) {
                getline(i_stream, str_conma_buf, ',');
                phi[i] = stod(str_conma_buf);
            }
            // ラジアンに変換
            for (int i = 0; i < PORT; ++i) {
                phi[i] *= M_PI / 180;
            }
        }

        // ゼロ電圧動作区間
        if (delta3_only) {
            getline(i_stream, str_conma_buf, ',');
            delta[2] = stod(str_conma_buf);
            delta[2] *= M_PI / 180;
            delta[0] = M_PI - (M_PI - delta[2]) * v[2] / v[0];
            delta[1] = M_PI - (M_PI - delta[2]) * v[2] / v[1];
        } else {
            for (int i = 0; i < PORT; ++i) {
                getline(i_stream, str_conma_buf, ',');
                delta[i] = stod(str_conma_buf);
            }
            // ラジアンに変換
            for (int i = 0; i < PORT; ++i) {
                delta[i] *= M_PI / 180;
            }
        }

        alt_simulator.set_condition(v, phi, delta);
        alt_simulator.calc();

        ofs_csv_file << str_conma_buf << ',';
        ofs_csv_file << std::endl;
    }
}


int main(int argc, char *argv[]) {
    // コマンドライン引数の処理
    cmdline::parser parser;
    parser.add<string>("file", 'f', "input file name", false, "input.csv");
    parser.add<string>("output", 'o', "output file name", false, "output.csv");
    parser.add("power", 'p', "power");
    parser.add("delta3", 'd', "delta3 only");
    parser.parse_check(argc, argv);

    if (!parser.parse(argc, argv) || parser.exist("help")) {
        std::cout << parser.error_full() << parser.usage();
        return 0;
    }

    AltSimulator alt_simulator;

    if (parser.exist("file")) {
        calc_by_input_csv(alt_simulator, parser.get<string>("file"), parser.get<string>("output"),
                          parser.exist("power"), parser.exist("delta3"));
    } else {
        calc_by_input_std(alt_simulator, parser.exist("power"), parser.exist("delta3"));
    }

    return 0;
}
