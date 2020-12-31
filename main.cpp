#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include "alt_simulator.h"
#include "cmdline.h"

using namespace std;
static const int MAX_DELTA = 180;
static const int WIDTH_DELTA = 1;

template<typename T>
T max_vector(vector<T> v) {
    double max_v = v[0];
    for (int i = 1; i < v.size(); ++i) {
        max_v = max(max_v, v[i]);
    }
    return max_v;
}

template<typename T>
T avg_vector(vector<T> v) {
    double avg_v = 0;
    for (int i = 0; i < v.size(); ++i) {
        avg_v += v[i];
    }
    avg_v /= v.size();
    return avg_v;
}

vector<double> phi_to_power(vector<double> v, vector<double> delta, vector<double> phi) {
    double L_s = L[0] * L[1] + L[1] * L[2] + L[2] * L[0];
    double A = (4 * T_PERIOD) / (pow(M_PI, 3) * L_s);
    double A_0 = A * L[0] * v[0] * v[1] * cos(delta[0] / 2) * cos(delta[1] / 2);
    double A_1 = A * L[1] * v[0] * v[2] * cos(delta[0] / 2) * cos(delta[2] / 2);
    double A_2 = A * L[2] * v[1] * v[2] * cos(delta[1] / 2) * cos(delta[2] / 2);
    vector<double> p(PORT);
    p[0] = A_0 * sin(phi[1]) + A_1 * sin(phi[2]);
    p[1] = A_0 * sin(phi[1]) - A_2 * sin(phi[2] - phi[1]);
    p[2] = A_1 * sin(phi[2]) + A_2 * sin(phi[2] - phi[1]);

    return p;
}

vector<double> power_to_phi(vector<double> v, vector<double> delta, vector<double> p) {
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

    return phi;
}

vector<vector<vector<double>>> calc_delta3_for_min(AltSimulator alt_simulator, vector<double> v, vector<double> p) {
    vector<vector<double>> min_avg_peak(3), min_max_peak(3), min_avg_rms(3), min_max_rms(3);
    // ex.) min_avg_peak[0]:delta, min_avg_peak[1]:il_peak, min_avg_peak[2]:il_rms
    min_avg_peak[1] = {1000, 1000, 1000};
    min_avg_peak[2] = {1000, 1000, 1000};
    min_max_peak[1] = {1000, 1000, 1000};
    min_max_peak[2] = {1000, 1000, 1000};
    min_avg_rms[1] = {1000, 1000, 1000};
    min_avg_rms[2] = {1000, 1000, 1000};
    min_max_rms[1] = {1000, 1000, 1000};
    min_max_rms[2] = {1000, 1000, 1000};

    for (int i = 0; i < MAX_DELTA; i += WIDTH_DELTA) {
        vector<double> delta(PORT), phi(PORT);
        delta[2] = i * M_PI / 180;
        delta[0] = M_PI - (M_PI - delta[2]) * v[2] / v[0];
        delta[1] = M_PI - (M_PI - delta[2]) * v[2] / v[1];
        phi = power_to_phi(v, delta, p);

        alt_simulator.set_condition(v, delta, phi);
        auto[il_peak, il_rms]=alt_simulator.calc(false);

        // 旧データとの比較
        if (avg_vector(il_peak) < avg_vector(min_avg_peak[1])) {
            min_avg_peak = {delta, il_peak, il_rms};
        }
        if (max_vector(il_peak) < max_vector(min_max_peak[1])) {
            min_max_peak = {delta, il_peak, il_rms};
        }
        if (avg_vector(il_rms) < avg_vector(min_avg_rms[2])) {
            min_avg_rms = {delta, il_peak, il_rms};
        }
        if (max_vector(il_rms) < max_vector(min_max_rms[2])) {
            min_max_rms = {delta, il_peak, il_rms};
        }
    }

    return {min_avg_peak, min_max_peak, min_avg_rms, min_avg_rms};
}

void output_parameter(vector<vector<double>> min_delta) {
    cout << "delta1:" << min_delta[0][0] * 180 / M_PI << ", delta2:" << min_delta[0][1] * 180 / M_PI
         << ", delta3:" << min_delta[0][2] * 180 / M_PI << endl;
    cout << "il1_peak:" << min_delta[1][0] << ", il2_peak:" << min_delta[1][1] << ", il3_peak:"
         << min_delta[1][2] << endl;
    cout << "il1_rms:" << min_delta[2][0] << ", il2_rms:" << min_delta[2][1] << ", il3_rms:"
         << min_delta[2][2] << endl;
}

string output_parameter_to_csv(vector<vector<double>> min_delta) {
    string output = to_string(min_delta[0][0] * 180 / M_PI) + "," + to_string(min_delta[0][1] * 180 / M_PI) + "," +
                    to_string(min_delta[0][2] * 180 / M_PI)
                    + "," + to_string(min_delta[1][0]) + "," + to_string(min_delta[2][0]) + "," +
                    to_string(min_delta[1][1]) + "," + to_string(min_delta[2][1])
                    + "," + to_string(min_delta[1][2]) + "," + to_string(min_delta[2][2]) + '\n';
    return output;
}

void calc_by_input_std(AltSimulator alt_simulator, bool power = false, bool delta3_only = false, bool min = false) {
    vector<double> v(PORT), delta(PORT), p(PORT), phi(PORT);

    // 電圧
    for (int i = 0; i < PORT; ++i) {
        cout << "v[" << i << "]:";
        cin >> v[i];
    }

    // 最小ピーク電流のdeltaを計算
    if (min) {
        for (int i = 0; i < PORT; ++i) {
            cout << "p[" << i << "]:";
            cin >> p[i];
        }
        vector<vector<vector<double>>> min_delta = calc_delta3_for_min(alt_simulator, v, p);
        cout << "min avg peak:" << endl;
        output_parameter(min_delta[0]);
        cout << "min max peak:" << endl;
        output_parameter(min_delta[1]);
        cout << "min avg rms:" << endl;
        output_parameter(min_delta[2]);
        cout << "min max rms:" << endl;
        output_parameter(min_delta[3]);
    } else {
        // ゼロ電圧動作区間
        if (delta3_only) {
            cout << "delta[" << 2 << "]:";
            cin >> delta[2];
            delta[2] *= M_PI / 180;
            delta[0] = M_PI - (M_PI - delta[2]) * v[2] / v[0];
            delta[1] = M_PI - (M_PI - delta[2]) * v[2] / v[1];
        } else {
            for (int i = 0; i < PORT; ++i) {
                cout << "delta[" << i << "]:";
                cin >> delta[i];
            }
            // ラジアンに変換
            for (int i = 0; i < PORT; ++i) {
                delta[i] *= M_PI / 180;
            }
        }

        // 電圧or位相差
        if (power) {
            for (int i = 0; i < PORT; ++i) {
                cout << "p[" << i << "]:";
                cin >> p[i];
            }
            // ラジアンに変換
            for (int i = 0; i < PORT; ++i) {
                p[i] *= M_PI / 180;
            }
            phi = power_to_phi(v, delta, p);
        } else {
            for (int i = 0; i < PORT; ++i) {
                cout << "phi[" << i << "]:";
                cin >> phi[i];
            }
            // ラジアンに変換
            for (int i = 0; i < PORT; ++i) {
                phi[i] *= M_PI / 180;
            }
        }

        alt_simulator.set_condition(v, delta, phi);
        alt_simulator.calc(true);
    }
}

void calc_by_input_csv(AltSimulator alt_simulator, const string &input_filename, bool power = false,
                       bool delta3_only = false, bool min = false) {
    string str_buf;
    string str_conma_buf;
    ifstream ifs_csv_file(input_filename);

    string output_filename = "output";
    ofstream ofs_csv_file(output_filename + ".csv");
    ofstream ofs_csv_file_avg_peak(output_filename + "_avg_peak.csv");
    ofstream ofs_csv_file_max_peak(output_filename + "_max_peak.csv");
    ofstream ofs_csv_file_avg_rms(output_filename + "_avg_rms.csv");
    ofstream ofs_csv_file_max_rms(output_filename + "_max_rms.csv");

    while (getline(ifs_csv_file, str_buf)) {
        vector<double> v(PORT), p(PORT), phi(PORT), delta(PORT);

        istringstream i_stream(str_buf);

        // 電圧
        for (int i = 0; i < PORT; ++i) {
            getline(i_stream, str_conma_buf, ',');
            if (min) {
                ofs_csv_file_avg_peak << str_conma_buf << ',';
                ofs_csv_file_max_peak << str_conma_buf << ',';
                ofs_csv_file_avg_rms << str_conma_buf << ',';
                ofs_csv_file_max_rms << str_conma_buf << ',';
            } else {
                ofs_csv_file << str_conma_buf << ',';
            }
            v[i] = stod(str_conma_buf);
        }

        // 最小ピーク電流のdeltaを計算
        if (min) {
            for (int i = 0; i < PORT; ++i) {
                getline(i_stream, str_conma_buf, ',');
                ofs_csv_file_avg_peak << str_conma_buf << ',';
                ofs_csv_file_max_peak << str_conma_buf << ',';
                ofs_csv_file_avg_rms << str_conma_buf << ',';
                ofs_csv_file_max_rms << str_conma_buf << ',';
                p[i] = stod(str_conma_buf);
            }
            vector<vector<vector<double>>> min_delta = calc_delta3_for_min(alt_simulator, v, p);
            ofs_csv_file_avg_peak << output_parameter_to_csv(min_delta[0]);
            ofs_csv_file_max_peak << output_parameter_to_csv(min_delta[1]);
            ofs_csv_file_avg_rms << output_parameter_to_csv(min_delta[2]);
            ofs_csv_file_max_rms << output_parameter_to_csv(min_delta[3]);
        } else {
            // ゼロ電圧動作区間
            for (int i = 0; i < PORT; ++i) {
                getline(i_stream, str_conma_buf, ',');
                delta[i] = stod(str_conma_buf);
            }
            // ラジアンに変換
            if (delta3_only) {
                delta[2] *= M_PI / 180;
                delta[0] = M_PI - (M_PI - delta[2]) * v[2] / v[0];
                delta[1] = M_PI - (M_PI - delta[2]) * v[2] / v[1];
            } else {
                for (int i = 0; i < PORT; ++i) {
                    delta[i] *= M_PI / 180;
                }
            }
            ofs_csv_file << delta[0] * 180 / M_PI << ',' << delta[1] * 180 / M_PI << ',' << delta[2] * 180 / M_PI
                         << ',';

            // 電圧or位相差
            if (power) {
                for (int i = 0; i < PORT; ++i) {
                    getline(i_stream, str_conma_buf, ',');
                    ofs_csv_file << str_conma_buf << ',';
                    p[i] = stod(str_conma_buf);
                }
                // ラジアンに変換
                for (int i = 0; i < PORT; ++i) {
                    p[i] *= M_PI / 180;
                }
                phi = power_to_phi(v, delta, p);
            } else {
                for (int i = 0; i < PORT; ++i) {
                    getline(i_stream, str_conma_buf, ',');
                    ofs_csv_file << str_conma_buf << ',';
                    phi[i] = stod(str_conma_buf);
                }
                // ラジアンに変換
                for (int i = 0; i < PORT; ++i) {
                    phi[i] *= M_PI / 180;
                }
            }

            alt_simulator.set_condition(v, delta, phi);
            auto[il_peak, il_rms] = alt_simulator.calc(false);

            ofs_csv_file << il_peak[0] << "," << il_rms[0] << "," << il_peak[1] << "," << il_rms[1] << "," << il_peak[2]
                         << "," << il_rms[2] << endl;
        }
    }
}

void phi_to_p_by_input() {
    vector<double> v(PORT), delta(PORT), phi(PORT), p(PORT);

    for (int i = 0; i < PORT; ++i) {
        cout << "v[" << i << "]:";
        cin >> v[i];
    }
    for (int i = 0; i < PORT; ++i) {
        cout << "delta[" << i << "]:";
        cin >> delta[i];
    }
    for (int i = 0; i < PORT; ++i) {
        cout << "phi[" << i << "]:";
        cin >> phi[i];
    }

    p = phi_to_power(v, delta, phi);

    cout << "POWER" << endl;
    cout << "PORT1:" << p[0] << ", PORT2: " << p[1] << ", PORT3: " << p[2] << endl;
}

void p_to_phi_by_input() {
    vector<double> v(PORT), delta(PORT), phi(PORT), p(PORT);

    for (int i = 0; i < PORT; ++i) {
        cout << "v[" << i << "]:";
        cin >> v[i];
    }
    for (int i = 0; i < PORT; ++i) {
        cout << "delta[" << i << "]:";
        cin >> delta[i];
    }
    for (int i = 0; i < PORT; ++i) {
        cout << "p[" << i << "]:";
        cin >> phi[i];
    }

    phi = power_to_phi(v, delta, p);

    cout << "phi" << endl;
    cout << "phi1:" << phi[0] * 180 / M_PI << ", phi2: " << phi[1] * 180 / M_PI << ", phi3: " << phi[2] * 180 / M_PI
         << endl;
}

int main(int argc, char *argv[]) {

    // コマンドの確認
    if (!strcmp(argv[1], "phitop")) {
        phi_to_p_by_input();
    } else if (!strcmp(argv[1], "ptophi")) {
        p_to_phi_by_input();
    } else if (!strcmp(argv[1], "calc")) {
        // オプションを定義
        cmdline::parser parser;
        parser.add<string>("file", 'f', "input file name", false, "input.csv");
        parser.add("power", 'p', "power");
        parser.add("delta3", 'd', "delta3 only");
        parser.add("min", 'm', "min");

        // オプションの処理
        parser.parse_check(argc, argv);
        if (!parser.parse(argc, argv) || parser.exist("help")) {
            std::cout << parser.error_full() << parser.usage();
            return 0;
        }

        AltSimulator alt_simulator;

        // ファイルor標準入力による分岐
        if (parser.exist("file")) {
            calc_by_input_csv(alt_simulator, parser.get<string>("file"), parser.exist("power"), parser.exist("delta3"),
                              parser.exist("min"));
        } else {
            calc_by_input_std(alt_simulator, parser.exist("power"), parser.exist("delta3"), parser.exist("min"));
        }
    } else {
        cout << "Please input command." << endl;
        cout << "- calc: Calculate il." << endl;
        cout << "- phitop: Calculate power from phase difference." << endl;
        cout << "- ptophi: Calculate phase difference from power." << endl;
    }

    return 0;
}
