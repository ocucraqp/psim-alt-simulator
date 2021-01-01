#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "alt_simulator.h"
#include "make_dataset.h"

void make_dataset() {
    string filename;
    vector<double> v_min(PORT), v_max(PORT), v_width(PORT), p_min(PORT), p_max(PORT), p_width(PORT);

    cout << "filename:";
    cin >> filename;
    for (int i = 0; i < PORT; ++i) {
        cout << "v[" << i << "] min  :";
        cin >> v_min[i];
        cout << "v[" << i << "] max  :";
        cin >> v_max[i];
        cout << "v[" << i << "] width:";
        cin >> v_width[i];
    }
    for (int i = 1; i < PORT; ++i) {
        cout << "p[" << i << "] min  :";
        cin >> p_min[i];
        cout << "p[" << i << "] max  :";
        cin >> p_max[i];
        cout << "p[" << i << "] width:";
        cin >> p_width[i];
    }

    ofstream ofs_csv_file(filename);
    for (double v0 = v_min[0]; v0 <= v_max[0]; v0 += v_width[0]) {
        for (double v1 = v_min[1]; v1 <= v_max[1]; v1 += v_width[1]) {
            for (double v2 = v_min[2]; v2 <= v_max[2]; v2 += v_width[2]) {
                for (double p1 = p_min[1]; p1 <= p_max[1]; p1 += p_width[1]) {
                    for (double p2 = p_min[2]; p2 <= p_max[2]; p2 += p_width[2]) {
                        double p0 = p1 + p2;
                        ofs_csv_file << v0 << ',' << v1 << ',' << v2 << ',' << p0 << ',' << p1 << ',' << p2 << ','
                                     << endl;
                    }
                }
            }
        }
    }
}