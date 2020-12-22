# alt-simulator

## variable memo
```c++
// v: vector<電圧>
vector<double> v(3, 0);
// v: vector<位相差>
vector<double> phi(3, 0);
// v: vector<>
vector<double> delta(3, 0);
// T: vector<pair<時間，ポート番号>>
vector<pair<double, int>> T;
// ui: pair<電圧，時間>
vector<vector<pair<double, double>>> ui(PORT);
// ut: vector<電圧>
vector<double> ut;
// il: vector<vector<電流>>
vector<vector<double>> il(PORT);
// il_peak: vector<電流>
vector<double> il_peak(PORT, 0);
```