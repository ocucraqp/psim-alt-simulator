# alt-simulator

## Usage

### calc: Calculate il.

```bash
$ ./alt_simulator calc -h
undefined short option: -h
usage: ./alt_simulator [options] ... 
options:
  -f, --file      input file name (string [=input.csv])
  -p, --power     power
  -d, --delta3    delta3 only
  -m, --min       min
  -?, --help      print this message

$ ./alt_simulator calc
v[0]:400
v[1]:400
v[2]:200
delta[0]:130
delta[1]:130
delta[2]:80
phi[0]:0
phi[1]:30
phi[2]:40
Input
V1:400, V2:400, V3:200
phi1:0, phi2:30, phi3:40
delta1:130, delta2:130, delta3:80

Output
il1_peak:31.2813, il2_peak:18.7688, il3_peak:17.5175
il1_rms:14.6058, il2_rms:7.66082, il3_rms:9.85616
```

### phitop: Calculate power from phase difference.

```bash
$ ./alt_simulator phitop
v[0]:400
v[1]:400
v[2]:200
delta[0]:130
delta[1]:130
delta[2]:80
phi[0]:0 
phi[1]:39
phi[2]:33
POWER
PORT1:1864.78, PORT2: 1202.39, PORT3: 662.386

```

### ptophi: Calculate phase difference from power.

```bash
$ ./alt_simulator ptophi
v[0]:400
v[1]:400
v[2]:200
delta[0]:130
delta[1]:130
delta[2]:80
p[0]:2000
p[1]:1300
p[2]:700
phi
phi1:0, phi2: 39.1761, phi3: 32.9122

```

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