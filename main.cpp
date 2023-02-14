#include <bits/stdc++.h>

using namespace std;
float inf = numeric_limits<float>::infinity();
using ll = long long;
const double radius = 100;
const double GAMMA = pow(10, 6.0 / 10.0);
const double end_time = 100;
constexpr double PI = 3.14159265358979323846264338;
const int num_channel = 6;
const int num_power = 4;
const double pass_loss_exponent = 2.5;
const double POWER = 1.0;
const double theta = pow(10, 0.4);
//const double delta = 2.0 / pass_loss_exponent;
const double lambda_BS = 1.0;

//出力ファイルを開く
ofstream outputfile;

struct terminal {
    pair<double, double> pos;
    //double to_BS[num_BS];
    double nearest_dst;
    int nearest_BS;
    //double gain;
    double coef;
    bool state; //True: in , False: out
};

struct BS {
    int channel[num_channel];
    int power[num_power];
};

//int num_terminal = 10000;
//vector<terminal> T_vec(num_terminal);


double urand(){
    double m, a;
    m = RAND_MAX + 1.0;
    a = (rand() + 0.5)/m;
    a = (rand() + a)/m;
    return (rand() + a)/m;
}

//擬似乱数発生
double my_rand(double Min, double Max) {
    mt19937 mt{ std::random_device{}() };
    //random_device rd;
    //default_random_engine eng(rd());
    uniform_real_distribution<double> distr(Min, Max);
    double g = distr(mt);
    return g;
}

//円形領域内の座標をランダムに取得
pair<double, double> coordinate() {
    double r = radius * sqrt(my_rand(0, 1));
    double theta = my_rand(-PI, PI);
    double x = r * cos(theta);
    double y = r * sin(theta);
    pair<double, double> pos = make_pair(x, y);
    return pos;
}

//2点間の距離を計算
double cal_dst(pair<double, double> pos1, pair<double, double> pos2) {
    return sqrt((pos1.first - pos2.first)*(pos1.first - pos2.first) + (pos1.second - pos2.second)*(pos1.second - pos2.second));
}

//正規分布乱数発生
double gauss_rand(double mu, double sig) {
    //mt19937 mt{ std::random_device{}() };
    random_device rd;
    default_random_engine eng(rd());
    normal_distribution<> dist(mu, sig);
    //uniform_real_distribution<double> distr(Min, Max);
    double g = dist(eng);
    return g;
}


double exp_dist(double lambda) {
    double g = my_rand(0, 1);
    double tau = - log(1 - g) / lambda;
    return tau;
}

double poisson_dist(double lambda) {
    double sum = 0, k;
    for (k = 0; sum < 1; k++) {
        sum += exp_dist(lambda);
    }
    return  k;
}

//ポアソン乱数
double poisson(double lambda) {
    random_device rd;
    default_random_engine eng(rd());
    poisson_distribution<> dist(lambda);
    return dist(eng);
}


double theo(double lambda_IoT, double lambda_BS, double alpha) {
    double delta = 2.0 / alpha;
    double A = pow(theta, delta) * PI * delta / sin(PI * delta);
    double I = 1.0 / (1.0 +  A * lambda_IoT / lambda_BS);
    return I;
}


double thp_theo(double lambda_IoT, double lambda_BS, double alpha) {
    double delta = 2.0 / alpha;
    double A = pow(theta, delta) * PI * delta / sin(PI * delta);
    double I = 1.0 / (A + lambda_BS / lambda_IoT);
    return I;
}




void sim_pr(double lambda_IoT, double alpha) {
    int num_IoT = poisson_dist(4 * radius * radius * lambda_IoT);
    int num_BS = poisson_dist(4 * radius * radius * lambda_BS);
    
    double success = 0;
    cout << "pass loss exp : " << alpha << endl;
    
//    bool p_flag[4] = {true, true, true, true};
//    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        //BSの座標を設定
        vector<pair<double, double>> BS_pos(num_BS);
        pair<double, double> origin = make_pair(0, 0);
        int index = 0;
        double nd = inf;
        for (int i = 0; i < num_BS; i++) {
            BS_pos.at(i) = coordinate();
            double dst = cal_dst(origin, BS_pos.at(i));
            if (dst < nd) {
                nd = dst; index = i;
            }
        }
        
        double SI = 0, coef = 0.0;
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            pair<double, double> pos;
            if (i == 0) {
                pos = origin;
            } else pos = coordinate();
            
            //端末-基地局間のフェージング係数を設定
            double H = gauss_rand(0, 1);
            double dst = cal_dst(pos, BS_pos.at(index));
            if (i == 0) coef = H / pow(dst, alpha);
            else SI += H / pow(dst, alpha);
        }
//        pair<double, double> SI = ini(num_IoT, num_BS, alpha);
        double SIR = coef / SI;
        if (SIR > theta) success++;
        
        //進捗状況を表示
//        double progress = i / end_time;
//        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
//        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
//        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        //else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    //cout << "...100%" << endl;
    
    double res = success / end_time;
    double I = theo(lambda_IoT, lambda_BS, alpha);
    cout << res << " " << I << endl << endl;
    outputfile << res << " " << I << " ";
    //outputfile << res;
}


void sim_thp(double lambda_IoT, double alpha) {
    int num_IoT = poisson(4 * radius * radius * lambda_IoT);
    int num_BS = poisson(4 * radius * radius * lambda_BS);
    
    double success = 0;
    cout << "pass loss exp : " << alpha << endl;
    
//    bool p_flag[4] = {true, true, true, true};
//    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        //BSの座標を設定
        vector<pair<double, double>> BS_pos(num_BS);
        pair<double, double> origin = make_pair(0, 0);
        BS_pos.at(0) = origin;
        for (int i = 1; i < num_BS; i++) {
            BS_pos.at(i) = coordinate();
        }
        
        double SI = 0;
        vector<terminal> device(num_IoT);
        vector<int> acl;
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            pair<double, double> pos = coordinate();
            device.at(i).nearest_dst = cal_dst(pos, origin);
            for (int j = 0; j < num_BS; j++) {
                double dst = cal_dst(pos, BS_pos.at(j));
                device.at(i).state = true;
                if (dst < device.at(i).nearest_dst) {
                    device.at(i).state = false;
                    break;
                }
                if (j == num_BS - 1) acl.push_back(i);
            }
            
            //端末-基地局間のフェージング係数を設定
            double H = gauss_rand(0, 1);
            double coef = H / pow(device.at(i).nearest_dst, alpha);
            if (device.at(i).state) device.at(i).coef = coef;
            SI += coef;
        }
        for (int i = 0; i < acl.size(); i++) {
            double P = device.at(acl.at(i)).coef;
            double SIR = P / (SI - P);
            if (SIR > theta) {success++; break;}
        }
        
        //進捗状況を表示
//        double progress = i / end_time;
//        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
//        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
//        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        //else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    //cout << "...100%" << endl;
    
    double res = success / end_time;
    double I = thp_theo(lambda_IoT, lambda_BS, alpha);
    cout << res << " " << I << endl << endl;
    outputfile << res << " " << I << " ";
}


int main() {
    string filename = "SPALOHA0.txt";
    outputfile.open(filename);
    for (double i = 0.5; i <= 5; i += 0.5) {
        outputfile << i << " ";
        cout << "lambda IoT : " << i << endl;
        for (double alpha = 2.5; alpha <= 4.5; alpha += 0.5) {
            //sim_pr(i, alpha);
            sim_thp(i, alpha);
            //cout << alpha << " " << thp_theo(i, lambda_BS, alpha) << " " << endl;
        }
        cout << endl;
        outputfile << endl;
    }
}
