#include <bits/stdc++.h>
//#include "MT.h"

using namespace std;
float inf = numeric_limits<float>::infinity();
using ll = long long;
const double radius = 100;
const double GAMMA = pow(10, 6.0 / 10.0);
const double end_time = 1000;
constexpr double PI = 3.14159265358979323846264338;
const double theta = pow(10, 0.4);
//const double delta = 2.0 / pass_loss_exponent;
const double lambda_BS = 1.;
const double POWER = 100;
const double Euler_constant = 0.5772156649;

//出力ファイルを開く
ofstream outputfile;

struct pos {
    double x;
    double y;
};

struct terminal {
    pair<double, double> pos;
    //double to_BS[num_BS];
    double dst_to_origin;
    int nearest_BS;
    //double gain;
    double coef;
    bool state; //True: in , False: out
    bool operator<( const terminal& right ) const {
            return coef == right.coef ? true : coef < right.coef;
    }
};


double urand(){
    double m, a;
    m = RAND_MAX + 1.0;
    a = (rand() + 0.5)/m;
    a = (rand() + a)/m;
    return (rand() + a)/m;
}

//乱数発生
double my_rand(double Min, double Max) {
    mt19937 mt{ random_device{}() };
    //random_device rd;
    //default_random_engine eng(rd());
    uniform_real_distribution<double> distr(Min, Max);
    return distr(mt);
}

//円形領域内の座標をランダムに取得
pair<double, double> coordinate() {
    //init_genrand(i);
    double r = radius * sqrt(urand()); //[0, 1)
    double theta = PI * (2. * urand() - 1.); //(-1, 1)
    while (theta < -PI || PI < theta) {
        theta = PI * (2. * urand() - 1.);
    }
    double x = r * cos(theta);
    double y = r * sin(theta);
    pair<double, double> pos = make_pair(x, y);
    return pos;
}

//2点間の距離を計算
double cal_dst(pair<double, double> pos1, pair<double, double> pos2) {
    return sqrt((pos1.first - pos2.first)*(pos1.first - pos2.first) + (pos1.second - pos2.second)*(pos1.second - pos2.second));
}

double distance(pos pos1, pos pos2) {
    return sqrt((pos1.x - pos2.x) * (pos1.x - pos2.x) + (pos1.y - pos2.y) * (pos1.y - pos2.y));
}

//正規分布乱数発生
double normrand(){
    return sqrt(-2*log(urand()))*cos(2*M_PI*urand());
}


double gauss_rand(double mu, double sig) {
    //mt19937 mt{ std::random_device{}() };
    random_device rd;
    default_random_engine eng(rd());
    normal_distribution<> dist(mu, sig);
    double g = dist(eng);
    return g;
}


double exp_dist(double lambda) {
    //double g = genrand_real3();
    double g = urand();
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


double theo2(double lambda_IoT, double lambda_BS, double alpha) {
    double delta = 2.0 / alpha;
    double n = lambda_IoT * PI * radius * radius;
    //double g = (n - log(n + Euler_constant)) * 4. / lambda_IoT;
    //cout << "g = " <<  g << endl;
    double A = pow(theta * Euler_constant, delta) * PI * delta / sin(PI * delta);
    double I = 1.0 / (1.0 +  A * lambda_IoT / lambda_BS);
    return I;
}

double thp_theo(double lambda_IoT, double lambda_BS, double alpha) {
    double delta = 2.0 / alpha;
    double A = pow(theta, delta) * PI * delta / sin(PI * delta);
    double I = 1.0 / (A + lambda_BS / lambda_IoT);
    return I;
}



void sim_pr2(double lambda_IoT, double alpha, double noise) {
    double success = 0;
    cout << "pass loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        int num_IoT = poisson_dist(4 * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(4 * radius * radius * lambda_BS);
        //BSの座標を設定
        //vector<pos> BS_pos(num_BS);
        pos BS_pos, nearest_pos;
        pos origin = {0, 0};
        double nd = inf;
        for (int i = 0; i < num_BS; i++) {
//            BS_pos.at(i).x = -radius + urand() * 2 * radius;
//            BS_pos.at(i).y = -radius + urand() * 2 * radius;
            BS_pos.x = -radius + urand() * 2 * radius;
            BS_pos.y = -radius + urand() * 2 * radius;
            double dst = distance(origin, BS_pos);
            if (dst < nd) {
                nd = dst;
                nearest_pos = BS_pos;
            }
        }
        
        double SI = 0, coef = 0.0;
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            pos IoT_pos;
            if (i == 0) {
                IoT_pos = origin;
            } else {
                IoT_pos.x = -radius + urand() * 2 * radius;
                IoT_pos.y = -radius + urand() * 2 * radius;
            }
            double dst = distance(origin, IoT_pos);
            if (radius < dst) continue;
            //端末-基地局間のフェージング係数を設定
            double H = exp_dist(1.0);//gauss_rand(0, 1);
            double dst2 = distance(IoT_pos, nearest_pos);
            if (i == 0) coef = H / pow(dst2, alpha);
            else SI += H / pow(dst2, alpha);
        }
//        pair<double, double> SI = ini(num_IoT, num_BS, alpha);
        double SIR = coef / (SI + noise);
        if (SIR > theta) success++;
        
        //進捗状況を表示
        double progress = i / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
    double res = success / end_time;
    double I = theo(lambda_IoT, lambda_BS, alpha);
    cout << res << " " << I << endl << endl;
    outputfile << res << " " << I << " ";
    //outputfile << res;
}


void sim_thp(double lambda_IoT, double alpha, double noise) {
    double success = 0.0, power_dif_sum = 0.0;
    cout << "pass loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        int num_IoT = poisson_dist(4 * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(4 * radius * radius * lambda_BS);
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
            double dst_to = cal_dst(pos, origin);
            device.at(i).dst_to_origin = dst_to;
            if (dst_to < 20) {
                for (int j = 0; j < num_BS; j++) {
                    double dst = cal_dst(pos, BS_pos.at(j));
                    device.at(i).state = true;
                    if (dst < device.at(i).dst_to_origin) {
                        device.at(i).state = false;
                        break;
                    }
                    if (j == num_BS - 1) acl.push_back(i);
                }
            }
            
            //端末-基地局間のフェージング係数を設定
            double H = exp_dist(1.0);//gauss_rand(0, 1);
            double coef = H / pow(device.at(i).dst_to_origin, alpha);
            if (device.at(i).state) device.at(i).coef = coef;
            SI += coef;
        }
        
        int s = (int)acl.size();
        for (int i = 0; i < s; i++) {
            double P = device.at(acl.at(i)).coef;
            double SIR = P / (SI - P + noise);
            if (SIR > theta) {success++; break;}
        }
        
        //進捗状況を表示
        double progress = i / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
    double res = success / end_time;
    double I = thp_theo(lambda_IoT, lambda_BS, alpha);
    cout << "Throughput : " <<  res << " " << I << endl << endl;
    //outputfile << res << " " << I << " ";
}


//NOMA
void NOMA_thp(double lambda_IoT, double alpha, double noise) {
    double success = 0.0;
    cout << "pass loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        int num_IoT = poisson_dist(4 * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(4 * radius * radius * lambda_BS);
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
            double dst_to = cal_dst(pos, origin);
            device.at(i).dst_to_origin = dst_to;
            if (dst_to < 20) {
                for (int j = 0; j < num_BS; j++) {
                    double dst = cal_dst(pos, BS_pos.at(j));
                    device.at(i).state = true;
                    if (dst < device.at(i).dst_to_origin) {
                        device.at(i).state = false;
                        break;
                    }
                    if (j == num_BS - 1) acl.push_back(i);
                }
            }
            
            //端末-基地局間のフェージング係数を設定
            double H = exp_dist(1.0);//gauss_rand(0, 1);
            double coef = H / pow(device.at(i).dst_to_origin, alpha);
            device.at(i).coef = coef;
            SI += coef;
        }
        
        int s = (int)acl.size();
        double big = 0, big2 = 0;
        for (int i = 0; i < s; i++) {
            double P = device.at(acl.at(i)).coef;
            if (big2 < P) {
                if (big < P) big = P;
                else big2 = P;
            }
        }
        
        SI -= big;
        double SINR = big / (SI + noise);
        if (SINR > theta) success++;
        SI -= big2;
        SINR = big2 / (SI + noise);
        if (SINR > theta) success++;
        
        //進捗状況を表示
        double progress = i / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
    double res = success / end_time;
    double I = thp_theo(lambda_IoT, lambda_BS, alpha);
    cout << "Throughput : " <<  res << " " << I << endl << endl;
    //outputfile << res << " " << I << " ";
}



int main() {
    cout << "1: Pr sim, 2: Thp sim, 3: NOMA : ";
    int key; cin >> key;
    cout << "Put start lambda IoT : ";
    double k; cin >> k; cout << endl;
    string filename;
    if (key == 1) filename = to_string(k) + "pr.txt";
    else filename = to_string(k) + "thp.txt";
    //outputfile.open(filename);
    for (double i = k; i <= 5; i += 0.5) {
        //outputfile << i << " ";
        cout << "lambda IoT : " << i << endl;
//        double noise = 0;
//        for (double alpha = 2.5; alpha <= 4.5; alpha += 0.5) {
//            if (key == 1) sim_pr2(i, alpha, noise);
//            else if (key == 2) sim_thp(i, alpha, noise);
//            else NOMA_thp(i, alpha, noise);
//        }
        //alpha固定、ノイズ変化
        double alpha = 4.0;
        for (double noise = 10; noise <= 100000; noise *= 10) {
            if (key == 1) sim_pr2(i, alpha, noise);
            else if (key == 2) sim_thp(i, alpha, noise);
            else NOMA_thp(i, alpha, noise);
        }
        cout << endl;
        //outputfile << endl;
    }
    
//    double k = 0.5;
//    for (double i = k; i <= 5; i += 0.5) {
//        cout << "lambda IoT : " << i << endl;
//        for (double alpha = 2.5; alpha <= 4.5; alpha += 0.5) {
//            double I = theo(i, lambda_BS, alpha);
//            double I2 =  theo2(i, lambda_BS, alpha);
//            cout << i * (I + I2) / lambda_BS << endl;
//        }
//        cout << endl;
//    }
    
    //max throughput
//    string filename = "max_thp2.txt";
//    outputfile.open(filename);
//    double i = 10.;
//    for (double alpha = 2; alpha <= 10; alpha += 0.5) {
//        outputfile << alpha << " ";
//        sim_thp(i, alpha);
//        //cout << alpha << " " << thp_theo(i, lambda_BS, alpha) << " " << endl;
//    }
//    cout << endl;
//    outputfile << endl;
    
    
//    random_device rd;
//    init_genrand(rd());
//    for (int i = 0; i < 1000; i++) {
//        pair<double, double> test = coordinate(i);
//        cout << test.first << " " << test.second << endl;
//    }
//
    
    
//    for (double alpha = 2; alpha <= 10; alpha += 0.5) {
//        double delta = 2.0 / alpha;
//        double A = pow(theta, delta) * PI * delta / sin(PI * delta);
//        cout << alpha << " " << 1/A << endl;
//        //cout << "alpha = " << alpha << " , A = " << A << ", 1/A = " << 1. / A << endl;
//    }
    
    //A 最大スループット
//    for (double i = 1; i < 1000000; i+=100) {
//        double alpha = 4.5;
//        double delta = 2.0 / alpha;
//        double A = pow(theta, delta) * PI * delta / sin(PI * delta);
//        double I = theo(i, lambda_BS, alpha);
//        double T = thp_theo(i, lambda_BS, alpha);
//        cout << i << " " << T << " " << I * i << " " << 1. / A << endl;
//    }
    
    
}
