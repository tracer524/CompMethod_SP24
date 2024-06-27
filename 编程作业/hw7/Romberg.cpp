#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

using namespace std;

class Romberg{
public:
    Romberg(function<double(double)> func, double a, double b, double precision, int max_iterations): func(func), a(a), b(b), precision(precision), max_iterations(max_iterations){
        double h = b - a;
        vector<double> R0, R1;
        R0.push_back((func(a)+func(b))*h/2);
        bool done = false;
        for(int k = 1; k < max_iterations; ++k){
            R0.resize(k); R1.resize(k+1);
            R1[0] = (R0[0] + (h/pow(2, k-1))*get_sum(k))/2;
            for(int j = 1; j <= k; ++j){
                R1[j] = R1[j-1] + (R1[j-1]-R0[j-1])/(pow(4, j)-1);
            }
            if(abs(R0.back()-R1.back()) < precision){
                reached_precision = true;
                result = R1.back();
                return;
            }
            R0.swap(R1);
        }
        result = R0.back();
    };
    double result;
    bool reached_precision = false;
private:
    function<double(double)> func;
    double a, b;
    double precision;
    int max_iterations;
    double get_sum(int k);
};

double a_x(double t){
    return sin(t)/(sqrt(t)+1);
}

double a_y(double t){
    return log(t+1)/(t+1);
}

class Solve{
public:
    Solve(double precision, int max_iterations, bool inspect): precision(precision), max_iterations(max_iterations), inspect(inspect){
        init_T();
        if(!inspect){
            for(double t : T){
                cout << x(t) << y(t) << "\n";
            }
        } else {
            for(double t : T){
                double x_ = x(t), y_ = y(t);
            }
            cout << "Ratio with max_iterations set as " << max_iterations << ": " << static_cast<double>(count_reach_precision)/count_total << endl;
        }
    }
private:
    vector<double> T;
    void init_T(){
        for(double i=0.1; i<=10; i+=0.1){
            T.push_back(i);
        }
    }
    double v_x(double t) {
        Romberg romberg(a_x, 0, t, precision, max_iterations);
        if(inspect){
            count_total++;
            if(romberg.reached_precision)   count_reach_precision++;
        }
        return romberg.result;
    }

    double v_y(double t) {
        Romberg romberg(a_y, 0, t, precision, max_iterations);
        if(inspect){
            count_total++;
            if(romberg.reached_precision)   count_reach_precision++;
        }
        return romberg.result;
    }

    double x(double t) {
        // 积分速度 v_x 得到位置 x
        Romberg romberg([this](double tau) { return this->v_x(tau); }, 0, t, precision, max_iterations);
        if(inspect){
            count_total++;
            if(romberg.reached_precision)   count_reach_precision++;
        }
        return romberg.result;
    }

    double y(double t) {
        // 积分速度 v_y 得到位置 y
        Romberg romberg([this](double tau) { return this->v_y(tau); }, 0, t, precision, max_iterations);
        if(inspect){
            count_total++;
            if(romberg.reached_precision)   count_reach_precision++;
        }
        return romberg.result;
    }

    double precision;
    int max_iterations;
    bool inspect;
    int count_total=0, count_reach_precision=0;
};

int main(){
    Solve(1e-6, 8, false);
    Solve(1e-6, 4, true);
    Solve(1e-6, 8, true);
    Solve(1e-6, 12, true);
    Solve(1e-6, 16, true);
    Solve(1e-6, 20, true);
    return 0;
}

double Romberg::get_sum(int k){
    double sum = 0.0;
    for(int i=1; i<=pow(2, k-2); ++i){
        sum += func(a+(2*i-1)*(b-a)/pow(2, k-1));
    }
    return sum;
}