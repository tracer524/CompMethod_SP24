#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <cmath>

#define PI M_PI

using namespace std;

class FFT{
public:
    FFT(const vector<complex<double> >& f): f(f){
        int n = f.size();
        g.resize(n);
        if(n==1){
            g = f;
            return;
        }
        complex<double> w_n = exp(-2*PI/n*complex<double>(0.0, 1.0));
        complex<double> w(1, 0);
        vector<complex<double> > f0, f1;
        for(int i = 0; i < n; ++i){
            if(i%2){
                f1.push_back(f[i]);
            }
            else{
                f0.push_back(f[i]);
            }
        }
        vector<complex<double> > g0 = FFT(f0).g, g1 = FFT(f1).g;
        for(int i=0; i<n/2; ++i){
            g[i] = (g0[i]+w*g1[i])/complex<double>(2, 0);
            g[i+n/2] = (g0[i]-w*g1[i])/complex<double>(2, 0);
            w *= w_n;
        }
    }
    vector<complex<double> > g;
private:
    vector<complex<double> > f;
};

class IFFT{
public:
    IFFT(const vector<complex<double> >& f): f(f){
        int n = f.size();
        g.resize(n);
        if(n==1){
            g = f;
            return;
        }
        complex<double> w_n = exp(2*PI/n*complex<double>(0.0, 1.0));
        complex<double> w(1, 0);
        vector<complex<double> > f0, f1;
        for(int i = 0; i < n; ++i){
            if(i%2){
                f1.push_back(f[i]);
            }
            else{
                f0.push_back(f[i]);
            }
        }
        vector<complex<double> > g0 = IFFT(f0).g, g1 = IFFT(f1).g;
        for(int i=0; i<n/2; ++i){
            g[i] = g0[i]+w*g1[i];
            g[i+n/2] = g0[i]-w*g1[i];
            w *= w_n;
        }
    }
    vector<complex<double> > g;
private:
    vector<complex<double> > f;
};

class Sample{
public:
    Sample(const int& exp, const bool& noise=false): exp(exp), noise(noise){
        int n = pow(2, exp);
        for(int i=0; i<n; ++i){
            double y = function(static_cast<double>(i) / n);
            if(noise){
                y += 0.3*random_number();
            }
            f.push_back(complex<double>(y, 0.0));
        }
    }
    vector<complex<double> > f;
private:
    bool noise;
    int exp;
    double function(double x){
        return 0.7*sin(4*PI*x)+sin(10*PI*x);
    }
    double random_number(){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0, 1);
        return dis(gen);
    }
};

class Solution{
public:
    Solution(const int& expn, const bool& noise, const bool& truncate): expn(expn), noise(noise), truncate(truncate){
        cout << "Under exponent of " << expn << ":" << '\n';
        vector<complex<double> > f = Sample(expn, noise).f;
        cout << "Values of f: " << '\n';
        for(auto i : f){
            cout << i.real() << ", ";
        }
        cout << '\n';
        vector<complex<double> > g = FFT(f).g;
        vector<double> g_abs;
        for(int i = 0; i <g.size(); ++i){
            g_abs.push_back(abs(g[i]));
        }
        cout << "Values of g_abs: " << '\n';
           for(auto i : g_abs){
            cout << i << ", ";
        }
        cout << '\n';
        vector<complex<double> > g_truncated;
        vector<complex<double> > rebuilded;
        if(truncate){
            for(int i=0; i<g.size(); ++i){
                if(i<g.size()*0.25){
                    g_truncated.push_back(g[i]);
                }
                else{
                    g_truncated.push_back(complex<double>(0.0, 0.0));
                }
            }
            rebuilded = IFFT(g_truncated).g;
        }
        else{
            rebuilded = IFFT(g).g;
        }
        cout << "Values of rebuilded vector(real part): " << '\n';
        for(auto i : rebuilded){
            cout << i.real() << ", ";
        }
        cout << '\n';
    }
private:
    bool noise;
    bool truncate;
    int expn;
};

int main(){
    cout << "---------------------------------------" << '\n';
    Solution(4, false, false);
    cout << "---------------------------------------" << '\n';
    Solution(7, false, false);
    cout << "---------------------------------------" << '\n';
    Solution(7, true, false);
    cout << "---------------------------------------" << '\n';
    Solution(7, true, true);
    cout << "---------------------------------------" << '\n';
    return 0;
}