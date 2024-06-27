#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

/*
class ThomasSolve {
public:
    ThomasSolve(const vector<vector<double>>& A, const vector<double>& f): A(A), f(f) {
        int n = A.size();
        result.resize(n);
        vector<double> c(n, 0.0), u(n), y(n);

        u[0] = A[0][0];
        y[0] = f[0] / u[0];

        // 追
        for (int i = 1; i < n; ++i){
            c[i] = A[i][i - 1]; // 下对角线
            double b = (i < n - 1) ? A[i][i + 1] : 0.0; // 上对角线
            u[i] = A[i][i] - c[i] * (A[i - 1][i] / u[i - 1]); // 主对角线减去调整项
            y[i] = (f[i] - c[i] * y[i - 1]) / u[i];
        }

        // 赶
        result[n - 1] = y[n - 1];
        for (int i = n - 2; i >= 0; --i){
            result[i] = y[i] - (A[i][i + 1] / u[i]) * result[i + 1];
        }
    }

    vector<double> result;
private:
    vector<vector<double>> A;
    vector<double> f;
};
*/
class ThomasSolve{
public:
    ThomasSolve(const vector<vector<double> > &A, const vector<double> &f): A(A), f(f){
        int n = A.size();
        result.resize(n);
        vector<double> a(n), b(n), c(n), u(n), v(n), y(n);

        // Initialize vectors
        for(int i = 0; i < n; ++i){
            a[i] = A[i][i];
            if(i > 0){
                c[i] = A[i][i-1];
            }
            if(i < n){
                b[i] = A[i][i+1];
            }
            c[0] = 0, b[n-1] = 0;
        }
        cout << '\n';

        // Solve
        for(int i = 0; i < n; ++i){
            if(i > 0){
                u[i] = a[i] - c[i]*v[i-1];
                y[i] = (f[i]-c[i]*y[i-1])/u[i];
            }
            else{
                u[i] = a[i];
                y[i] = f[i]/u[i];
            }
            v[i] = b[i]/u[i];
        }
        for(int i = n-1; i >= 0; --i){
            if(i < n-1){
                result[i] = y[i] - v[i]*result[i+1];
            }
            else{
                result[i] = y[i];
            }
        }
        for(auto i : u){
            cout << i << " ";
        }
    }
    vector<double> result;
private:
    vector<vector<double> > A;
    vector<double> f;
};

int main(){
    vector<vector<double> > A(5, vector<double>(5, 0.0));
    vector<double> f(5, 3);
    for(int i = 0; i < 5; ++i){
        if(i==0){
            A[i][i+1]=1;
        }
        else if(i==4){
            A[i][i-1]=1;
        }
        else{
            A[i][i-1]=1;
            A[i][i+1]=1;
        }
        A[i][i]=1;
    }
    f[0]=2; f[4]=2;
    ThomasSolve TS(A, f);
    for(auto res : TS.result){
        cout << res << endl;
    }
    return 0;
}