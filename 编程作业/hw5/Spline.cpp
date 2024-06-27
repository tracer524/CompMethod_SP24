#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

class ReadFile {
public:
    ReadFile(const string& filename) {
        ifstream file(filename);
        string line;
        if (!file.is_open()) {
            cerr << "Failed to open file: " << filename << endl;
            return;
        }

        while (getline(file, line)) {
            istringstream iss(line);
            double xValue, yValue;

            if (!(iss >> xValue >> yValue)) {
                cerr << "Failed to parse line: " << line << endl;
                continue; // Skip bad lines
            }

            x.push_back(xValue);
            y.push_back(yValue);
        }
        file.close();
    }
    vector<double> x;
    vector<double> y;
};

class GaussianElimination {
public:
    GaussianElimination(const vector<vector<double> > A, const vector<double> b) : A(A), b(b) {
        cout << "Gaussian Elimination" << " ";
        auto start = chrono::high_resolution_clock::now();
        if (A.size() != b.size())
            throw runtime_error("Size mismatch.");
        result.resize(A.size());
        Solve();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
        cout << "Time:" << duration.count() << "ms" << endl;
    }
    vector<double> result;

private:
    vector<vector<double> > A;
    vector<double> b;

    int SelectMax(int i) {
        int maxIndex = i;
        double maxVal = fabs(A[i][i]);
        for (int k = i + 1; k < A.size(); k++) {
            if (fabs(A[k][i]) > maxVal) {
                maxVal = fabs(A[k][i]);
                maxIndex = k;
            }
        }
        return maxIndex;
    }

    void Solve() {
        int n = A.size();
        for (int i = 0; i < n; i++) {
            int maxRow = SelectMax(i);
            if (i != maxRow) {
                swap(A[i], A[maxRow]);
                swap(b[i], b[maxRow]);
            }

            for (int k = i + 1; k < n; k++) {
                double c = A[k][i] / A[i][i];
                for (int j = i; j < n; j++) {
                    A[k][j] -= c * A[i][j];
                }
                b[k] -= c * b[i];
            }
        }

        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * result[j];
            }
            result[i] = (b[i] - sum) / A[i][i];
        }
    }
};


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

class ConstructMatrix{
public:
    ConstructMatrix(const vector<double> lambda, const vector<double> mu, const vector<double> d): lambda(lambda), mu(mu), d(d){
        A.resize(lambda.size()-1, vector<double>(lambda.size()-1, 0.0));
        f.resize(lambda.size()-1);
        
        // Initialize A
        for(int i = 0; i < A.size(); ++i){
            if(i == 0){
                A[0][1] = lambda[1];
            }
            else if(i == A.size()-1){
                A[i][i-1] = mu[i+1];
            }
            else{
                A[i][i-1] = mu[i+1];
                A[i][i+1] = lambda[i+1];
            }
            A[i][i] = 2;
        }

        // Initialize f
        for(int i = 0; i < A.size(); ++i){
            f[i] = d[i+1];
        }
    }
    vector<vector<double> > A;
    vector<double> f;
private:
    vector<double> lambda, mu, d;
};

class Solution{
public:
    Solution(const string& filename){
        ReadFile RF(filename);
        int num_points = RF.x.size();
        vector<double> h(num_points-1), lambda(num_points-1), mu(num_points-1), d(num_points-1);
        Initialize(h, lambda, d, mu, RF);
        ConstructMatrix CM(lambda, mu, d);
        ThomasSolve TS(CM.A, CM.f);
        vector<double> M(num_points, 0.0);
        for(int i = 1; i < num_points-1; ++i){
            M[i] = TS.result[i-1];
        }
        vector<vector<double> > coef(num_points-1, vector<double>(4));
        for(int i = 0; i < num_points-1; ++i){
            cout << "第" << i+1 << "个多项式的系数：";
            coef[i][0] = (M[i+1]-M[i])/(6*h[i]);
            cout << coef[i][0] << " ";
            coef[i][1] = (RF.x[i+1]*M[i]-RF.x[i]*M[i+1]) / (2*h[i]);
            cout << coef[i][1] << " ";
            coef[i][2] = (-pow(RF.x[i+1], 2)*M[i]+pow(RF.x[i], 2)*M[i+1]) / (2*h[i]) + (RF.y[i+1]-RF.y[i]) / h[i] - h[i]*(M[i+1]-M[i]) / 6;
            cout << coef[i][2] << " ";
            coef[i][3] = (pow(RF.x[i+1], 3)*M[i]-pow(RF.x[i], 3)*M[i+1]) / (6*h[i]) + (RF.x[i+1]*RF.y[i]-RF.x[i]*RF.y[i+1]) / h[i] - h[i]*(RF.x[i+1]*M[i]-RF.x[i]*M[i+1]) / 6;
            cout << coef[i][3] << " " << endl;
        }
}
private:
    void Initialize(vector<double> &h, vector<double> &lambda, vector<double> &d, vector<double> &mu, ReadFile RF) {
        for(int i = 0; i < RF.x.size()-1; i++){
            h[i] = RF.x[i+1] - RF.x[i];
        }
        for(int i = 1; i < RF.x.size()-1; i++){
            lambda[i] = h[i]/(h[i]+h[i-1]);
            mu[i] = 1-lambda[i];
            d[i] = 6*((RF.y[i+1]-RF.y[i])/h[i]-(RF.y[i]-RF.y[i-1])/h[i-1])/(h[i]+h[i-1]);
        }
    }
};

int main(){
    Solution sol1("point.txt");
    Solution sol2("point2.txt");
    return 0;
}