#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <chrono>

using namespace std;

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
        cout << "Time:" << duration.count() << "microsec" << endl;
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

class GaussSeidel{
public:
    GaussSeidel(const vector<vector<double> > A, const vector<double> b, const double tolerance): A(A), b(b), tolerance(tolerance){
        cout << "Gauss-Seidel" << " ";
        auto start = chrono::high_resolution_clock::now();
        if(A.size() != b.size())
            throw runtime_error("Size mismatch.");
        result.resize(A.size(), 0.0);
        Solve();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
        cout << "Time:" << duration.count() << "microsec" << endl;
    }
    vector<double> result;
private:
    vector<vector<double> > A;       // 系数矩阵
    vector<double> b;
    double tolerance;
    void Solve() {
        int n = A.size();
        vector<double> prevResult(n, 0.0);
        bool continueIteration = true;

        while (continueIteration) {
            continueIteration = false;

            for (int i = 0; i < n; ++i) {
                double sum = b[i];
                for (int j = 0; j < n; ++j) {
                    if (i != j) {
                        sum -= A[i][j] * result[j];
                    }
                }
                double newValue = sum / A[i][i];

                if (fabs(newValue - result[i]) > tolerance) {
                    continueIteration = true;
                }

                result[i] = newValue;
            }

            double norm = 0.0;
            for (int i = 0; i < n; ++i) {
                norm += pow(result[i] - prevResult[i], 2);
                prevResult[i] = result[i];
            }
            norm = sqrt(norm);

            if (norm < tolerance) {
                break;
            }
        }
    }
};

class Solution{
public:
    Solution(const double epsilon, const double a, const int n) : epsilon(epsilon), a(a), n(n) {
        Init_A();
        Init_exact_answer();
        Init_b();
        GaussianElimination GE(A, b);
        result_GaussianElimination = GE.result;
        GaussSeidel GS(A, b, 0.0000001);
        result_GaussSeidel = GS.result;
        
        for(int i = 0; i < A.size(); i++){
            cout << result_GaussianElimination[i] << '\t' << exact_answer[i] << '\t' << "误差：" << 100*(result_GaussianElimination[i]-exact_answer[i])/exact_answer[i] << "%" << endl;
        }
        cout << "--------------------------------" << endl;
        for(int i = 0; i < A.size(); i++){
            cout << result_GaussSeidel[i] << '\t' << exact_answer[i] << '\t' << "误差：" << 100*(result_GaussianElimination[i]-exact_answer[i])/exact_answer[i] << "%" << endl;
        }
    }
    vector<double> result_GaussianElimination;
    vector<double> result_GaussSeidel;
    vector<double> exact_answer;     // 对于每一个切分值的精确解
private:
    double epsilon, a;
    int n;
    vector<vector<double> > A;       // 系数矩阵
    vector<double> b;

    double ExactFunc(double x){
        // 方程的精确解
        return (1-a)*(1-exp(-x/epsilon))/(1-exp(-1.0/epsilon))+a*x;
    }

    // 初始化函数
    
    void Init_A(){
        // 初始化系数矩阵
        A.resize(n-1, vector<double>(n-1, 0));
        for(int i = 0; i < n-1; i++){
            for(int j = 0; j < n-1; j++){
                if(i == j){
                    A[i][j] = -(2 * epsilon + 1.0/n);
                }
                else if(i == j-1){
                    A[i][j] = epsilon + 1.0/n;
                }
                else if(i == j+1){
                    A[i][j] = epsilon;
                }
            }
        }
    }

    void Init_b(){
        b.resize(n-1, a/(n*n));
        b[b.size()-1] -= epsilon + 1.0/n;
    }

    void Init_exact_answer(){
        // 计算精确解
        exact_answer.resize(n-1);
        for(int i = 1; i < n; i++){
            exact_answer[i-1] = ExactFunc(static_cast<double>(i) / n);
        }
    }
};

int main(){
    cout << "Testcase 1" << endl;
    Solution sol_1(1, 0.5, 100);
    cout << "Testcase 2" << endl;
    Solution sol_2(0.1, 0.5, 100);
    cout << "Testcase 3" << endl;
    Solution sol_3(0.001, 0.5, 100);
    cout << "Testcase 4" << endl;
    Solution sol_4(0.0001, 0.5, 100);
    return 0;
}
