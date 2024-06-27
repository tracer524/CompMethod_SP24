#include <iostream>

using namespace std;

void PrintMatrix(vector<vector<double> > matrix){
    int row_num = matrix.size();
    int col_num = matrix[0].size();
    cout << "--------------------------------" << endl;
    for(int i = 0; i < row_num; i++){
        for(int j = 0; j < col_num; j++){
            cout << matrix[i][j] << '\t';
        }
        cout << '\n';
    }
    cout << "--------------------------------" << endl;
}

void PrintVector(vector<double> vector){
    cout << "--------------------------------" << endl;
    for(int i = 0; i < vector.size(); i++){
        cout << vector[i] << '\t';
    }
    cout << endl;
    cout << "--------------------------------" << endl;
}

class Doolittle{
public:
    Doolittle(const vector<vector<double> > A): A(A){
        int n = A.size();
        L.resize(n, vector<double>(n, 0));
        U.resize(n, vector<double>(n, 0));
        for(int k=0; k<n; k++){
            // Compute U[k][]
            for(int j=k; j<n; j++){
                U[k][j] = A[k][j];
                for(int r=0; r<k; r++)
                    U[k][j] -= L[k][r]*U[r][j];
            }

            // Compute L[][k]
            L[k][k] = 1;
            for(int i=k+1; i<n; i++){
                L[i][k] = A[i][k];
                for(int r=0; r<k; r++)
                    L[i][k] -= L[i][r]*U[r][k];
                L[i][k] /= U[k][k];
            }
        }
    }
    vector<vector<double> > L;
    vector<vector<double> > U;
private:
    vector<vector<double> > A;
};

class LUSolve{
public:
    LUSolve(const vector<vector<double> > L, const vector<vector<double> > U, const vector<double> b): L(L), U(U), b(b){
        int n = L.size();
        solution.resize(n);
        vector<double> temp(n);

        // Solve LY=b
        for(int i = 0; i < n; i++){
            temp[i] = b[i];
            for(int j=0; j<i; j++)
                temp[i] -= L[i][j]*temp[j];
        }

        // Solve UX=Y
        for(int i=n-1; i >= 0; i--){
            solution[i] = temp[i];
            for(int j=i+1; j<n; j++)
                solution[i] -= U[i][j]*solution[j];
            solution[i] /= U[i][i];
        }
    }
    vector<double> solution;
private:
    vector<vector<double> > L;
    vector<vector<double> > U;
    vector<double> b;
};

class ReversedPowerMethod{
public:
    ReversedPowerMethod(const vector<vector<double> > A, const double tolerance) : A(A), tolerance(tolerance){
        int n = A.size();
        Doolittle decomposion(A);
        vector<double> X_temp(n, 1.0);
        X.push_back(X_temp);
        vector<double> Y_temp(n, 1.0);
        Eigenvalues.push_back(1.0);
        double error = 1.0;
        int i;
        for(i = 0; error > tolerance; i++){

            // Get Y
            for(int j=0; j<n; j++){
                Y_temp[j] = X[i][j]*Eigenvalues[i];
            }
            Y.push_back(Y_temp);
            cout << "Y^(" << i << "):" << endl;
            PrintVector(Y_temp);

            // Get next X
            LUSolve solve(decomposion.L, decomposion.U, Y_temp);
            X_temp = solve.solution;
            X.push_back(X_temp);
            cout << "X^(" << i+1 << "):" << endl;
            PrintVector(X_temp);

            // Find the eigenvalue
            double eigenvalue = 0.0;
            for(double &x : X_temp){
                if(fabs(x)>eigenvalue)
                    eigenvalue = fabs(x);
            }
            Eigenvalues.push_back(1.0/eigenvalue);
            cout << "eigenvalue at iteration " << i+1 << ": " << 1.0/eigenvalue << endl;

            // Error evaluation
            error = fabs(Eigenvalues[i]-Eigenvalues[i-1]);
        }
        Eigenvector = Y[i-1];
        iterations = i-1;
    }
    vector<vector<double> > X;
    vector<vector<double> > Y;
    vector<double> Eigenvalues;
    vector<double> Eigenvector;
    int iterations;
private:
    vector<vector<double> > A;
    double tolerance;
};

int main(){
    vector<vector<double> > A1(5, vector<double>(5));
    vector<vector<double> > A2 = {{4, -1, 1, 3}, {16, -2, -2, 5}, {16, -3, -1, 7}, {6, -4, 2, 9}};
    for(int i=0; i<5; i++)
        for(int j=0; j<5; j++)
            A1[i][j] = 1.0/(9-i-j);
    ReversedPowerMethod sol_1(A1, 0.00001);
    ReversedPowerMethod sol_2(A2, 0.00001);
    return 0;
}
