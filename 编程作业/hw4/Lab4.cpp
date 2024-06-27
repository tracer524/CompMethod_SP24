#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

class Matrix{
public:
    int cols, rows;
    vector<vector<double> > data;
    Matrix(){}
    Matrix(vector<vector<double> > &data): data(data){
        rows = data.size();
        cols = data[0].size();
    }
    Matrix T(){
        vector<vector<double> > t(cols, vector<double>(rows));
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                t[j][i] = data[i][j];
            }
        }
        return Matrix(t);
    }
    Matrix dot(Matrix A){
        vector<vector<double> > result(rows, vector<double>(A.cols, 0.0));
        if(cols != A.rows)
            throw runtime_error("Size mismatch.");
        for(int i=0; i<rows; ++i){
            for(int j=0; j<A.cols; ++j){
                for(int k=0; k<cols; ++k){
                    result[i][j] += data[i][k]*A.data[k][j];
                }
            }
        }
        return Matrix(result);
    }
    Matrix Multiply(double x){
        for(int i=0; i<rows; ++i){
            for(int j=0; j<cols; ++j){
                data[i][j] = x*data[i][j];
            }
        }
        return Matrix(data);
    }
    void Print(){
        cout << "--------------------------------" << endl;
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                cout << data[i][j] << ' ';
            }
            cout << '\n';
        }
        cout << "--------------------------------" << endl;
    }
};

class Jacobi{
public:
    vector<double> eigenvalue;
    Matrix eigenvectors;
    vector<int> indices;
    Jacobi(Matrix& M, const double tolerance, bool print): M(M), tolerance(tolerance), backup(M), print(print){
        if(M.rows<M.cols) dim = M.rows;
        else dim = M.cols;
        vector<vector<double> > temp(dim, vector<double>(dim, 0));
        for(int i = 0; i < dim; ++i){
            temp[i][i] = 1;
        }
        eigenvectors = Matrix(temp);
        indices.resize(dim);
    }
    void Iterate(){
        int iterations = 0;
        while(Error() > tolerance && iterations < 10000){
            // Find the maximum value(except for diagonal)
            int p = 0, q = 1;
            double max = 0;
            for(int i=0; i<M.rows; ++i){
                for(int j=1; j<M.cols && i!=j; ++j){
                    if(abs(M.data[i][j]) > max){
                        p = i; q = j; max = abs(M.data[i][j]);
                    }
                }
            }

            // Determine the angle
            double s = (M.data[q][q]-M.data[p][p]) / (2*M.data[p][q]);
            double t;
            if(s == 0)
                t = 1.0;
            else{
                double t1 = -s-sqrt(pow(s, 2)+1), t2 = -s+sqrt(pow(s, 2)+1);
                if(abs(t1)>abs(t2))
                    t = t2;
                else
                    t = t1;
            }
            double c = 1/sqrt(1+pow(t, 2)), d = t*c;
            /*
            for(int i = 0; i<n && i!=p && i!=q; ++i){
                M.data[p][i] = c*M.data[p][i] - d*M.data[q][i];
                M.data[i][p] = M.data[p][i];
                M.data[q][i] = c*M.data[q][i] + d*M.data[p][i];
                M.data[i][q] = M.data[q][i];
            }
            M.data[p][p] -= t*M.data[p][q];
            M.data[q][q] += t*M.data[p][q];
            */
            vector<vector<double> > Q_data(M.rows, vector<double>(M.cols, 0.0));
            for(int i = 0; i < dim; ++i){
                Q_data[i][i] = 1;
            }
            Q_data[p][p] = c;
            Q_data[p][q] = d;
            Q_data[q][q] = c;
            Q_data[q][p] = -d;
            Matrix Q(Q_data);
            M = Q.T().dot(M).dot(Q);
            eigenvectors = eigenvectors.dot(Q);

            if(print){
                cout << "Error after iteration" << iterations+1 << ": " << Error()/2 << endl;
            }
            iterations++;
        }
        eigenvalue.resize(dim);
        for(int i=0; i<dim; ++i){
            eigenvalue[i] = M.data[i][i];
        }

        // Create the indices
        for(int i = 0; i < eigenvalue.size(); ++i){
            indices[i] = i;
        }

        // Sort the eigenvalues in descending order
        sort(indices.begin(), indices.end(), [this](int a, int b){
            return eigenvalue[a] > eigenvalue[b];
        });
        if(print){
            cout << "Eigenvalues:" << endl;
            for(double index : indices){
                cout << eigenvalue[index] << '\t';
            }
            cout << endl;
        }
    }
private:
    Matrix M;
    Matrix backup;
    double tolerance;
    bool print;
    int dim;

    double Error(){
        double result = 0.0;
        for(int i = 0; i < M.rows; ++i)
            for(int j = 0; j<M.cols && i!=j; ++j)
                result += pow(M.data[i][j], 2);
        return result;
    }

};

class SVD_Decomposition{
public:
    SVD_Decomposition(Matrix& A, const double tolerance): A(A), tolerance(tolerance){
        Matrix ATA = A.T().dot(A);
        Jacobi jacobi_ATA(ATA, tolerance, true);
        jacobi_ATA.Iterate();
        vector<vector<double> > V_data(ATA.rows, vector<double>(ATA.rows, 0.0));
        for(int i=0; i<ATA.rows; ++i){
            for(int j=0; j<ATA.rows; ++j){
                V_data[i][j] = jacobi_ATA.eigenvectors.data[j][i];
            }
        }
        V = Matrix(V_data);
        cout << "V:" << endl;
        V.Print();

        Matrix AAT = A.dot(A.T());
        Jacobi jacobi_AAT(AAT, tolerance, false);
        jacobi_AAT.Iterate();
        vector<vector<double> > U_data(AAT.rows, vector<double>(AAT.rows, 0.0));
        for(int i=0; i<AAT.rows; ++i){
            for(int j=0; j<AAT.rows; ++j){
                U_data[i][j] = jacobi_AAT.eigenvectors.data[j][i];
            }
        }
        U = Matrix(U_data);
        cout << "U:" << endl;
        U.Print();

        vector<vector<double> > Sigma_data(A.rows, vector<double>(A.rows, 0.0));
        for(int i=0; i<A.rows; ++i){
            Sigma_data[i][i] = sqrt(jacobi_AAT.eigenvalue[jacobi_AAT.indices[i]]);
        }
        Sigma = Matrix(Sigma_data);
        cout << "Sigma:" << endl;
        Sigma.Print();
    }
private:
    Matrix A;
    Matrix U;
    Matrix V;
    Matrix Sigma;
    double tolerance;
};

class PCA{
public:
    PCA(vector<vector<double> >& data_vectors, const double tolerance, const int(ndim)): data_vectors(data_vectors), tolerance(tolerance), ndim(ndim){
        if(data_vectors[0].size()<ndim){
            throw runtime_error("ndim too large");
        }
        data_vectors = Decentralize(data_vectors);
        coordinates.resize(ndim);
        vector<vector<double> > data = data_vectors;
        for(auto& row : data){
            row.pop_back();
        }
        X = Matrix(data).T();
        Cov = X.dot(X.T()).Multiply(1.0/data.size());

        cout << "Covariance matrix:" << endl;
        Cov.Print();

        Jacobi jacobi_cov(Cov, tolerance, false);
        jacobi_cov.Iterate();
        vector<double> eigVec;
        eigVec.resize(jacobi_cov.eigenvectors.data.size());
        for(int i=0; i<2; ++i){
            int index = jacobi_cov.indices[i];
            for(int j=0; j<jacobi_cov.eigenvectors.data.size(); ++j){
                eigVec[j] = jacobi_cov.eigenvectors.data[j][index];
            }
            base.push_back(eigVec);
        }
        GetCoordinates();
        writeDataToFile("0.txt", coordinates[0]);
        writeDataToFile("1.txt", coordinates[1]);
        writeDataToFile("2.txt", coordinates[2]);
    }
private:
    vector<vector<double> > data_vectors;
    Matrix X;
    Matrix Cov;
    double tolerance;
    int ndim;
    vector<vector<double> > base;
    vector<vector<vector<double> > > coordinates;
    vector<vector<double> > Decentralize(vector<vector<double> > data){
        vector<double> average(data[0].size()-1, 0.0);
        for(int i = 0; i < data.size(); ++i){
            for(int j = 0; j < data[0].size()-1; ++j){
                average[j] += (data[i][j]);
            }
        }
        for(auto& avg : average){
            avg /= data.size();
        }
        for(int i = 0; i < data.size(); ++i){
            for(int j = 0; j < data[0].size()-1; ++j){
                data[i][j] -= average[j];
            }
        }
        return data;
    }
    void GetCoordinates(){
        coordinates.resize(3);
        for(auto data_vector : Decentralize(data_vectors)){
            double x1 = 0;
            double x2 = 0;
            for(int i=0; i<data_vector.size()-1; ++i){
                x1 += data_vector[i]*base[0][i];
                x2 += data_vector[i]*base[1][i];
            }
            vector<double> temp;
            temp.push_back(x1);
            temp.push_back(x2);
            coordinates[data_vector[4]].push_back(temp);
        }
    }
    void writeDataToFile(const string& filename, const vector<vector<double>>& data) {
        ofstream file(filename);
        if (!file) {
            cerr << "Unable to open file for writing: " << filename << endl;
            return;
        }
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i < row.size() - 1) file << ",";
            }
            file << "\n";
        }
        file.close();
    }
};

class Solution1{
public:
    Solution1(const double tolerance): tolerance(tolerance){
        cout << "Question 1:" << endl;
        vector<vector<double> > data = GenerateMatrix(4, 3);
        Matrix A(data);
        cout << "Randomly generated matrix:" << endl;
        A.Print();
        SVD_Decomposition svd(A, tolerance);
    }
private:
    vector<vector<double> > GenerateMatrix(int m, int n){
        vector<vector<double> > result(m, vector<double>(n));
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0.0, 1.0);
        for(int i = 0; i < m; ++i){
            for(int j = 0; j < n; ++j){
                result[i][j] = dis(gen);
            }
        }
        return result;
    }
    double tolerance;
};

class Solution2{
public:
    Solution2(const double tolerance):tolerance(tolerance){
        cout << "Question 2:" << endl;
        vector<vector<double> > data = readData("iris.txt");
        PCA pca(data, tolerance, 2);
    }
private:
    double tolerance;
    vector<vector<double>> readData(const string& filename) {
        vector<vector<double>> data;
        ifstream file(filename);
        string line;
        if (!file) {
            cerr << "Unable to open file: " << filename << endl;
            return data;
        }
        while (getline(file, line)) {
            stringstream ss(line);
            string value;
            vector<double> row;
            while (getline(ss, value, ',')) { 
                row.push_back(stod(value)); 
            }
            data.push_back(row);
        }
        file.close();
        return data;
    }
};

int main(){
    Solution1 solution1(0.000001);
    Solution2 solution2(0.000001);
    return 0;
}
