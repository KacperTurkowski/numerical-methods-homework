#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
#include <iomanip>

const double eps = 1e-15;
typedef double type;
const int size =6;
const double divider=12;
using namespace std;
vector<type> sub_vector(vector<type>& vector1, vector<type> vector2) {
    vector<type> temp(size);
    for(int i = 0; i < size; i++) temp[i] = vector1[i] - vector2[i];
    return temp;
}
void vectorxconstant(vector<type>& vector1, type constant) {
    for(int i = 0; i < size; i++) vector1[i] *= constant;
}
type Dot_product(vector<type> vector1, vector<type> vector2) {
    double result = 0;
    for(int i = 0; i < size; i++)result += vector1[i] * vector2[i];
    return result;
}
type Norm(vector<type> &vector1) {
    double sum = 0;
    for(double i : vector1)sum += pow(i, 2);
    return sqrt(sum);
}
vector<type> matrixXvector(vector<vector<type> >& matrix, vector<type>& vector1) {
    vector<type> bufor(size);
    for(int i = 0; i < size; i++) {
        bufor[i] = 0;
        for(int j = 0; j < size; j++) bufor[i] += matrix[i][j] * vector1[j];
    }
    return bufor;
}
vector<type> vectorXconstant(vector<type>& vector1, type constant) {
    vector<type> temp(size);
    for(int i = 0; i < size; i++)temp[i] = vector1[i] * constant;
    return temp;
}
void vector_normalization(vector<type>& vector1) {
    vectorxconstant(vector1, 1 / Norm(vector1));
}
vector<type> orthogonal_vector(vector<type> vector1) {
    vector < type > temp(size);
    double sum = 0;
    for(int i = 0; i < size-1; i++) {
        temp[i] = rand() % 100 + 1;
        sum += temp[i] * vector1[i];
    }
    temp[size - 1] = -sum / vector1[size - 1];
    return temp;
}
vector<type> Random_vector() {
    vector<type> vector;
    vector.reserve(size);
    for(int i = 0; i <size; i++) vector.push_back(rand()%10+1);
    vector_normalization(vector);
    return vector;
}
void Print_Vector(vector<type>& vector) {
    cout<<"[ ";
    for(double i:vector) cout<<i<<", ";
    cout<<"]"<<endl;
}
void Print(type lambda,vector<type>& vector){
    cout<<setprecision(1)<<fixed;
    cout << "Eigenvalue: " << lambda / divider << endl;
    cout<<setprecision(14)<<fixed;
    cout<<"Eigenvector: "<<endl;
    Print_Vector(vector);
}
void Power_method(vector<vector<type> >& matrix, vector<type> Y ) {
    vector<type> z,v1, v2, sub, ZN;
    type Lambda1, Lambda2;
    do {
        z = matrixXvector(matrix, Y);
        ZN = z;
        vector_normalization(z);
        sub = sub_vector(Y, z);
        Y = z;
    } while(Norm(sub) >= eps);
    Lambda1 = Norm(ZN);
    v1 = z;
    Print(Lambda1,v1);
    Y = orthogonal_vector(v1);
    do {
        z = matrixXvector(matrix, Y);
        z = sub_vector(z, vectorXconstant(v1, Dot_product(v1, z)));
        ZN = z;
        vector_normalization(z);
        sub = sub_vector(Y, z);
        Y = z;
    } while(Norm(sub) >= eps);
    Lambda2 = Norm(ZN);
    v2 = z;
    Print(Lambda2,v2);
}
int main() {
    srand(time(nullptr));
    vector<vector<type>> A = {{19,  13,  10, 10, 13,  -17},
                              {13,  13,  10, 10, -11, 13},
                              {10,  10,  10, -2, 10,  10},
                              {10,  10,  -2, 10, 10,  10},
                              {13,  -11, 10, 10, 13,  13},
                              {-17, 13,  10, 10, 13,  19} };
    Power_method(A, Random_vector());
    return 0;
}