#include <iostream>

using namespace std;

const int size=8;
struct Matrix{
    double a[size][size];
};
void PrintVector(double *vector){
    for(int i=1;i<=size-1;i++) cout << vector[i] << endl;
    cout<<endl;
}

void Thomas(Matrix matrix){
    double beta[size];
    double gamma[size];
    double VectorX[size];
    beta[1]= -(matrix.a[1][2]) / matrix.a[1][1];
    gamma[1]= 1 / matrix.a[1][1];
    for (int i = 2; i <= size-1; i++) {
        if(i==size-1) beta[i]=0;
        else beta[i]= -(matrix.a[i][i + 1]) / (matrix.a[i][i - 1] * beta[i - 1] + matrix.a[i][i]);
        gamma[i]= (i - matrix.a[i][i - 1] * gamma[i - 1]) / (matrix.a[i][i - 1] * beta[i - 1] + matrix.a[i][i]);
    }
    for(int n=size-1;n>=1;n--){
        if(n==7) VectorX[n]=gamma[n];
        else VectorX[n]= beta[n] * VectorX[n + 1] + gamma[n];
    }
    PrintVector(VectorX);
}

int main() {
    cout.setf( ios::fixed );
    Matrix matrix{};
    for(int i=1;i<=size-1;i++){
        matrix.a[i][i]=4;
        matrix.a[i - 1][i]=1;
        matrix.a[i][i - 1]=1;
    }
    Thomas(matrix);
}