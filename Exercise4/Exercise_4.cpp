#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
const int size=5;
struct Matrix{
    double a[size][size];
};
struct Vector{
    double a[size];
};
void CompleteMatrix(Matrix &matrix) {
    matrix.a[0][0] = 2;
    matrix.a[0][1] = -1;
    matrix.a[0][4] = 1;

    matrix.a[1][0] = -1;
    matrix.a[1][2] = 1;
    matrix.a[1][1] = 2;

    matrix.a[2][1] = 1;
    matrix.a[2][3] = 1;
    matrix.a[2][2] = 1;

    matrix.a[3][2] = 1;
    matrix.a[3][4] = -1;
    matrix.a[3][3] = 2;

    matrix.a[4][4] = 2;
    matrix.a[4][3] = -1;
    matrix.a[4][0] = 1;
}
Matrix SubMatrix(Matrix &matrix,Matrix &matrix1){
    Matrix solve{};
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            solve.a[i][j]=matrix.a[i][j]-matrix1.a[i][j];
        }
    }
    return solve;
}
double Euclidean_Norm(Vector &vector){
    double Norm=0;
    for(double i : vector.a)
        Norm+=i*i;
    return sqrt(Norm);
}
Vector ConstantXVector(double constant,Vector &vector1){
    Vector vector{};
    for(int i=0;i<size;i++)
        vector.a[i]=vector1.a[i]*constant;
    return vector;
}
void PrintVector(Vector &vector){
    for(double i : vector.a)
        cout<<i<<endl;
}
void LUFactorization(Matrix &matrix)
{
    int i, j, k;
    for(k = 0; k < size - 1; k++)
    {
        for(i = k + 1; i<size; i++)
            matrix .a[ i ][ k ] /= matrix .a[ k ][ k ];
        for(i = k + 1; i < size; i++)
            for(j = k + 1; j < size; j++)
                matrix .a[ i ][ j ] -= matrix .a[ i ][ k ] * matrix .a[ k ][ j ];
    }
}
void LUsolve(Matrix &matrix, Vector &vector,Vector &a)
{
    int i, j;
    double s;
    a.a[0]=vector.a[0];
    for(i = 1; i < size; i++){
        s = 0;
        for(j = 0; j < i; j++){
            s+= matrix.a[i][j] * a.a[j];
        }
        a.a[i]= vector.a[i] - s;
    }
    a.a[size - 1]/=matrix.a[size - 1][size - 1];
    for(i = size - 2; i >= 0; i--){
        s = 0;
        for(j = i + 1; j < size; j++) s += matrix.a[i][j] * a.a[ j ];
        a.a[i]= (a.a[i]-s) / matrix.a[i][i];
    }
}
int main(){
    cout<<setprecision(10)<<fixed;

    Matrix A{};
    CompleteMatrix(A);

    double Tau=0.38197;//Numeryczne przybliżenie wartości własnej
    Matrix TauI{};
    for(int i=0;i<size;i++)TauI.a[i][i]=Tau;

    Vector Y{};
    for(int i=0;i<size;i++)
        Y.a[i]=1;
    Y=ConstantXVector((1.0/Euclidean_Norm(Y)),Y);//Wektor Y unormowany

    Vector Z{};

    Matrix temp{};
    temp= SubMatrix(A, TauI);
    LUFactorization(temp);

    for(int i=0;i<1000;i++){ //Algorytm
        LUsolve(temp,Y,Z);
        Y = ConstantXVector((1 / Euclidean_Norm(Z)), Z);
    }

    PrintVector(Y);
    return 0;
}
