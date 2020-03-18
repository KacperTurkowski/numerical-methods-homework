#include <iostream>
const int size = 8;
const double eps = 1e-12;
using namespace std;
class Vector{
public:
    double X[size]{};
    double Y[size]{};
    Vector(){
        X[0]=0.062500;
        Y[0]=0.687959;
        X[1]=0.187500;
        Y[1]=0.073443;
        X[2]=0.312500;
        Y[2]=-0.517558;
        X[3]=0.437500;
        Y[3]=-1.077264;
        X[4]=0.562500;
        Y[4]=-1.600455;
        X[5]=0.687500;
        Y[5]=-2.080815;
        X[6]=0.812500;
        Y[6]=-2.507266;
        X[7]=0.937500;
        Y[7]=-2.860307;
    }
};
class coefficients{
public:
    double a[size]{};
    coefficients(){
        for(double & i : a)
            i=0;
    }
};
class Matrix{
public:
    double a[size][size]{};
    explicit Matrix(Vector &vector){
        for(int j=size-1;j>=0;j--)
            for(int i=size-1;i>=0;i--){
                if(j==size-1) a[i][j]=1;
                else a[i][j]= a[i][j+1] * vector.X[i];
            }
    }
    Matrix(){
        for(int j=size-1;j>=0;j--){
            for(int i=0;i<size;i++){
                a[j][i]=0;
            }
        }
    }
};
void LUfactorization(Matrix &matrix){
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
void LUsolve(Matrix &matrix, Vector &vector, coefficients &a)
{
    int i, j;
    double s;
    a.a[0]=vector.Y[0];
    for(i = 1; i < size; i++){
        s = 0;
        for(j = 0; j < i; j++){
            s+= matrix.a[i][j] * a.a[j];
        }
        a.a[i]= vector.Y[i] - s;
    }
    a.a[size - 1]/=matrix.a[size - 1][size - 1];
    for(i = size - 2; i >= 0; i--){
        s = 0;
        for(j = i + 1; j < size; j++) s += matrix.a[i][j] * a.a[ j ];
        a.a[i]= (a.a[i]-s) / matrix.a[i][i];
    }
}
int main() {
    Vector vector{};
    coefficients a{};
    Matrix matrix{vector};
    LUfactorization(matrix);
    LUsolve(matrix,vector,a);
    for(int i=0;i<size; i++){
        if(size-1-i==0) cout<<"+"<<a.a[i];
        else{
            if(a.a[i]>0)    cout<<"+"<<a.a[i]<<"*x^"<<size-1-i<<" ";
            else    cout<<a.a[i]<<"*x^"<<size-1-i<<" ";
        }
    }
    return 0;
}