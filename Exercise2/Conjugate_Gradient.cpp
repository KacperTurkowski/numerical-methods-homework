#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
const int size = 128;
const double eps=0.00001;
class Vector{
public:
    double a[size]{};
    Vector(){
        for(double & i : a){
            i=0;
        }
    }
};
Vector MatrixXVector(Vector &vector){
    Vector Result{};
    for(int i=0; i < size; i++){
        if(i==0)Result.a[i]=(4*vector.a[i])+vector.a[i+1]+vector.a[i+4];
        else if(i>0&&i<=3)Result.a[i]=vector.a[i-1]+(4*vector.a[i])+vector.a[i+1]+vector.a[i+4];
      	else if(i>3&&i<=size-5)Result.a[i]=vector.a[i+4]+vector.a[i+1]+(4*vector.a[i])+vector.a[i-1]+ vector.a[i-4];
        else if(i>size-5&&i<=size-2)Result.a[i]=vector.a[i+1]+(4*vector.a[i])+vector.a[i-1]+vector.a[i-4];
        else if(i > size - 2) Result.a[i]= (4 * vector.a[i]) + vector.a[i - 1] + vector.a[i - 4];
    }
    return Result;
}
double vectorXvectorT(Vector &normalvector, Vector &transposedvector){
    double sum=0;
    for(int i=0; i < size; i++){
        sum+= normalvector.a[i] * transposedvector.a[i];
    }
    return sum;
}
Vector constantXvector(Vector &vector, double constant){
    Vector result{};
    for(int i=0; i < size; i++){
        result.a[i]= constant * vector.a[i];
    }
    return result;
}
Vector subvector(Vector &vector1,Vector &vector2){
    Vector result{};
        for(int i=0; i < size; i++){
            result.a[i]= vector1.a[i] - vector2.a[i];
        }
        return result;
}
Vector sumvector(Vector &vector1, Vector &vector2){
    Vector result{};
    for(int i=0; i < size; i++){
        result.a[i]= vector1.a[i] + vector2.a[i];
    }
    return result;
}
double euklidesnorm(Vector &vector){
    double norma=0;
    for(double i : vector.a)norma+= i * i;
    return sqrt(norma);
}
void PrintVector(Vector &vector){
    for(int i=0;i<128;i++) cout<<"i: "<<i<<" "<<vector.a[i]<<endl;
}
void GS(){
    Vector R[2]{};
    Vector b{};
    Vector X[2]{};
    Vector P[2]{};
    for(double & i : X[0].a)i=1;
    Vector temp{};
    for(double & j : b.a) j=1;
    temp= MatrixXVector(X[0]);
    R[0]= subvector(b, temp);
    P[0]=R[0];
    double alfa;
    double beta;
    while(euklidesnorm(R[0]) > eps){
        temp= MatrixXVector(P[0]);//A * P
        alfa= (vectorXvectorT(R[0], R[0])) / (vectorXvectorT(P[0], temp));
        temp= constantXvector(temp, alfa);//alfa *A *p
        R[1]= subvector(R[0], temp);//r - alfa * A * p
        beta= vectorXvectorT(R[1], R[1]) / vectorXvectorT(R[0], R[0]);
        temp= constantXvector(P[0], beta);
        P[1]= sumvector(R[1], temp);
        temp= constantXvector(P[0], alfa);
        X[1]= sumvector(X[0], temp);
        R[0]=R[1];
        P[0]=P[1];
        X[0]=X[1];
    }
    PrintVector(X[1]);
}
int main(){
    cout<<setprecision(10)<<fixed;
    GS();
    return 0;
}
