#include <iostream>

using namespace std;

const int size=8;
struct Matrix{
    double diagonal[size];
    double underdiagonal[size];
};
struct Vector{
    double a[size];
};
void print(Vector &x){
    for(int i=1;i<=size-1;i++){
        cout<<x.a[i]<<endl;
    }
    cout<<endl;
}
void Sherman_Morrison(Matrix &A){
    Vector Result{};
    for(int i=0;i<=size-1;i++)
        Result.a[i]=i;
    Vector beta{};
    beta.a[1]= -A.underdiagonal[1] / A.diagonal[1];
    for(int i=2;i<=size-1;i++){
        if(i==size-1) beta.a[i]=0;
        else beta.a[i]= -(A.underdiagonal[i]) / (A.underdiagonal[i] * beta.a[i - 1] + A.diagonal[i]);
    }
    Vector gamma{};
    gamma.a[1]=1/A.diagonal[1];
    for(int i=2;i<=size-1;i++){
        gamma.a[i]= (i-(A.underdiagonal[i] * gamma.a[i - 1])) / (A.underdiagonal[i] * beta.a[i - 1] + A.diagonal[i]);
    }
    Vector Z{};
    Z.a[size-1]=gamma.a[size-1];
    for(int i=size-2;i>=1;i--){
        Z.a[i]=beta.a[i]*Z.a[i+1]+gamma.a[i];
    }
    gamma.a[1]=1/A.diagonal[1];
    for(int i=2;i<=size-2;i++){
        gamma.a[i]= (0-(A.underdiagonal[i] * gamma.a[i - 1])) / (A.underdiagonal[i] * beta.a[i - 1] + A.diagonal[i]);
    }
    gamma.a[size-1]= (1-(A.underdiagonal[size-1] * gamma.a[size-2])) / ((A.underdiagonal[size-1] * beta.a[size-2] + A.diagonal[size-1]));
    Vector Q{};
    Q.a[size-1]=gamma.a[size-1];
    for(int i=size-2;i>=1;i--){
        Q.a[i]=beta.a[i]*Q.a[i+1]+gamma.a[i];
    }
    Vector X{};
    double numerator= Z.a[1] + Z.a[size - 1];
    double denominator= Q.a[1] + Q.a[size - 1] + 1;
    for(int i=1;i<=size-1;i++){
        X.a[i]=Z.a[i]- (numerator / denominator) * Q.a[i];
    }
    print(X);
}
int main() {
    Matrix A{};
    for(int i=1;i<=size-1;i++){
        if(i==7||i==1)A.diagonal[i]=3;
        else A.diagonal[i]=4;
        A.underdiagonal[i]=1;
    }
    Sherman_Morrison(A);
    return 0;
}