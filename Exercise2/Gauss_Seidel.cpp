#include <iostream>
#include <iomanip>

using namespace std;
const int size = 129;
struct Vector{
    double x[size];
};
void CompleteVector(Vector &vector){
    for(double & i : vector.x)i=0;
}
void PrintVector(Vector &vector){
    for(int i=1; i < size; i++)
    cout << i << ": " << vector.x[i] << endl;
}
void Gauss_Seidel(Vector &vector){
    for(int j=0;j<=100;j++) {
        for (int i=1; i<size; i++) {
        if (i == 1) { vector.x[i] = (1 - vector.x[i + 1] - vector.x[i + 4]) / 4; }
        else if(i>1&&i<=4){vector.x[i]=(1-vector.x[i-1]-vector.x[i+1]-vector.x[i+4])/4;}
    	else if(i>4&&i<=size-5){vector.x[i]=(1-vector.x[i-4]-vector.x[i-1]-vector.x[i+1]-vector.x[i+4])/4;}
        else if(i>size-5&&i<=size-2){vector.x[i]=(1-vector.x[i-4]-vector.x[i-1]-vector.x[i+1])/4;}
        else if (i == size - 1) { vector.x[i] = (1 - vector.x[i - 4] - vector.x[i - 1]) / 4; }
        }
    }
}
int main(){
    cout<<setprecision(10)<<fixed;
    Vector vector{};
    CompleteVector(vector);
    Gauss_Seidel(vector);
    PrintVector(vector);
    return 0;
}
