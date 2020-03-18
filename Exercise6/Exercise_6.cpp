#include <iostream>
#include <cmath>
using namespace std;
const double e=M_E;
const double pi=M_PI;
const int size = 20;
struct Vector{
    double a[size + 1];
};
int main()
{
    Vector R[size + 1]{};
    double multiplier[size];
    for (int i=1; i < size + 1; i++) {
        multiplier[i] = 17.0 / pow(2, i - 1);
    }
    R[1].a[1]= multiplier[1] / 2 * ((sin(pi * (1 + sqrt(0)) / (1 + 0 * 0))) * pow(e, -0) + (sin(pi * (1 + sqrt(17)) / (1 + 17 * 17))) * pow(e, -17));
    for (int i = 2; i < size + 1; ++i) {
        double coefficient = 0;
        for (int k=1;k<=pow(2,i-2);k++) {
            coefficient+= (sin(pi * (1 + sqrt(((2 * k - 1) * multiplier[i]))) / (1 + ((2 * k - 1) * multiplier[i]) * ((2 * k - 1) * multiplier[i])))) * pow(e, -((2 * k - 1) * multiplier[i]));
        }
        R[i].a[1] = 0.5*(R[i-1].a[1] + multiplier[i - 1] * coefficient);
    }
    for (int i=2; i <= size; i++) {
        for (int j=2;j<=i;j++){
            R[i].a[j]=R[i].a[j-1]+(R[i].a[j-1]-R[i-1].a[j-1])/(pow(4,j-1)-1);
        }
    }
    cout << R[size].a[size] << endl;
    return 0;
}


