#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>
using namespace std;

typedef complex<double> Complex;
typedef vector<Complex> Polynomial;
const double EPS = 0.00000001;
const int Max_iteration =100;
pair<Polynomial, Complex> horner(Polynomial &polynomial, Complex Z) {
    int size = polynomial.size();
    Polynomial temp = Polynomial(max(1, size - 1));
    for(int i= size - 1; i > 0; i--){
        temp[i - 1]=polynomial[i];
        if(i < size - 1) temp[i - 1]+= temp[i] * Z;
    }
    return make_pair(temp, polynomial[0] + temp[0] * Z);
}
void Complete_Polynomial(Polynomial &polynomial){
    polynomial.push_back(1);
    polynomial.push_back(-6);
    polynomial.push_back(3);
    polynomial.push_back(2);
    polynomial.push_back(1);
}
Polynomial derivative(Polynomial &polynomial) {
    int n = polynomial.size();
    Polynomial temp=Polynomial(max(1, n - 1));
    for(int i = 1; i < n; i++)
        temp[i - 1]= polynomial[i] * Complex(i);
    return temp;
}
int Compare(Complex x, Complex y) {
    double sub = abs(x) - abs(y);
    if(sub < -EPS) return -1;
    else{
        if(sub > EPS) return 1;
        else return 0;
    }
}
Complex First_Zero_Place(Polynomial &polynomial, Complex x) {
    int n = (int)polynomial.size() - 1;
    Polynomial temp1 = derivative(polynomial);
    Polynomial temp2 = derivative(temp1);
    for (int iteration = 0; iteration < Max_iteration; iteration++) {
        Complex Y=horner(polynomial, x).second;
        if (Compare(Y, 0) == 0) break;
        Complex G = horner(temp1, x).second / Y;
        Complex sub = sqrt(Complex(n - 1) * ((G * G - horner(temp2, x).second - Y) * Complex(n) - G * G));
        Complex divider1 = G - sub;
        Complex divider2 = G + sub;
        Complex a=Complex(n);
        if(Compare(divider1, divider2) > 0) a/=divider1;
        else a/=divider2;
        x -= a;
        if (Compare(a, 0) == 0) break;
    }
    return x;
}
vector<Complex> Zero_Places(Polynomial &polynomial) {
    vector<Complex> result;
    Polynomial temp = polynomial;

    while (temp.size() > 2) {
        Complex z(0.0, 0.0);
        z = First_Zero_Place(temp, z);
        temp = horner(temp, z).first;
        result.push_back(z);
    }
    result.push_back(-temp[0] / temp[1]);
    return result;
}

int main() {
    Polynomial p;
    Complete_Polynomial(p);
    vector<Complex> result = Zero_Places(p);
    cout<<setprecision(10)<<fixed;
    for(int i=0; i < result.size(); i++) {
        if(result[i].imag() > 0)
            cout << "X" <<i+1 << ": " << result[i].real() << " +" << result[i].imag() << " i" << endl;
        else if(result[i].imag() < 0)
            cout << "X" <<i+1 << ": " << result[i].real() << " " << result[i].imag() << " i" << endl;
        else
            cout << "X" <<i+1 << ": " << result[i].real() << endl;
    }
    return 0;
}
