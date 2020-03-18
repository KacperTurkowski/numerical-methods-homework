#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>
using namespace std;

typedef complex<double> Complex;
typedef vector<Complex> Polynomial;
typedef vector<complex<double>> Polynomial_Z;
const double EPS = 1e-13;
const int Max_iteration =1000000;
pair<Polynomial, Complex> horner(Polynomial &polynomial, Complex Z) {
    int size = polynomial.size();
    Polynomial temp = Polynomial(max(1, size - 1));
    for(int i= size - 1; i > 0; i--){
        temp[i - 1]=polynomial[i];
        if(i < size - 1) temp[i - 1]+= temp[i] * Z;
    }
    return make_pair(temp, polynomial[0] + temp[0] * Z);
}
void complete_polynomial(Polynomial &polynomial){
    polynomial.push_back(16);
    polynomial.push_back(-72);
    polynomial.push_back(-28);
    polynomial.push_back(558);
    polynomial.push_back(-990);
    polynomial.push_back(783);
    polynomial.push_back(-486);
    polynomial.push_back(243);
}
void complete_polynomialB(Polynomial &polynomial){
    polynomial.push_back(-4);
    polynomial.push_back(-4);
    polynomial.push_back(-12);
    polynomial.push_back(-8);
    polynomial.push_back(-11);
    polynomial.push_back(-3);
    polynomial.push_back(-1);
    polynomial.push_back(2);
    polynomial.push_back(3);
    polynomial.push_back(1);
    polynomial.push_back(1);
}
void complete_polynomialC(Polynomial_Z &polynomialZ){
    Complex a1(1, 0);//z4 + iz3 − z2 − iz + 1 = 0
    polynomialZ.push_back(a1);
    Complex a2(0, -1);
    polynomialZ.push_back(a2);
    Complex a3(-1, 0);
    polynomialZ.push_back(a3);
    Complex a4(0, 1);
    polynomialZ.push_back(a4);
    Complex a5(1, 0);
    polynomialZ.push_back(a5);
}
Polynomial derivative(Polynomial &polynomial) {
    int n = polynomial.size();
    Polynomial temp=Polynomial(max(1, n - 1));
    for(int i = 1; i < n; i++)
        temp[i - 1]= polynomial[i] * Complex(i);
    return temp;
}
Polynomial derivative_Z(Polynomial_Z &polynomialZ) {
    int n = polynomialZ.size();
    Polynomial temp=Polynomial(max(1, n - 1));
    for(int i = 1; i < n; i++)
        temp[i - 1]= polynomialZ[i] * Complex(i);
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
Complex first_zero_placeZ(Polynomial_Z &polynomialZ, Complex x){
    int n = (int)polynomialZ.size() - 1;
    Polynomial_Z temp1 = derivative_Z(polynomialZ);
    Polynomial_Z temp2 = derivative_Z(temp1);
    for (int iteration = 0; iteration < Max_iteration; iteration++) {
        Complex Y=horner(polynomialZ, x).second;
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
Complex first_zero_place(Polynomial &polynomial, Complex x) {
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
vector<Complex> zeroZ(Polynomial_Z &polynomialZ) {
    vector<Complex> result;
    Polynomial_Z temp = polynomialZ;
    while (temp.size() > 2) {
        Complex z(0.0, 0.0);
        z = first_zero_place(temp, z);
        temp = horner(temp, z).first;
        result.push_back(z);
    }
    result.push_back(-temp[0] / temp[1]);
    return result;
}
vector<Complex> zero(Polynomial &polynomial) {
    vector<Complex> result;
    Polynomial temp = polynomial;

    while (temp.size() > 2) {
        Complex z(0.0, 0.0);
        z = first_zero_place(temp, z);
        temp = horner(temp, z).first;
        result.push_back(z);
    }
    result.push_back(-temp[0] / temp[1]);
    return result;
}
int main(){
    Polynomial A;
    complete_polynomial(A);
    Polynomial B;
    complete_polynomialB(B);
    Polynomial_Z C;
    complete_polynomialC(C);
    cout<<setprecision(10)<<fixed;
    cout<<"Zero(of function a): "<<endl;
    vector<Complex> resultA = zero(A);
    for(int i=0; i < resultA.size(); i++) {
        if(resultA[i].imag() > 0)
            cout << "X" <<i+1 << ": " << resultA[i].real() << " +" << resultA[i].imag() << " i" << endl;
        else if(resultA[i].imag() < 0)
            cout << "X" <<i+1 << ": " << resultA[i].real() << " " << resultA[i].imag() << " i" << endl;
        else
            cout << "X" <<i+1 << ": " << resultA[i].real() << endl;
    }

    cout<<endl<<"Zero (of function b): "<<endl;
    vector<Complex> resultB = zero(B);
    for(int i=0; i < resultB.size(); i++) {
        if(resultB[i].imag() > 0)
            cout << "X" <<i+1 << ": " << resultB[i].real() << " +" << resultB[i].imag() << " i" << endl;
        else if(resultB[i].imag() < 0)
            cout << "X" <<i+1 << ": " << resultB[i].real() << " " << resultB[i].imag() << " i" << endl;
        else
            cout << "X" <<i+1 << ": " << resultB[i].real() << endl;
    }

    cout<<endl<<"Zero (of function c): "<<endl;
    vector<Complex> resultC= zeroZ(C);
    for(int i=0; i < resultC.size(); i++) {
        if(resultC[i].imag() > 0)
            cout << "X" <<i+1 << ": " << resultC[i].real() << " +" << resultC[i].imag() << " i" << endl;
        else if(resultC[i].imag() < 0)
            cout << "X" <<i+1 << ": " << resultC[i].real() << " " << resultC[i].imag() << " i" << endl;
        else
            cout << "X" <<i+1 << ": " << resultC[i].real() << endl;
    }
    return 0;
}
