#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f (double x){
    return (x+1) * sin(x);
}

double f_derN (double x, int N){ // похідна 2N-ого порядку
    if(N%2 == 0) return (x+1) * sin(x) - 2 * N * cos(x);
    else return 2 * N * cos(x) - (x+1) * sin (x);
}

int factorial (int x){
    int out = 1;
    for(int i = 2; i <=x; i++){
        out*=i;
    }
    return out;
}

double getMax (double a, double b, int m){
    int N = 100000;
    double maximum = max(abs(f_derN(a,m)), abs(f_derN(b,m)));
    for(int i = 1; i < N ; i++){
        if(abs(f_derN(a + i * (b-a) / N,m)) > maximum)  maximum = f_derN(a + i * (b-a) / N,m);
    }
    return maximum;
}

double Simpson (double a, double b, int n){
    double result = f(a);
    for(int i = 1; i <= n; i++){
        result += 4 * f(a + (b-a) * (2*i-1)/ 2 / n);
        result += 2 * f(a + (b-a) * i/n);
    }
    result -= f(b);
    result *= (b-a) / 6 / n;
    return result;
}

int error_Simpson (double a, double b, double wanted_error){    // виводить необхідне 2n для заданої точності
    double maximum = getMax(a,b,2);
    int n = ceil(pow(maximum * pow((b-a), 2) / 180 / wanted_error , 0.25));
    if(n%2 ==1) n++;
    return n;
}

double Gaus (double a, double b){
    double x[4] = {-0.861136311594, -0.339981043584856, 0.339981043584856, 0.861136311594}; // з таблиці коренів поліномів Лежандра
    double w[4] = {0.34785484513745, 0.6521451548625, 0.6521451548625, 0.34785484513745};
    double result = 0;
    for(int i = 0; i < 4; i++){
        double y = x[i] * (b-a)/2 + (a+b)/2;
        result += w[i] * f(y);
    }
    return result * (b-a) / 2;
}

int error_Gaus (double a, double b, double wanted_error){
    int n = 2;
    double R = wanted_error + 1;
    while (R > wanted_error){
        R = (pow(factorial(n),4) * pow((b-a), 2*n+1)) / ((2*n+1) * pow(factorial(2*n), 3)) * (getMax(a,b,n));
        n++;
    }
    return n;
}

int main()
{
    double lower = 0.8;
    double upper = 1.6;
    double precision = 0.0001;
    cout << "Метод Сімпсона: \nЗа формулою похибки знадено необхідну кількість підінтервалів: " << 2 * error_Simpson(lower,upper,precision);
    cout << "\nАналітична похибка для методу: " << pow((upper - lower),2) / 180 / pow(2 * error_Simpson(lower,upper,precision), 4) * getMax(lower,upper,2);
    cout << "\nЗнайдений інтеграл : " << setprecision(11) << Simpson(lower,upper, error_Simpson(lower,upper,precision));
    int n = error_Gaus(lower, upper, precision);
    cout << "\n\nМетод Гауса: \nЗа формулою похибки знайдено необхідну кількість підінтервалів: " << n;
    cout << "\nАналітична похибка для методу: " <<  (pow(factorial(n),4) * pow((upper-lower), 2*n+1)) / ((2*n+1) * pow(factorial(2*n), 3)) * (getMax(lower,upper,n));
    cout << "\nЗнайдений інтеграл: " << setprecision(11) << Gaus(lower,upper) << "\n";
    return 0;
}
