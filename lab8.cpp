#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>

using namespace std;

double arr_f[4] = {0.0, 0.0, 0.0, 0.0};
double arr_y[4] = {0.0, 0.0, 0.0, 0.0};
vector<double> y_Runge;
vector<double> y_Adams;

double f (double x, double y){
    double a = 1 + 0.4 * 3;
    return cos(a*x + y) + (x-y);
}

void Runge_Kutt (double a, double b, double h, bool out){
    double i = a;
    double y = 0; // значення функції в х=а
    while(i < b){
        if(i < 4*h){
            int v = i/h;
            arr_y[v] = y;
            arr_f[v] = h * f(i,y);
        }
        double k0 = h * f(i,y);
        double k1 = h * f(i + h/2, y + k0/2);
        double k2 = h * f(i + h/2, y + k1/2);
        double k3 = h * f(i + h, y + k2);
        double useful = k0 + 2*k1 + 2 * k2 + k3;

        if(out){
        double tau = abs((k2-k1)/(k1-k0));
        double epsilon = (y - y_Runge[(int)(round(i/h))]) / 31;
        cout << fixed << setprecision(2) << "x = " << i;
        cout << setprecision(7) << ";  y(x) = " << y << ";  ";
        cout << setprecision(10) << "τ = " << tau << ";  ε = "<< epsilon <<"\n";
        }
        else{
            int v = round(i/h);
            if(v % 2 == 0) y_Runge.push_back(y);
        }
        y += (useful/6);
        i += h;
    }
}

void Adams (double a, double b, double h, bool out){
    double i = a + 4 * h;
    while(i < b){
    double y = arr_y[3] + h/24 * (55 * arr_f[3] - 59 * arr_f[2] + 37 * arr_f[1] - 9 * arr_f[0]);
    if(out){
        cout << setprecision(2) << "x = " << i;
        cout << setprecision(7) << "  y ext = " << y;
    }
    for(int j = 0; j < 3; j++){
        arr_f[j] = arr_f[j+1];
        arr_y[j] = arr_y[j+1];
    }
    arr_y[3] = y;
    arr_f[3] = f(i,y);
    y = arr_y[2] + h/24 * (9 * arr_f[3] + 19 * arr_f[2] - 5 * arr_f[1] + arr_f[0]);
    arr_y[3] = y;
    arr_f[3] = f(i,y);
    if(out){
        cout << setprecision(7) << "  y int = " << y;
        cout << setprecision(6) << " ε = " << (y - y_Adams[round(i/h)]) / 15 << "\n";
    }
    else{
        int v = round(i/h);
        if(v % 2 == 0) y_Adams.push_back(y);
    }
    i += h;
    }
}

int main()
{
    double h = 0.02;
    double a = 0;
    double b = 6;
    cout << "Метод Рунге-Кутта:\n\n";
    Runge_Kutt(a,b,h/2, false);
    Adams(a,b+ 2*h,h/2, false);
    //cout << y_Runge.size() << "  " << y_Adams.size() << endl;
    Runge_Kutt(a,b,h, true);
    cout << "\nМетод Адамса:\n\n";
    Adams(a,b,h, true);
    cout << y_Runge.size();
    return 0;
}
