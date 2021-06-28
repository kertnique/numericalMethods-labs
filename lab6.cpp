#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

struct Bound {
    double low;
    double up;
};

const int num = 5; //степінь поліному + 1
const double coefs[num] = {-6, -2, 1, -3, 1}; // коефи поліному при х^i
double coefs_der[num-1]; //коефи похідної
vector<Bound> bounds; // динамічний масив для меж проміжків
const double precision = 0.0000001; //точність

double f (double x){ // вхідна функція
    double out = 0;
    for(int i = 0; i < num; i++){
        out +=  coefs[i] * pow(x, i);
    }
    return out;
}

double f_der (double x){ // значення похідної
    double out = 0;
    for(int i = 0;i < num-1; i++){
        out += coefs_der[i] * pow(x,i);
    }
    return out;
}

void getDer (void){ // отримуємо похідну з функції
    for(int i = 0; i < num-1; i++){
        coefs_der[i] = coefs[i+1] * (i+1);
    }
}

int rootsSturm (double x[num + 1][num + 1], double a, double b){ // різниця к-ті змін знаку у послідовності для а і b
    double prev_a, prev_b, temp_a, temp_b;
    int counter_a = 0, counter_b = 0;
    for(int i = 0; i <= num; i++){
        temp_a = 0;
        temp_b = 0;
        for(int j = 0; j <= x[i][num]; j++){
            temp_a += x[i][j] * pow(a, j);
            temp_b += x[i][j] * pow(b, j);
        }
        if(i!=0){
            if(temp_a * prev_a < 0) counter_a++;
            if(temp_b * prev_b < 0) counter_b++;
        }
        prev_a = temp_a;
        prev_b = temp_b;
    }
    return abs(counter_a - counter_b);
}

void Sturm (double lower, double upper){ // теорема Штурма
    double seq [num+1][num+1] = {0.0};
    for(int i = 0; i < num-1; i++){
        seq[0][i] = coefs[i];
        seq[1][i] = coefs_der[i];
    }
    seq[0][num - 1] = coefs[num-1];
    seq[0][num] = num-1;
    seq[1][num] = num-2;
    for(int i = 2; i <= num; i++){
            if(seq[i-2][num]>=seq[i-1][num] && seq[i-1][num] > 0){ // ділення многочленів - знаходження остачі
               for(int j = 0; j <= num; j++){
                    seq[i][j] = seq[i-2][j];
               }
               int desired = seq[i-1][num] - 1;// ступінь многочлена, який ми хочемо
               while(seq[i][num] > desired){
                    int ab = seq[i][num] - seq[i-1][num];
                    int pow_1 = seq[i][num];
                    int pow_2 = seq[i-1][num];
                    for(int j = 0; j <= seq[i-1][num]; j++){ // отримуємо остачу
                        seq[i][j+ab] -= seq[i-1][j] * seq[i][pow_1] / seq[i-1][pow_2];
                    }
                    for(int j = num - 1; j >= 0; j--){
                        if(seq[i][j] == 0) seq[i][num] = j-1;
                        else break;
                    }
               }
            }
            else{
            if(seq[i-1][num] > 0){
                for(int j = 0; j <= num; j++){
                    seq[i][j] = seq[i-2][j];
                    }
                }
                else break;
            }
        for(int j = 0; j < num; j++){
            seq[i][j] = -seq[i][j];
        }
    }

    int roots = rootsSturm(seq,lower,upper);
    if(roots == 1){
        Bound x;
        x.low = lower;
        x.up = upper;
        bounds.push_back(x);
    }
    if(roots > 1){
        for(int i = 2; i < 10; i++){
            if(f((lower+upper)/i) != 0){
                Sturm(lower, (lower+upper)/2);
                Sturm((lower+upper)/2, upper);
                break;
            }
        }

    }
}

void borders (void){
    double A1 = abs(coefs[0]);
    double B = abs(coefs[1]);
    double C0 = -1, C1= -1, C2 = -1, C3 = -1;
    double m0 = -1, m1 = -1, m2 = -1, m3 = -1;
    double coefs1[num], coefs2[num], coefs3[num];
    for(int i = 0; i < num; i++){
        coefs1[i] = coefs[num - 1 - i];
        if(i % 2 == 0)   coefs2[i] = -coefs[i];
        else            coefs2[i] = coefs[i];
    }
    for(int i = 0; i < num ; i++)coefs3[i] = coefs2[num - 1 - i];
    if(coefs1[num-1]<0){
        for(int i = 0; i < num ; i++) coefs1[i] = -coefs1[i];
    }
    if(coefs2[num-1]<0){
        for(int i = 0; i < num ; i++) coefs2[i] = -coefs2[i];
    }
    if(coefs3[num-1]<0){
        for(int i = 0; i < num ; i++) coefs3[i] = -coefs3[i];
    }
    for(int i = 0; i < num; i++){
        if(i != num-1 && abs(coefs[i]) > A1) A1 = abs(coefs[i]);
        if(i != 0 && abs(coefs[i]) > B) B = abs(coefs[i]);
        if(coefs[i] < 0) {
            if(abs(coefs[i]) > C0) C0 = abs(coefs[i]);
            m0 = i;
        }
        if(coefs1[i] < 0) {
            if(abs(coefs1[i]) > C1) C1 = abs(coefs1[i]);
            m1 = i;
        }
        if(coefs2[i] < 0) {
            if(abs(coefs2[i]) > C2) C2 = abs(coefs2[i]);
            m2 = i;
        }
        if(coefs3[i] < 0) {
            if(abs(coefs3[i]) > C3) C3 = abs(coefs3[i]);
            m3 = i;
        }
    }
    double lower, upper; // межі додатних коренів
    lower = abs(coefs[0]) / (B + abs(coefs[0])); // т. про границю всіх комплексних коренів рівняння
    upper = (abs(coefs[num-1]) + A1) / abs(coefs[num-1]);
    cout << "За теоремою про границі усіх коренів рівняння:\n X є [ -" << upper << " ; -" << lower << "] u [" << lower << " ; " << upper << "].\n";
    double neg_lower = - upper;
    double neg_upper = -lower;

    double R0 = 1 + pow((C0/coefs[num-1]), (1/(num - 1 - m0))); // т. про верхню межу додатніх коренів
    double R1 = 1 + pow((C1/coefs1[num-1]), (1/(num - 1 - m1))); // т. про верхню межу додатніх коренів
    double R2 = 1 + pow((C2/coefs2[num-1]), (1/(num - 1 - m2))); // т. про верхню межу додатніх коренів
    double R3 = 1 + pow((C3/coefs3[num-1]), (1/(num - 1 - m3))); // т. про верхню межу додатніх коренів
    cout << "За теоремою про верхню межу додатніх коренів:\n X є [ -" << R2 << " ; -" << 1/R3 << "] u [" << 1/R1 << " ; " << R0 << "].\n\n";
    if(R0 < upper) upper = R0;
    if(1/R1 > lower) lower = 1/R1;
    if(-R2 > neg_lower) neg_lower = -R2;
    if(-1/R3 < neg_upper) neg_upper = -1/R3;
    cout << "Межі проміжків, у яких знаходяться корені:\n X є [ " << neg_lower << " ; " << neg_upper << "] u [" << lower << " ; " << upper << "].\n\n";

    int Gua = 0; // теорема Гюа
    for(int i = 1; i < num - 1; i++){
        if(pow(coefs[i],2) < coefs[i-1] * coefs[i+1]){
            Gua = 1;
            break;
        }
    }
    if(Gua == 0){
        int numRoots = num - 1;
        cout << "За теоремою Гюа не визначено існування комплексних коренів, отже максимальна к-ть дійсних коренів : " << numRoots << ".\n\n";
    }
    else{
        int numRoots = num - 3;
        cout << "За теоремою Гюа, існує хоча б одна пара комплексних коренів, отже максимальна к-ть дійсних коренів : " << numRoots << ".\n\n";
    }
    Sturm(neg_lower,neg_upper); //від'ємні проміжки
    Sturm(lower,upper); // додатні проміжки
    cout << "За допомогою теореми Штурма, знайдено інтервали, у кожному з яких є по одному кореню:\n";
    for(int i = 0; i < bounds.size(); i++){
        cout << "X є ( " << bounds[i].low << " ; " << bounds[i].up << ");\n";
    }
    cout << "\n";
}

void Bicection (double a, double b){
    double c = (a + b)/2;
    int steps = 0;
    while (b - a > precision || abs(f(c)) > precision){
        steps ++;
        if(f(c) * f(a) > 0)   a = c;
        else                  b = c;
        c = (a+b) / 2;
    }
    cout << "Проведено " << steps << " ітерацій, знайдено х = " << setprecision(10) << c << ", f(x) = " << f(c) << "\n";
}

void Chordes(double a, double b){
    double c =(a*f(b)-b*f(a))/(f(b)-f(a));
    double c_prev = c + 100; //від балди, але щоб різнилося точно
    int steps = 0;
        while (abs(c - c_prev) > precision || abs(f(c)) > precision){
        steps ++;
        c_prev = c;
        if(f(c) * f(a) > 0)   a = c;
        else                  b = c;
        c = (a*f(b)-b*f(a))/(f(b)-f(a));
    }
    cout << "Проведено " << steps << " ітерацій, знайдено х = " << setprecision(10) << c << ", f(x) = " << f(c) << "\n";
}

void Newton (double a, double b){
    int steps = 0;
    double x = (a+b) / 2; //всеодно яке, головне всередині проміжку
    double x_prev = x + 100; // всеодно яке, головне інше, ніж х
    while(abs(x-x_prev) > precision || abs(f(x)) > precision){
        steps ++;
        x_prev = x;
        x = x_prev - (f(x_prev) / f_der(x_prev));
    }
    cout << "Проведено " << steps << " ітерацій, знайдено х = " << setprecision(10) << x << ", f(x) = " << f(x) <<"\n";
}

int main()
{
    getDer();
    borders();

    cout << "Метод бісекцій:\n";
    for(int i = 0; i < bounds.size(); i++){
        Bicection(bounds[i].low, bounds[i].up);
    }
    cout << "\nМетод хорд:\n";
    for(int i = 0; i < bounds.size(); i++){
        Chordes(bounds[i].low, bounds[i].up);
    }
    cout << "\nМетод Ньотона:\n";
    for(int i = 0; i < bounds.size(); i++){
        Newton(bounds[i].low, bounds[i].up);
    }
    return 0;
}
