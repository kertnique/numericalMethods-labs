#include <iostream>
#include <string>
#include <cmath>


using namespace std;

const int num = 4;

void outMatrix (float a[num][num], float b[num]){
    cout << "Отримано матрицю з діагональною перевагою:\n\n";
    for(int i = 0; i < num; i++){
        for(int j = 0; j <= num; j++){
            if(j == num) cout << " | " << to_string(b[i]);
            else cout<< to_string(a[i][j]) << " ";
            }
        cout<< "\n";
        }
    cout<< "\n";
}

void outAnswers (float x[num], int i){ // виведення вектору шуканих змінних
    cout << "Результат " << i << "-ого циклу:\n";
    for(int j = 0; j < num; j++){
        cout << "X" << j+1 << " = " << x[j] <<";\n";
    }
    cout << "\n";
}

void outNeviazka (float a[num][num], float b[num], float x[num]){ //виведення вектору незв'язки
    cout<<"Знайдено вектор незв'язки:\n";
    for(int i = 0;i < num; i++){
        float output = 0;
        for(int j = 0; j < num; j++){
            output -= a[i][j] * x[j];
        }
        output += b[i];
        cout<< output << ";\n";
    }
    cout << "\n";
    }

int exit (int a[num]){
    for(int i = 0; i < num; i++){
        if(a[i] != 1) return 0;
    }
    return 1;
}

void zeidel (float a[num][num],float b[num], const float precision){
    cout << "Обрано методу Зейделя\n";
    for(int i = 0; i < num; i++){
        for(int j = 0; j < num; j++){
            if(i != j)  a[i][j] /= a[i][i];
        }
        b[i] /= a[i][i];
        a[i][i] = 1;
    }
    float answers[num] = {1.0}; // початкові значення іксів  - одиниця та нулі
    int step = 0;
    int ifexit[num] = {0}; //чи вже досягла і-та змінна критерії зупинки
    outAnswers(answers,step);
    while(exit(ifexit) == 0){
        step++;
        for(int i = 0; i < num; i++){
            float temp = b[i];
            for(int j = 0; j < num ; j++){
                if(i != j) temp -= (a[i][j] * answers[j]);
            }
           if(answers[i] - temp <= precision && answers[i] - temp >= -precision)    ifexit[i] = 1;
           answers[i] = temp;
        }
        if(step <= 3)  {
        outAnswers (answers, step);
        outNeviazka(a,b,answers);
        }
    }
    outAnswers(answers, step);
    outNeviazka(a,b,answers);
}

int main()
{
    const float a = 0.4;
    const float b = 0.4;
    float matrix[num][num] = {
    {8.30, 2.62 + a, 4.10, 1.90},
    {3.92, 8.45, 8.78 - a, 2.46},
    {3.77, 7.21 + a, 8.04, 2.28},
    {2.21, 3.65 - a, 1.69, 6.99}};
    float vektor[num] = {-10.65 + b, 12.21, 15.45 - b, -8.35};
    const float precision = 0.0000001;
    float matrixNew[num][num] = {0};
    float vektorNew[num] = {0};
    for(int i=0; i< num; i++){
        for(int j=0; j< num; j++){
            for(int k = 0; k < num; k++){
                matrixNew[i][j] += matrix[k][i] * matrix[k][j];
            }
            vektorNew[i] += matrix[j][i] * vektor[j];
        }
    }
    zeidel(matrix, vektor, precision);
    return 0;
}
