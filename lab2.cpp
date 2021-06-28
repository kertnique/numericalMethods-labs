#include <iostream>
#include <cmath>
#include <string>

using namespace std;

const int num = 5; //розмірність матриці

bool isSymmetrical (const float matrix[num][num]){
    bool symmetry = true;
    for(int i = 0; i < num; i++){
        for(int j = 0; j < i; j++){
            if(matrix[i][j] != matrix[j][i]) symmetry = false;
        }
    }
    return symmetry;
}

void outForwardSqrt (float a[num][num]){ //виведення ходу
    cout << "Результат прямого ходу: створено трикутну матрицю:\n";
    for(int i = 0; i<num; i++){
        for(int j=0; j<num; j++){
            string output = to_string(a[i][j]);
            cout << output << "  ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void outReverseSqrt (float y[num]){ //виведення ходу
    cout << "Результат зворотнього ходу: створено вектор коефіцієнтів у:\n";
    for(int i = 0; i < num; i++){
        string output = to_string(y[i]);
        cout << output << " ";
    }
    cout << "\n\n";
}

void outFinal (float y[num]){ //виведення результатів
    cout << "Знайдено розв'язки СЛАР\n";
    for(int i = 0; i < num; i++){
        cout << "x" << i+1 << " = " << y[i] << "\n";
    }
}

void outGauss (float a[num][num + 1]){ //виведення результатів

}

void nezviazka (const float a[num][num],const float b[num], float x[num]){ //виведення вектору незв'язки
    cout<<"\nЗнайдено вектор незв'язки:\n";
    for(int i = 0;i < num; i++){
        float output = 0;
        for(int j = 0; j < num; j++){
            output+= a[i][j] * x[j];
        }
        output -= b[i];
        cout<< output << ";\n";
    }
}

void squareRootMethod (const float matrix[num][num], const float vektor[num]){
    cout << "Матриця симетрична, обрано методу квадратних коренів.\n\n";
    float triang[num][num]; // трикутна матриця
    for(int i = 0; i<num; i++){ //починаєм прямий хід
        for(int j = 0; j < num; j++){
            if(i <= j){
                triang[i][j] = matrix[i][j];
                for(int k =0 ; k < i; k++)  triang[i][j] -= triang[k][i] * triang[k][j];
            if(i == j)  triang[i][j] = sqrt(triang[i][j]);
            else  triang[i][j] /= triang[i][i];
            }
            if(i > j) triang[i][j] = 0;
        }
    }
    outForwardSqrt(triang); // прямий хід закінчився
    float coefs[num]; // починаємо зворотній хід
    for(int i = 0; i < num; i++){
        coefs[i] = vektor[i];
        for(int k = 0; k < i; k++) coefs[i] -= triang[k][i] * coefs[k];
        coefs[i] /= triang[i][i];
    }
    outReverseSqrt(coefs); // зворотній хід закінчився
    float answers[num]; // останній крок, знаходимо відповіді
    for(int i = num - 1; i >= 0; i--){
        answers[i] = coefs[i];
        for(int k = num - 1; k > i; k--) answers[i] -= triang[i][k] * answers[k];
        answers[i] /= triang[i][i];
    }
    outFinal(answers);
    nezviazka(matrix,vektor,answers);
}


void gaussMethod (const float matrix[num][num], const float vektor[num]){
    cout << "Матриця несиметрична, обрано методу Гауса.\n\n";
    float nMatrix[num][num + 1]; //створення розширеної матриці
    for(int i = 0; i < num; i++){
        for(int j = 0 ; j < num; j++){
            nMatrix[i][j] = matrix[i][j];
        }
        nMatrix[i][num] = vektor [i];
    }
    for(int i = 0; i < num; i++){    //прямий хід
        for(int j = num; j >= i; j--){
            nMatrix[i][j] /= nMatrix[i][i];
        }
        for(int k = i+1; k < num ; k++){
            for(int j = num; j >= i; j--){
                nMatrix[k][j] -= nMatrix[k][i] * nMatrix[i][j];
            }
        }
        outGauss(nMatrix);
    } // прямий хід закінчився
    for(int i = num-1; i >= 0; i--){ //початок зворотнього ходу
        for(int k = i - 1; k >= 0; k--){
            for(int j = num; j >= i; j--){
                nMatrix[k][j] -= nMatrix[i][j] * nMatrix[k][i];
            }
        }
        outGauss(nMatrix);
    } //кінець зворотнього ходу
    float answers[num]; //відповіді
    for(int i = 0;i< num ;i++) answers[i] = nMatrix[i][num];
    outFinal(answers);
    nezviazka(matrix,vektor,answers);
}

int main()
{
    const float a = 0.5;
    const float b = 0.7;
    const float matrix[num][num] = {
        {5.18 + a, 1.12, 0.95, 1.32, 0.83},
        {1.12, 4.28 - a, 2.12, 0.57, 0.91},
        {0.95, 2.12, 6.13 + a, 1.29, 1.57},
        {1.32, 0.57, 1.29, 4.57 - a, 1.25},
        {0.83, 0.91, 1.57, 1.25, 5.21 + a}
    };
    const float vektor[num] = {6.19 + b, 3.21, 4.28 - b, 6.25, 4.95 + b};
    if(isSymmetrical(matrix) == true) squareRootMethod(matrix,vektor); //перевіряємо на симетричність
    else gaussMethod(matrix,vektor);
    return 0;
}
