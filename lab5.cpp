#include <iostream>
#include <cmath>

using namespace std;

const int num = 5; // number of points


float f (float x){ // given function
    return x/2*cos(3*x);
}

float u_2 (float x1, float x2, float y1, float y2){ //роздільна різниця двох
    return (y1 - y2)/(x1 - x2);
}

float u_3 (float x1, float x2, float x3, float y1, float y2, float y3){ // роздільна різниця трьох
    return (u_2(x1, x2, y1, y2) - u_2(x2, x3, y2, y3)) / (x1 - x3);
}

void outBase (float l[num][num]){
    cout << "Базисні поліноми:\n";
    for(int i = 0; i< num; i++){
        cout << "l" << i+1 << " = ";
        for(int j = 0; j < num ; j ++){
            cout << l[i][j] << " * x^" << j;
            if(j != num - 1) cout << " + ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void outPolinome(float x[num]){
    cout << "Поліном Лагранжа:\n";
    for(int i = 0; i < num; i++){
        cout << x[i] << " * x^" << i;
        if(i != num - 1) cout<< " + ";
    }
    cout << "\n\n";
}

void outPochybka (float polinome[num], float x[num]){ //тут х - це коефи поліному
    cout << "Похибки знайденого поліному для заданих значень:\n";
    for(int i = 0; i < num; i++){
       float pochybka  = f(x[i]);
       for(int j = 0; j < 4; j++){
            pochybka -= polinome[j] * pow(x[i], j);
       }
       cout << "e" << i+1 << " = " <<   abs(pochybka) << ";\n";
    }
    cout << "\n";
}

void Lagrange (float x[num], float y[num]){
    float l[num][num] = {0.0}; //базисні поліноми для кожного ікса, заповнений нулями
    for(int i = 0; i < num; i++){ // знаходимо і-товий базисний поліном
        float L_i[num] = {0};
        L_i[0] = 1; // створили тимчасовий базовий поліном l = 1
        for(int j = 0; j < num; j++){
            if(i != j){
                float iter[2] = {-x[j]/(x[i] - x[j]), 1/(x[i] - x[j])};
                float temp_L[num];
                for(int a = 0; a < num; a++){ //створюємо temp_L - копію L_i
                    temp_L[a] = L_i[a];
                    L_i[a] = 0;
                }
                for(int a = 0; a < num; a++){ //множення поліномів temp_L i iter
                    for(int b = 0; b < 2; b++){
                        if(a+b < num) L_i[a+b] += temp_L[a] * iter[b];
                    }
                }
            }
        }
        for(int j = 0; j < num; j++){
            l[i][j] = L_i[j];
        }
    }
    //outBase(l);
    float polinome[num] = {0}; // фінальний поліном
    for(int i = 0; i < num; i++){
        for(int j = 0; j < num; j++){
            polinome[j] += y[i]*l[i][j];
        }
    }
    outPolinome(polinome);
    outPochybka(polinome, x);
    cout << "Похибки в точках для поліному Лагранжа:\n";
    for(int i = 0; i <= num * 8; i++){
        float input = x[0] + i * (x[num - 1] - x[0]) / 8 / num;
        float out = f(input);
        for(int j = 0; j < 5; j++){
            out -= polinome[j] * pow(input, j);
        }
        cout << "x = " << input << ", e = " << out << ";\n";
    }
    cout << "\n";
}

void outSplines (float a[num-1][4], float x[num]){
    cout << "Сплайни на проміжках:\n";
    for(int i = 0; i < num - 1; i++){
        cout << "Сплайн між точками: " << x[i] << " , "<< x[i+1] << " є поліномом:\n";
        for(int j = 0; j < 4; j++){
        cout << a[i][j] << " * x^" << j;
        if(j != num - 2) cout<< " + ";
    }
    cout << "\n\n";
    }
}

void Splines (float x[num], float y[num]){
    float coefs[num - 1][4] = {0}; // фінальні коефіцієнти
    float nMatrix[num - 2][num - 1]; // знаходимо коефи с
    for(int i = 0; i < num-2; i++){
        for(int j = 0; j < num-1; j++){
            if(j == i - 1)   nMatrix[i][j] = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
            if(j == i)       nMatrix[i][j] = 2;
            if(j == i + 1)   nMatrix[i][j] = (x[i+2] - x[i+1]) / (x[i+2] - x[i]);
            if(j == num - 2) nMatrix[i][j] = 6 * u_3(x[i],x[i+1],x[i+2],y[i],y[i+1],y[i+2]);
        }
    } // тут починається Гаус
    for(int i = 0; i < num - 2; i++){    //прямий хід
        for(int j = num - 2; j >= i; j--){
            nMatrix[i][j] /= nMatrix[i][i];
        }
        for(int k = i+1; k < num - 2; k++){
            for(int j = num - 2; j >= i; j--){
                nMatrix[k][j] -= nMatrix[k][i] * nMatrix[i][j];
            }
        }
    } // прямий хід закінчився
    for(int i = num - 3; i >= 0; i--){ //початок зворотнього ходу
        for(int k = i - 1; k >= 0; k--){
            for(int j = num - 2; j >= i; j--){
                nMatrix[k][j] -= nMatrix[i][j] * nMatrix[k][i];
            }
        }
    } //кінець зворотнього ходу
    float outs[num - 2]; //відповіді - коефи с
    for(int i = 0;i < num - 2;i++) outs[i] = nMatrix[i][num - 2];
    // тут закінчився Гаус
    for(int i = 0; i < num - 1; i++){
        coefs[i][0] = y[i];
        if(i != num - 2) coefs[i][2] = outs[i];
        if(i == 0)  {
        coefs[i][3] = coefs[i][2]/(x[i+1] - x[i]);
        coefs[i][1] = (x[i+1] - x[i])*coefs[i][2]/3 + u_2(x[0],x[1],y[0],y[1]);
        }
        else{
        coefs[i][3] = (coefs[i][2] - coefs[i-1][2]) / (x[i+1] - x[i]);
        coefs[i][1] = (x[i+1] - x[i]) / 6 * (coefs[i][2] * 2 + coefs[i-1][2]) + u_2(x[i],x[i+1],y[i],y[i+1]);
        }
    }
    float polinome[num-1][4] = {0}; //вихідні поліноми
    for(int i = 0; i < num - 1; i++){
        polinome[i][1] = coefs[i][1] - 2 * x[i+1] * coefs[i][2] / 2 + 3 * pow(x[i+1],2) * coefs[i][3] / 6;
        polinome[i][2] = coefs[i][2] / 2 - 3 * x[i+1] * coefs[i][3] / 6;
        polinome[i][3] = coefs[i][3] / 6;
        polinome[i][0] = y[i+1] - polinome[i][1] * x[i+1] - polinome[i][2] * pow(x[i+1], 2) - polinome[i][3] * pow(x[i+1], 3);
    }
    outSplines(polinome, x);
    cout << "Похибки в точках для сплайнів:\n";
    for(int i = 0; i <= 8 * num ; i++){
        float input = x[0] + i * (x[num - 1] - x[0]) / 8 / num;
        float output = f(input);
        for(int j = 0; j < 4; j++){
            int pol = (i-1)/8; //цілочисельне обов'язково
            if(pol > 0) pol = pol - 1;
            output -= polinome[pol][j] * pow(input, j);
        }
        cout << "x = " << input << "; e = " << output << ";\n";
    }
}

int main()
{
    float args[num] = {-1,1,3,5,7}; // inputs
    float outputs[num];
    cout << "Значення початкової функції:\n";
    for(int i = 0; i < num; i++){
        outputs[i] = f(args[i]);
        cout << "x" << i+1 << ": "<< args[i] << ", y" << i+1 << ": " << outputs[i] << ";\n";
        // outputs like: "x1 = -4, y1 = 0.291;"
    }
    cout << "\n";
    Lagrange(args,outputs);
    Splines(args, outputs);
    return 0;
}
