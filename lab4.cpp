#include <iostream>

using namespace std;

const int num = 4;

void outMatrix (float a[num][num], string str, int step){
    cout << "Отримано матрицю " << str << " " << step << ":\n\n";
    for(int i = 0; i < num; i++){
        for(int j = 0; j < num; j++){
             cout<< to_string(a[i][j]) << " ";
            }
        cout<< "\n";
        }
    cout<< "\n";
}



void toFrobenius (float a[num][num]){
    for(int k = num-1; k > 0 ; k--){
        if(a[k][k-1] != 0){
        float M[num][num];
        float M1[num][num];
        for(int i = 0; i < num; i++){
        for(int j = 0; j < num; j++){
            if(i == k - 1){ // заповнення матриць М і М^-1 для к-того кроку
                if(j == i) M[i][j] = 1/a[k][k-1];
                else M[i][j] = -a[k][j]/a[k][k-1];
                    M1[i][j] = a[k][j];
                    }
                else{
                    if(j == i)  M[i][j] = 1;
                    else M[i][j] = 0;
                    M1[i][j] = M[i][j];
                }
            }
        }
            outMatrix(M, "M для кроку", num-k);
            outMatrix(M1, "M^-1 для кроку", num-k);

            float mult[num][num];
            for(int i = 0; i < num; i++){ // M^-1 * A
                for(int j = 0; j < num; j++){
                    float c = 0;
                    for(int z = 0; z < num; z++){
                        c+= M1[i][z] * a[z][j];
                    }
                    mult[i][j] = c;
                }
            }
             for(int i = 0; i < num; i++){ // A * M
                for(int j = 0; j < num; j++){
                    float c = 0;
                    for(int z = 0; z < num; z++){
                        c+= mult[i][z] * M[z][j];
                    }
                    a[i][j] = c;
                }
            }
        }
        else{
            cout << "На "<< k <<" кроці виникла проблема: ділення на нуль.\n";
            break;
        }
    }
    a[1][3] = 0;
    outMatrix(a, "у нормальній формі Фробеніюса", 0);
}

int main()
{
    const float a = 0.33;
    const float b = -0.04;
    const float g = -0.04;
    const float d = 0.045;
    float matrix[num][num] = {
    {6.26 + a, 1.10 - b, 0.97 + g, 1.24 - d},
    {1.10 - b, 4.16 - a, 1.30, 0.16},
    {0.97 + g, 1.30, 5.44 + a, 2.10},
    {1.24 - d, 0.16, 2.10, 6.10 - a}};
    toFrobenius(matrix);
    return 0;
}
