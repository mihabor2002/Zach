#include "binary.h"
//211-Borovikov-Mikhail-2021
void binary::funk(Field* TheField, float alph)//funkciya sozdayet Find_Clasters s binarnoi matritsei so znacheniem alpha
{
    Find_Clasters TheFind_Clasters;//find_clasters dlya dannogo razbieniya
    B_Distant_matrix MatrB;//binarnaya matritsa
    int n = (*TheField).MatrD.n;
    MatrB.n = n;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            MatrB.matrx[j].push_back((*TheField).MatrD.matrx[i][j] < alph);
        }
        std::vector <bool> v;
        for (int j = 0; j < i; j++)
        {
            v.push_back(MatrB.matrx[j][i]);
        }
        v.push_back(1);
        MatrB.matrx.push_back(v);
    }
    TheFind_Clasters.MatrB = MatrB;
    TheFind_Clasters.tip_alg = 'b';
    (*TheField).Find_Clastersarr.push_back(TheFind_Clasters);
}