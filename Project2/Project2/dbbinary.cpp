#include "dbbinary.h"
//211-Borovikov-Mikhail-2021
void dbbinary::funk(Field* TheField, float alph, int k)//funkciya snachala zapuskayet DBSCAN, a potom sozdayet Find_Clasters s matritsey MatrB, postroennoy po vnutrennim i pereferiinim tochkam
{
    B_Distant_matrix MatrxB;//binarnaya matritsa vsego polya dlya raspredeleeniya tochek
    for (int i = 0; i < (*TheField).numb_Point; i++)
    {
        for (int j = 0; j < i; j++)
        {
            MatrxB.matrx[j].push_back((*TheField).MatrD.matrx[i][j] < alph);
        }
        std::vector <bool> v;
        for (int j = 0; j < i; j++)
        {
            v.push_back(MatrxB.matrx[j][i]);
        }
        v.push_back(1);
        MatrxB.matrx.push_back(v);
    }
    for (int i = 0; i < (*TheField).numb_Point; i++)
    {
        int t = 0;
        for (int j = 0; j < (*TheField).numb_Point; j++)
        {
            t += MatrxB.matrx[i][j];
        }
        if (t >= k)
        {
            (*TheField).Pointarr[i].DB = 1;
        }
    }
    for (int i = 0; i < (*TheField).numb_Point; i++)
    {
        int t = 0;
        if ((*TheField).Pointarr[i].DB != 1)
        {
            for (int j = 0; j < (*TheField).numb_Point; j++)
            {
                if (MatrxB.matrx[i][j] == 1 && (*TheField).Pointarr[j].DB == 1)
                {
                    (*TheField).Pointarr[i].DB = 2;
                    break;
                }
            }
        }
    }
    int n = 0, s = 0;
    std::vector <int> dbpoint;
    B_Distant_matrix MatrB;//matritsa tochek,kotorie ne pil'
    for (int i = 0; i < (*TheField).numb_Point; i++)
    {
        if ((*TheField).Pointarr[i].DB != 0)
        {
            n++;
            std::vector <int> vspom;
            dbpoint.push_back(i);
            int t = 0;
            for (int j = 0; j < i; j++)
            {
                if ((*TheField).Pointarr[j].DB != 0)
                {
                    MatrB.matrx[t].push_back((*TheField).MatrD.matrx[j][i] < alph);
                    t++;
                }
            }
            std::vector <bool> v;
            std::cout << t << " ";
            t = 0;
            std::cout << t;
            for (int j = 0; j < i; j++)
            {
                if ((*TheField).Pointarr[j].DB != 0)
                {
                    v.push_back(MatrB.matrx[t][n - 1]);
                    t++;
                }
            }
            v.push_back(1);
            MatrB.matrx.push_back(v);
        }
    }
    MatrB.n = n;
    Find_Clasters TheFind_Clasters;
    TheFind_Clasters.tip_alg = 'd';
    TheFind_Clasters.MatrB = MatrB;
    TheFind_Clasters.dbpoint = dbpoint;
    (*TheField).Find_Clastersarr.push_back(TheFind_Clasters);
}