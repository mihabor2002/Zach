#include "dbwave.h"
//211-Borovikov-Mikhail-2021
void dbwave::funk(Field* TheField)//funkcia sozdayet massiv clasterov s pomosch'u volnovogo algoritma v poslednem sozdannom dbbinary Find_Clusters'e 
{
    while (1)
    {
        if ((*TheField).Find_Clastersarr[(*TheField).last_Find_Clasters_db].tip_alg == 'd')
        {
            break;
        }
        (*TheField).last_Find_Clasters_db++;
    }
    Find_Clasters TheFind_Clasters;
    TheFind_Clasters = (*TheField).Find_Clastersarr[(*TheField).last_Find_Clasters_db];
    B_Distant_matrix MatrB;
    MatrB = TheFind_Clasters.MatrB;
    std::vector <int> dbpoint;
    dbpoint = TheFind_Clasters.dbpoint;
    int n = MatrB.n;
    int t = 1;
    for (int i = 0; i < n; i++)
    {
        if ((*TheField).Pointarr[dbpoint[i]].a == 0)
        {
            std::vector <int> numb;
            (*TheField).Pointarr[dbpoint[i]].a = t;
            numb.push_back(dbpoint[i]);
            for (int j = 0; j < n; j++)
            {
                if (MatrB.matrx[i][j] == 1 && i != j && (*TheField).Pointarr[dbpoint[j]].a == 0)
                {
                    (*TheField).Pointarr[dbpoint[j]].b = t + 1;
                }
            }
            int u = 0;
            do
            {
                u = 0;
                t++;
                for (int j = 0; j < n; j++)
                {
                    if ((*TheField).Pointarr[dbpoint[j]].b == t)
                    {
                        (*TheField).Pointarr[dbpoint[j]].a = t;
                        numb.push_back(dbpoint[j]);
                        for (int v = 0; v < n; v++)
                        {
                            if (MatrB.matrx[j][v] == 1 && (*TheField).Pointarr[dbpoint[v]].a == 0 && (*TheField).Pointarr[dbpoint[v]].b == 0 && v != j)
                            {
                                u++;
                                (*TheField).Pointarr[dbpoint[v]].b = t + 1;
                            }
                        }
                    }
                }
            } while (u != 0);
            Claster Cl;
            for (int j = 0; j < numb.size(); j++)
            {
                Cl.arr.push_back((*TheField).Pointarr[numb[j]]);
            }
            (*TheField).Find_Clastersarr[(*TheField).last_Find_Clasters_db].Clasterarr.push_back(Cl);
            numb.clear();
        }
        t++;

    }
    for (int i = 0; i < (*TheField).numb_Point; i++)
    {
        (*TheField).Pointarr[i].a = 0;
        (*TheField).Pointarr[i].b = 0;
        (*TheField).Pointarr[i].DB = 0;
    }
    (*TheField).last_Find_Clasters_db++;
}