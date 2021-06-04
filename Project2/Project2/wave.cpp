#include "wave.h"
void wave::funk(Field* TheField)//funkcia sozdayet massiv clasterov s pomosch'u volnovogo algoritma v poslednem sozdannom binary Find_Clusters'e 
{
    while (1)
    {
        if ((*TheField).Find_Clastersarr[(*TheField).last_Find_Clasters_bin].tip_alg == 'b')
        {
            break;
        }
        (*TheField).last_Find_Clasters_bin++;
    }
    B_Distant_matrix MatrB = (*TheField).Find_Clastersarr[(*TheField).last_Find_Clasters_bin].MatrB;
    int t = 1;
    for (int i = 0; i < (*TheField).numb_Point; i++)
    {
        if ((*TheField).Pointarr[i].a == 0)
        {
            Claster Cl;
            Cl.t_nach = t;
            std::vector <int> numb;
            (*TheField).Pointarr[i].a = t;
            numb.push_back(i);
            for (int j = 0; j < (*TheField).numb_Point; j++)
            {
                if (MatrB.matrx[i][j] == 1 && i != j && (*TheField).Pointarr[j].a == 0)
                {
                    std::cout << (*TheField).Pointarr[j].x << " " << (*TheField).Pointarr[j].y << std::endl;
                    (*TheField).Pointarr[j].b = t + 1;
                }
            }
            std::cout << std::endl;
            int u = 0;
            do
            {
                u = 0;
                t++;
                for (int j = 0; j < (*TheField).numb_Point; j++)
                {
                    if ((*TheField).Pointarr[j].b == t)
                    {
                        (*TheField).Pointarr[j].a = t;
                        numb.push_back(j);
                        for (int v = 0; v < (*TheField).numb_Point; v++)
                        {
                            if (MatrB.matrx[j][v] == 1 && (*TheField).Pointarr[v].a == 0 && (*TheField).Pointarr[v].b == 0 && v != j)
                            {
                                u++;
                                (*TheField).Pointarr[v].b = t + 1;
                            }
                        }
                    }
                }
            } while (u != 0);
            for (int j = 0; j < numb.size(); j++)
            {
                Cl.arr.push_back((*TheField).Pointarr[numb[j]]);
            }
            Cl.t_kon = t;
            (*TheField).Find_Clastersarr[(*TheField).last_Find_Clasters_bin].Clasterarr.push_back(Cl);
            numb.clear();
        }
        t++;
    }
    for (int i = 0; i < (*TheField).numb_Point; i++)
    {
        (*TheField).Pointarr[i].a = 0;
        (*TheField).Pointarr[i].b = 0;
    }
    (*TheField).last_Find_Clasters_bin++;
}
//211-Borovikov-Mikhail-2021