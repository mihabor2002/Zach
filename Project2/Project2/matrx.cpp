#include "matrx.h"
//211-Borovikov-Mikhail-2021
void matrx::funk(Field* TheField)
{
    (*TheField).MatrD.n = (*TheField).numb_Point;
    int n = (*TheField).numb_Point;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            (*TheField).MatrD.matrx[j].push_back(dist((*TheField).Pointarr[i], (*TheField).Pointarr[j]));
        }
        std::vector <float> v;
        for (int j = 0; j < i; j++)
        {
            v.push_back((*TheField).MatrD.matrx[j][i]);
        }
        v.push_back(0);
        (*TheField).MatrD.matrx.push_back(v);
    }
}