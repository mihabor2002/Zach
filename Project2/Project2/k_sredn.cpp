#include "k_sredn.h"
//211-Borovikov-Mikhail-2021
void k_sredn::funk(Field* TheField, int k)
{
    std::vector <Point> centr;
    for (int i = 0; i < k; i++)
    {
        centr.push_back((*TheField).Pointarr[i * ((*TheField).numb_Point / k)]);
    }
    int u;
    do
    {
        u = 0;
        std::vector <float> mass(k);
        int min;
        float dist_min;
        for (int i = 0; i < (*TheField).numb_Point; i++)
        {
            for (int j = 0; j < k; j++)
            {
                mass[j] = dist((*TheField).Pointarr[i], centr[j]);
            }
            dist_min = mass[0];
            min = 0;
            for (int j = 1; j < k; j++)
            {
                if (mass[j] < dist_min)
                {
                    dist_min = mass[j];
                    min = j;
                }
            }
            if ((*TheField).Pointarr[i].for_k_sr != min)
            {
                u++;
                (*TheField).Pointarr[i].for_k_sr = min;
            }
        }
        for (int i = 0; i < k; i++)
        {
            int c = 0;
            float x = 0, y = 0;
            Point ThePoint;
            for (int j = 0; j < (*TheField).numb_Point; j++)
            {
                if ((*TheField).Pointarr[j].for_k_sr == i)
                {
                    c++;
                    x += (*TheField).Pointarr[j].x;
                    y += (*TheField).Pointarr[j].y;
                }
            }
            x = x / c;
            y = y / c;
            ThePoint.x = x;
            ThePoint.y = y;
            centr[i] = ThePoint;
        }
    } while (u != 0);
    std::vector <Claster> clas(k);
    for (int i = 0; i < k; i++)
    {
        Claster TheClaster;
        clas[i] = TheClaster;
    }
    for (int i = 0; i < (*TheField).numb_Point; i++)
    {
        clas[(*TheField).Pointarr[i].for_k_sr].arr.push_back((*TheField).Pointarr[i]);
    }
    for (int i = 0; i < k; i++)
    {
        clas[i].arr.push_back(centr[i]);
    }
    Find_Clasters TheFind_Clasters;
    TheFind_Clasters.tip_alg = 'k';
    for (int i = 0; i < k; i++)
    {
        TheFind_Clasters.Clasterarr.push_back(clas[i]);
    }
    (*TheField).Find_Clastersarr.push_back(TheFind_Clasters);
}