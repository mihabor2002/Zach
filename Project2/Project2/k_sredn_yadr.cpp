#include "k_sredn_yadr.h"
//211-Borovikov-Mikhail-2021
void k_sredn_yadr::funk(Field* ThisField, int k, int p)
{
    std::vector <Point> centr;
    for (int i = 0; i < k; i++)
    {
        centr.push_back((*ThisField).Pointarr[i * ((*ThisField).numb_Point / k)]);
    }
    for (int i = 0; i < (*ThisField).numb_Point; i++)
    {
        std::vector <float> mass(k);
        for (int j = 0; j < k; j++)
        {
            mass[j] = dist((*ThisField).Pointarr[i], centr[j]);
        }
        float dist_min = mass[0];
        int min = 0;
        for (int j = 1; j < k; j++)
        {
            if (mass[j] < dist_min)
            {
                dist_min = mass[j];
                min = j;
            }
        }
        if ((*ThisField).Pointarr[i].for_k_sr != min)
        {
            (*ThisField).Pointarr[i].for_k_sr = min;
        }
        mass.clear();
    }
    for (int q = 0; q < 100; q++)
    {
        std::vector<std::vector<Point>> yadr(k);
        std::vector<std::vector<Point>> clast(k);
        for (int i = 0; i < (*ThisField).numb_Point; i++)
        {
            clast[(*ThisField).Pointarr[i].for_k_sr].push_back((*ThisField).Pointarr[i]);
        }
        for (int c = 0; c < k; c++)
        {
            Field TheField;
            for (int u = 0; u < clast[c].size(); u++)
            {
                TheField.Pointarr.push_back(clast[c][u]);
                TheField.numb_Point++;
            }
            k_sredn Thek_sredn;
            Thek_sredn.funk(&TheField, p);
            for (int u = 0; u < p; u++)
            {
                yadr[c].push_back(TheField.Find_Clastersarr[0].Clasterarr[u].arr[TheField.Find_Clastersarr[0].Clasterarr[u].arr.size() - 1]);
            }
        }
        for (int i = 0; i < (*ThisField).numb_Point; i++)
        {
            std::vector<std::vector<float>> rasst(k);
            for (int c = 0; c < k; c++)
            {
                for (int u = 0; u < p; u++)
                {
                    rasst[c].push_back(dist((*ThisField).Pointarr[i], yadr[c][u]));
                }
            }
            int cl = 0;
            float min = rasst[0][0];
            for (int c = 0; c < k; c++)
            {
                for (int u = 0; u < p; u++)
                {
                    if (rasst[c][u] < min)
                    {
                        cl = c;
                        min = rasst[c][u];
                    }
                }
            }
            (*ThisField).Pointarr[i].for_k_sr = cl;
            rasst.clear();
        }
        std::ofstream outpu;
        outpu.open("output_for_ya" + std::to_string(q + 2) + ".txt");
        for (int i = 0; i < (*ThisField).numb_Point; i++)
        {
            outpu << (*ThisField).Pointarr[i].x << " " << (*ThisField).Pointarr[i].y << " " << (*ThisField).Pointarr[i].for_k_sr << std::endl;
        }
        outpu.close();
        outpu.open("output_for_ya_ya" + std::to_string(q + 2) + ".txt");
        for (int c = 0; c < k; c++)
        {
            for (int u = 0; u < p; u++)
            {
                outpu << yadr[c][u].x << " " << yadr[c][u].y << std::endl;
            }
        }
        yadr.clear();
        clast.clear();
    }
}