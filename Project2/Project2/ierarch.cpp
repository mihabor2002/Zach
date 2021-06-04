#include "ierarch.h"
//211-Borovikov-Mikhail-2021
void ierarch::funk(Field* ThisField)
{
    std::vector <std::vector<Point>> clast;
    std::vector <std::vector<float>> rasst;
    std::vector <Point> sered;
    std::vector <std::vector<Point>> otrez;
    std::ofstream output;
    output.open("output_for_pointi.txt");
    for (int i = 0; i < (*ThisField).numb_Point; i++)
    {
        output << (*ThisField).Pointarr[i].x << " " << (*ThisField).Pointarr[i].y << std::endl;
    }
    output.close();
    std::vector <Point> kek;
    for (int i = 0; i < (*ThisField).numb_Point; i++)
    {
        kek.push_back((*ThisField).Pointarr[i]);
        clast.push_back(kek);
        kek.clear();
    }
    std::cout << "kek";
    for (int i = 0; i < (*ThisField).numb_Point; i++)
    {
        for (int j = 0; j < i; j++)
        {
            rasst[j].push_back(dist((*ThisField).Pointarr[i], (*ThisField).Pointarr[j]));
        }
        std::vector <float> vspom;
        for (int j = 0; j < i; j++)
        {
            vspom.push_back(dist((*ThisField).Pointarr[i], (*ThisField).Pointarr[j]));
        }
        vspom.push_back(0);
        rasst.push_back(vspom);
        vspom.clear();
    }
    int u = 0;
    std::cout << "pep";
    while (clast.size() != 1)
    {
        int mass[2];
        float k = rasst[0][1];
        mass[0] = 0;
        mass[1] = 0;
        for (int i = 0; i < rasst.size(); i++)
        {
            for (int j = 0; j < rasst[i].size(); j++)
            {
                if (rasst[i][j] < k && i != j)
                {
                    k = rasst[i][j];
                    mass[0] = i;
                    mass[1] = j;
                }
            }
        }
        Point C1, C2, C3;
        C1.x = 0;
        C1.y = 0;
        C2.x = 0;
        C2.y = 0;
        for (int i = 0; i < clast[mass[0]].size(); i++)
        {
            C1.x += clast[mass[0]][i].x;
            C1.y += clast[mass[0]][i].y;
        }
        C1.x /= clast[mass[0]].size();
        C1.y /= clast[mass[0]].size();
        for (int i = 0; i < clast[mass[1]].size(); i++)
        {
            C2.x += clast[mass[1]][i].x;
            C2.y += clast[mass[1]][i].y;
        }
        C2.x /= clast[mass[1]].size();
        C2.y /= clast[mass[1]].size();
        C3.x = 0;
        C3.y = 0;
        C3.x = (C1.x + C2.x) / 2;
        C3.y = (C1.y + C2.y) / 2;
        std::vector <Point> zm;
        zm.push_back(C1);
        zm.push_back(C2);
        otrez.push_back(zm);
        sered.push_back(C3);
        output.open("output_for_ier_otr" + std::to_string(u) + ".txt");
        for (int i = 0; i < otrez.size(); i++)
        {
            output << otrez[i][0].x << " " << otrez[i][0].y << std::endl << otrez[i][1].x << " " << otrez[i][1].y << std::endl << std::endl;
        }
        output.close();
        output.open("output_for_ier_ser" + std::to_string(u) + ".txt");
        for (int i = 0; i < sered.size(); i++)
        {
            output << sered[i].x << " " << sered[i].y << std::endl;
        }
        output.close();
        zm.clear();
        u++;
        std::vector <Point> new_clast;
        for (int i = 0; i < clast[mass[0]].size(); i++)
        {
            new_clast.push_back(clast[mass[0]][i]);
        }
        for (int i = 0; i < clast[mass[1]].size(); i++)
        {
            new_clast.push_back(clast[mass[1]][i]);
        }
        if (mass[0] < mass[1])
        {
            clast.erase(clast.begin() + mass[0]);
            clast.erase(clast.begin() + mass[1] - 1);
        }
        else
        {
            clast.erase(clast.begin() + mass[0]);
            clast.erase(clast.begin() + mass[1]);
        }
        clast.push_back(new_clast);
        std::vector <Point> clast_sr;
        for (int i = 0; i < clast.size(); i++)
        {
            Point Tochk;
            Tochk.x = 0;
            Tochk.y = 0;
            for (int j = 0; j < clast[i].size(); j++)
            {
                Tochk.x += clast[i][j].x;
                Tochk.y += clast[i][j].y;
            }
            Tochk.x /= clast[i].size();
            Tochk.y /= clast[i].size();
            clast_sr.push_back(Tochk);
        }
        rasst.clear();
        for (int i = 0; i < clast.size(); i++)
        {
            for (int j = 0; j < i; j++)
            {
                rasst[j].push_back(dist(clast_sr[i], clast_sr[j]));
            }
            std::vector <float> vspom;
            for (int j = 0; j < i; j++)
            {
                vspom.push_back(dist(clast_sr[i], clast_sr[j]));
            }
            vspom.push_back(0);
            rasst.push_back(vspom);
            vspom.clear();
        }
        clast_sr.clear();
        new_clast.clear();
    }
}