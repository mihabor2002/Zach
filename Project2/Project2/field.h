#pragma once
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "point.h"
#include "claster.h"
#include "b_matr.h"
#include "d_matr.h"
#include "clowd.h"
#include "find_cl.h"
#include "Delaunay_triangulation.h"
using namespace std;
//211-Borovikov-Mikhail-2021
int vspom(int k, std::vector <int>& v);
class Field
{
public:
    int numb_Point;//kolichestvo tochek v polye
    int numb_Clowd;//kolichestvo oblakov v pole
    bool sostoyanie_polya;//=0, esli polye v sostoyanii analiza, 1 v obratnom sluchae
    Field()
    {
        sostoyanie_polya = 1;
        numb_Point = 0;
        numb_Clowd = 0;
    }
    void change_sost()//funkcia izmeneniya sostoyaniya polya
    {
        if (sostoyanie_polya > 0)
        {
            sostoyanie_polya = 0;
        }
        else
        {
            sostoyanie_polya = 1;
        }
    }
    std::vector <Clowd> Clowdarr;//spisok oblakov polya
    std::vector <Point> Pointarr;//spisok tochek polya
    std::vector <Find_Clasters> Find_Clastersarr;//spisok Find_Clasterov
    int last_Find_Clasters_bin = 0;
    int last_Find_Clasters_db = 0;
    D_Distant_matrix MatrD;//martitsa rasstoyanii mezhdu tochkami polya
    void find_min_y(vector <Point>& points);
    void polar_angle_sort(vector <Point>& points);
    float left_rotate(Point p, Point q, Point r);
    vector <Point> convex_hull(vector <Point> points);
    Delaunay_triangulation generate_delaunay_trinagulation(vector <Point>& points);
    void Voronoy()
    {
        std::vector <Line> sp;
        std::vector <Point> st2;
        std::vector <Point> st3;
        std::vector <std::vector<Point>> st4;
        Line L1;
        for (int i = 0; i < numb_Point; i++)
        {
            for (int j = 0; j < numb_Point; j++)
            {
                if (i != j)
                {
                    L1 = serper(Pointarr[i], Pointarr[j]);
                    sp.push_back(L1);
                }
            }
            for (int j = 0; j < numb_Point; j++)
            {
                for (int v = 0; v < numb_Point; v++)
                {
                    if (v != j)
                    {
                        if (kolvo_peresech(sp[j], sp[v]) == 1)
                        {
                            st2.push_back(peresech(sp[j], sp[v]));
                        }
                    }
                }
            }
            for (int j = 0; j < st2.size(); j++)
            {
                int T = 0;
                for (int v = 0; v < numb_Point - 1; v++)
                {
                    if (!po_odnu_storonu(st2[j], Pointarr[i], sp[v]))
                    {
                        T = 1;
                        break;
                    }
                }
                if (T == 0)
                {
                    st3.push_back(st2[j]);
                }
            }
            st4.push_back(st3);
            sp.clear();
            st2.clear();
        }
        std::ofstream output;
        output.open("output_vor.txt");
        for (int i = 0; i < numb_Point - 1; i++)
        {
            for (int j = 0; j < st4[i].size() - 1; j++)
            {
                output << st4[i][j].x << " " << st4[i][j].y << std::endl;
                output << st4[i][j + 1].x << " " << st4[i][j + 1].y << std::endl << std::endl << std::endl;
            }
            output << st4[i][st4[i].size() - 1].x << " " << st4[i][st4[i].size() - 1].y << std::endl;
            output << st4[i][0].x << " " << st4[i][0].y << std::endl << std::endl << std::endl;
        }
    }
};

