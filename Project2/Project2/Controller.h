#pragma once
#include "field.h"
#include "matrx.h"
#include "binary.h"
#include "wave.h"
#include "dbbinary.h"
#include "dbwave.h"
#include "k_sredn.h"
#include "k_sredn_yadr.h"
#include "em_alg.h"
#include "forel.h"
#include "ierarch.h"
#include "create_sp.h"
#include "Interpolation.h"
//211-Borovikov-Mikhail-2021
#include <vector>
#include <iostream>
#include <string>
float sum_200(void);
class Controller
{
public:
    std::string filename_log;//nazvanie faila loga controllera
    std::ofstream log;//fail s logom
    Field TheField;//polye kontrollera
    Controller()
    {
        filename_log = "log_cntrl.txt";
        log.open(filename_log);
        Field TheField;
    }
    Controller(const Controller& c) { filename_log = c.filename_log; TheField = c.TheField; }

    const Controller& operator=(const Controller& c)
    {
        filename_log = c.filename_log;
        TheField = c.TheField;
        return *this;
    }
    void generate_clowd(float mX, float mY, float sX, float sY, int n)//finkciya generazii oblaka v polye
    {
        log << "poluchen zapros na vipolnenie programmi generazii polya s mX=" << mX << ", mY = " << mY << ", sX = " << sX << ", sY = " << sY << ", n = " << n << std::endl;
        Clowd TheClowd;
        for (int i = 0; i < n; i++)
        {
            Point ThePoint;
            ThePoint.x = mX + (sum_200() / 200) * sX;
            ThePoint.y = mY + (sum_200() / 200) * sY;
            ThePoint.ID = TheField.numb_Point;
            (TheClowd.arr).push_back(ThePoint);
            TheField.Pointarr.push_back(ThePoint);
            TheField.numb_Point++;
            TheClowd.length++;
        }
        (TheField.Clowdarr).push_back(TheClowd);
        TheField.numb_Clowd++;
        log << "programma generazii polya vipolnena" << std::endl;
    }
    void print_field()//pechat polya v fail "output_clowds.txt"
    {
        log << "poluchen zapros na vipolnenie programmi pechati polya" << std::endl;
        int i, j;
        std::ofstream output;
        output.open("output_clowds.txt");
        for (i = 0; i < TheField.numb_Clowd; i++)
        {
            for (j = 0; j < TheField.Clowdarr[i].length; j++)
            {
                output << ((TheField.Clowdarr[i].arr)[j]).x << " " << ((TheField.Clowdarr[i].arr)[j]).y << " " << i + 1 << std::endl;
            }
        }
        output.close();
        log << "programma pechati polya vipolnena" << std::endl;
    }
    void change_sost()//funkciya menyaet sostoyanie polya
    {
        log << "poluchen zapros na vipolnenie programmi izmeneniya sostoyaniya polya" << std::endl;
        TheField.change_sost();
        log << "programma izmeneniya sostoyaniya polya vipolnena" << std::endl;
    }
    bool get_sost()//funkcia vozvraschaet sostoyanie polya
    {
        log << "poluchen zapros na tekuschee sostoyanie polya" << std::endl;
        return TheField.sostoyanie_polya;
        log << "otvet dan" << std::endl;
    }
    void Matrx()//funkcia zapuaskaet funkciuy Matrx dlya polya
    {
        log << "poluchen zapros na generaziu binarnoy martistsi" << std::endl;
        matrx Thematrx;
        Thematrx.funk(&TheField);
    }
    void Binary(float alph)
    {
        log << "poluchen zapros na zapusk komandi binary" << std::endl;
        binary Thebinary;
        Thebinary.funk(&TheField, alph);
        log << "komanda binary vipolnena" << std::endl;
    }
    void Voronoy()
    {
        log << "poluchen zapros na zapusk komandi voronoy" << std::endl;
        TheField.Voronoy();
        log << "komanda vipolnena";
    }
    void Wave()
    {
        log << " poluchen zapros na zapusk komandi wave";
        wave Thewave;
        Thewave.funk(&TheField);
        log << "komanda wave vipolnena";
    }
    void Dbbinary(float alph, int k)
    {
        log << "poluchen zapros na zapusk komandi dbbinary";
        dbbinary Thedbbinary;
        Thedbbinary.funk(&TheField, alph, k);
        log << "komanda dbbinary vipolnena";
    }
    void Dbwave()
    {
        log << "poluchen zapros na zapusk komandi dbwave";
        dbwave Thedbwave;
        Thedbwave.funk(&TheField);
        log << "komanda dbwave vipolnena";
    }
    void K_sredn(int k)
    {
        log << "poluchen zapros na zapusk komandi k_sredn";
        k_sredn Thek_sredn;
        Thek_sredn.funk(&TheField, k);
        log << "komanda k_sredn vipolnena";
    }
    void EM_alg(int k)
    {
        log << "poluchen zapros na zapusk komandi EM-algoritma";
        em_alg Theem_alg;
        Theem_alg.funk(&TheField, k);
        log << "komanda EM_alg vipolnena";
    }
    void K_sredn_yadr(int k, int p)
    {
        log << "poluchen zapros na zapusk komandi k_sredn_yadr s k = " << k << ", p=" << p << std::endl;
        k_sredn_yadr Thek_sredn_yadr;
        Thek_sredn_yadr.funk(&TheField, k, p);
        log << "komanda k_sredn_yadr vipolnena" << std::endl;
    }
    void print_clasters1()
    {
        log << "poluchen zapros na zapusk komandi print_clasters1";
        std::ofstream output_for_clast;
        std::ofstream output_for_point;
        std::ofstream output_gnu;
        for (int u = 0; u < TheField.Find_Clastersarr.size(); u++)
        {
            output_for_clast.open("output_claster_exp" + std::to_string(u) + ".txt");
            output_gnu.open("output_gnu_exp" + std::to_string(u) + ".txt");
            if (TheField.Find_Clastersarr[u].tip_alg == 'b')
            {
                for (int i = 0; i < TheField.Find_Clastersarr[u].Clasterarr.size(); i++)
                {
                    for (int j = TheField.Find_Clastersarr[u].Clasterarr[i].t_nach; j < TheField.Find_Clastersarr[u].Clasterarr[i].t_kon + 1; j++)
                    {
                        for (int v = 0; v < TheField.Find_Clastersarr[u].Clasterarr[i].arr.size(); v++)
                        {
                            if (TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].a == j)
                            {
                                for (int w = 0; w < TheField.Find_Clastersarr[u].Clasterarr[i].arr.size(); w++)
                                {
                                    if (TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].b == j + 1 && TheField.Find_Clastersarr[u].MatrB.matrx[TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].ID][TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].ID] == 1)
                                    {
                                        output_for_clast << TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].x << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].y << " " << 1 << std::endl;
                                        output_for_clast << TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].x << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].y << " " << 1 << std::endl;
                                        output_for_clast << std::endl << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if (TheField.Find_Clastersarr[u].tip_alg == 'k')
            {
                for (int i = 0; i < TheField.Find_Clastersarr[u].Clasterarr.size(); i++)
                {
                    for (int j = 0; j < TheField.Find_Clastersarr[u].Clasterarr[i].arr.size() - 1; j++)
                    {
                        output_for_clast << TheField.Find_Clastersarr[u].Clasterarr[i].arr[TheField.Find_Clastersarr[u].Clasterarr[i].arr.size() - 1].x << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[TheField.Find_Clastersarr[u].Clasterarr[i].arr.size() - 1].y << " " << 1 << std::endl;
                        output_for_clast << TheField.Find_Clastersarr[u].Clasterarr[i].arr[j].x << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[j].y << " " << 1 << std::endl;
                        output_for_clast << std::endl << std::endl;
                    }
                }
            }
            else
            {
                std::vector <int> numb;
                for (int i = 0; i < TheField.numb_Point; i++)
                {
                    numb.push_back(i);
                }
                for (int i = 0; i < TheField.Find_Clastersarr[u].dbpoint.size(); i++)
                {
                    std::cout << TheField.Find_Clastersarr[u].dbpoint[i] << " ";
                }
                for (int i = 0; i < TheField.Find_Clastersarr[u].Clasterarr.size(); i++)
                {
                    for (int j = TheField.Find_Clastersarr[u].Clasterarr[i].t_nach; j < TheField.Find_Clastersarr[u].Clasterarr[i].t_kon + 1; j++)
                    {
                        for (int v = 0; v < TheField.Find_Clastersarr[u].Clasterarr[i].arr.size(); v++)
                        {
                            if (TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].a == j)
                            {
                                for (int w = 0; w < TheField.Find_Clastersarr[u].Clasterarr[i].arr.size(); w++)
                                {
                                    std::cout << TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].ID << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].ID << std::endl;
                                    if (TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].b == j + 1 && TheField.Find_Clastersarr[u].MatrB.matrx[vspom(TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].ID, TheField.Find_Clastersarr[u].dbpoint)][vspom(TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].ID, TheField.Find_Clastersarr[u].dbpoint)] == 1)
                                    {
                                        output_for_clast << TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].x << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].y << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].DB << std::endl;
                                        output_for_clast << TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].x << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].y << " " << TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].DB << std::endl;
                                        output_for_clast << std::endl << std::endl;
                                        numb[TheField.Find_Clastersarr[u].Clasterarr[i].arr[w].ID] = -1;
                                        numb[TheField.Find_Clastersarr[u].Clasterarr[i].arr[v].ID] = -1;
                                    }
                                }
                            }
                        }
                    }
                }
                for (int i = 0; i < TheField.numb_Point; i++)
                {
                    if (numb[i] != -1)
                    {
                        output_for_clast << TheField.Pointarr[numb[i]].x << " " << TheField.Pointarr[numb[i]].y << " " << 0 << std::endl;
                    }
                }
            }
            output_for_clast.close();
            output_gnu << "set multiplot" << std::endl;
            output_gnu << "set yrange[-20:20]" << std::endl << "set xrange [-20:20]" << std::endl;
            output_gnu << "set palette defined (0 \"brown\",1 \"blue\", 2 \"green\",3 \"red\", 4 \"orange\")" << std::endl;
            output_gnu << "plot \"output_claster_exp" << u << ".txt\"w l palette" << ", ";
            output_gnu << "\"output_claster_exp" << u << ".txt\" w p palette" << std::endl;
            output_gnu.close();
        }
        log << "komanda print_clasters1 vipolnena";
    }
    void Forel(float R)
    {
        log << "poluchen zapros na vipolnenie programmi algoritma forel";
        forel Theforel;
        Theforel.funk(&TheField, R);
        log << "programma algoritma forel vipolnena";
    }
    void Ierarch()
    {
        log << "poluchen zapros na vipolnenie programmi ierarchicheskogo algoritma";
        ierarch Theierarch;
        Theierarch.funk(&TheField);
        log << "programma ierarchicheskogo algoritma vipolnena";
    }
    void gist()
    {
        create_sp Thecr;
        Thecr.funk(&TheField);
    }
    void Triang()
    {
        Delaunay_triangulation T;
        vector <Point> arr;
        for (int i = 0; i < TheField.numb_Point; i++)
        {
            arr.push_back(TheField.Pointarr[i]);
        }
        T = TheField.generate_delaunay_trinagulation(arr);
        T.print_Delaunay_triangulation();
    }
    void Inter(Point Z)
    {
        Interpolation TheInt;
        TheInt.Start(Z,&TheField);
    }
    ~Controller()
    {
        log.close();
    }
};

