#include "forel.h"
//211-Borovikov-Mikhail-2021
void forel::funk(Field* ThisField, float R)
{
    std::ofstream output1, output2, output3, output;
    output.open("output_for_pointf.txt");
    for (int i = 0; i < (*ThisField).numb_Point; i++)
    {
        output << (*ThisField).Pointarr[i].x << " " << (*ThisField).Pointarr[i].y << std::endl;
    }
    output.close();
    Field TheField;
    for (int i = 0; i < (*ThisField).numb_Point; i++)
    {
        TheField.Pointarr.push_back((*ThisField).Pointarr[i]);
        TheField.numb_Point++;
    }
    int j = 0;
    std::vector<Point> deleted;
    std::vector <float> tokr;
    Point C;
    C.x = 100;
    C.y = 100;
    std::vector <std::vector<float>> okr;
    deleted.push_back(C);
    while (TheField.numb_Point != 0)
    {
        int u = 1;
        Point Tek_toch, Centr_toch;
        Tek_toch = TheField.Pointarr[0];
        Centr_toch.x = 0;
        Centr_toch.y = 0;
        std::vector <Point> sosedn;
        for (int i = 0; i < TheField.numb_Point; i++)
        {
            if (dist(Tek_toch, TheField.Pointarr[i]) < R)
            {
                sosedn.push_back(TheField.Pointarr[i]);
            }
        }
        for (int i = 0; i < sosedn.size(); i++)
        {
            Centr_toch.x += sosedn[i].x;
            Centr_toch.y += sosedn[i].y;
        }
        std::cout << R << std::endl;
        Centr_toch.x /= sosedn.size();
        Centr_toch.y /= sosedn.size();
        std::cout << Centr_toch.x << " " << Centr_toch.y << " " << Tek_toch.x << " " << Tek_toch.y << std::endl;
        j++;
        std::cout << j << " " << R << std::endl;
        output1.open("foutput_for_for" + std::to_string(j) + ".txt");
        for (int i = 0; i < okr.size(); i++)
        {
            for (int l = 0; l < 4; l++)
            {
                output1 << okr[i][l] << " ";
            }
            output1 << std::endl;
        }
        output1 << Centr_toch.x << " " << Centr_toch.y << " " << 2 * R << " " << 2 * R;
        output1.close();
        output2.open("foutput_for_for_sos" + std::to_string(j) + ".txt");
        for (int i = 0; i < sosedn.size(); i++)
        {
            output2 << sosedn[i].x << " " << sosedn[i].y << std::endl;
        }
        output2.close();
        output3.open("foutput_for_for_ud" + std::to_string(j) + ".txt");
        for (int i = 0; i < deleted.size(); i++)
        {
            output3 << deleted[i].x << " " << deleted[i].y << std::endl;
        }
        output3.close();
        while (!P_eq(Centr_toch, Tek_toch))
        {
            std::ofstream output;
            Tek_toch = Centr_toch;
            Centr_toch.x = 0;
            Centr_toch.y = 0;
            sosedn.clear();
            for (int i = 0; i < TheField.numb_Point; i++)
            {
                if (dist(Tek_toch, TheField.Pointarr[i]) < R)
                {
                    sosedn.push_back(TheField.Pointarr[i]);
                }
            }
            for (int i = 0; i < sosedn.size(); i++)
            {
                Centr_toch.x += sosedn[i].x;
                Centr_toch.y += sosedn[i].y;
            }
            Centr_toch.x /= +sosedn.size();
            Centr_toch.y /= sosedn.size();
            j++;
            std::cout << j << " " << R << std::endl;
            output1.open("foutput_for_for" + std::to_string(j) + ".txt");
            for (int i = 0; i < okr.size(); i++)
            {
                for (int l = 0; l < 4; l++)
                {
                    output1 << okr[i][l] << " ";
                }
                output1 << std::endl;
            }
            output1 << Centr_toch.x << " " << Centr_toch.y << " " << 2 * R << " " << 2 * R;
            output1.close();
            output2.open("foutput_for_for_sos" + std::to_string(j) + ".txt");
            for (int i = 0; i < sosedn.size(); i++)
            {
                output2 << sosedn[i].x << " " << sosedn[i].y << std::endl;
            }
            output2.close();
            output3.open("foutput_for_for_ud" + std::to_string(j) + ".txt");
            for (int i = 0; i < deleted.size(); i++)
            {
                output3 << deleted[i].x << " " << deleted[i].y << std::endl;
            }
            std::cout << deleted.size() << "kek" << std::endl;
            output3.close();
        }
        for (int i = 0; i < sosedn.size(); i++)
        {
            TheField.Pointarr.erase(TheField.Pointarr.begin() + sosedn[i].ID);
            for (int j = 0; j < sosedn.size(); j++)
            {
                if (sosedn[j].ID > sosedn[i].ID)
                {
                    sosedn[j].ID--;
                }
            }
            TheField.numb_Point--;
            deleted.push_back(sosedn[i]);
        }
        for (int i = 0; i < TheField.numb_Point; i++)
        {
            TheField.Pointarr[i].ID = i;
        }
        tokr.clear();
        tokr.push_back(Centr_toch.x);
        tokr.push_back(Centr_toch.y);
        tokr.push_back(2 * R);
        tokr.push_back(2 * R);
        okr.push_back(tokr);
    }
}