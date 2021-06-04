#pragma once
#include <vector>
#include "point.h"
//211-Borovikov-Mikhail-2021
class Claster
{
public:
    int numb_Point;//kolichestvo tochek v klastere
    std::vector <Point> arr;//vektor s tochkami etogo klastera
    int t_nach;
    int t_kon;
};

