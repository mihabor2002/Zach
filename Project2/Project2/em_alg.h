#pragma once
#include "field.h"
#include <vector>
#include <string>
//211-Borovikov-Mikhail-2021
float opred(std::vector<std::vector<float>>& matrx);
std::vector<float> sobst_znach(std::vector<std::vector<float>>& matrx);
std::vector<std::vector<float>> sobstv_vector(std::vector<std::vector<float>>& matrx);
class em_alg
{public:

    void funk(Field* TheField, int k);
};
