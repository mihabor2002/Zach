#pragma once
#include "field.h"
#include <vector>
using namespace std;
class Interpolation
{
public:
	int Start(Point& P, Field *TheField);
	float Evaluation(vector <Point>& All_Point);
};

