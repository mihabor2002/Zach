#pragma once
#include <vector>
#include "point.h"
//211-Borovikov-Mikhail-2021
class Clowd
{
public:
	int length;//kolichestvo tochek v oblake
	Clowd()
	{
		length = 0;
	}
	std::vector <Point> arr;//vector tochek v oblake
};

