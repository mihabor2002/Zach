#pragma once
#include <vector>
#include "b_matr.h"
#include "claster.h"
//211-Borovikov-Mikhail-2021
class Find_Clasters
{
public:
	B_Distant_matrix MatrB;//binarnaya matritsa dannogo ob'ekta Find_clasters
	std::vector <int> dbpoint;//esli zapuschen DBSCAN, to dbpoint[i]=j, gde i-aya tochka v matritse MatrB imeet nomer j v polye
	std::vector <Claster> Clasterarr;//massiv claserov
	char tip_alg;//tip algoritma,dlya kotorogo ispolzuetsya dannii obekt
	Find_Clasters()
	{

	}
};

