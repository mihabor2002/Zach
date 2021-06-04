#pragma once
#include <vector>
#include <fstream>
//211-Borovikov-Mikhail-2021
class B_Distant_matrix
{
public:
	std::vector <std::vector <bool>> matrx;//matritsa takaya, chto matrx[i][j] - 1,esli rasstoyanie mezhdu i-oi i j-oi tochkoi < alpha, i 0 v obratnom sluchae 
	int n;//razmer matritsi
	void print_matrxB()//pechat' matritsi v fail "output_for_martxB.txt" 
	{
		std::ofstream out;
		out.open("output_for_martxB.txt");
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				out << matrx[i][j] << " ";
			}
			out << std::endl;
		}
		out.close();
	}
};

