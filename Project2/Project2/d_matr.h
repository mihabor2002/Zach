#pragma once
#include <vector>
#include <fstream>
//211-Borovikov-Mikhail-2021
class D_Distant_matrix
{
public:
	std::vector <std::vector<float>> matrx;//matritsa rasstoyanii mezhdu tochkami v polye,matrx[i][j]=rasstoyanie mezhdu i-oi i j-oi tochkoi
	int n;//razmer matritsi
	void print_matrxD()//raspechatat' matritsu v fail output_for_martxD.txt
	{
		std::ofstream out;
		out.open("output_for_martxD.txt");
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

