#include "Interpolation.h"
#define NNNN 10
#include "Delaunay_triangulation.h"
#include "k_sredn.h"

int Interpolation::Start(Point& P, Field *TheField)
{
    
	/*k_sredn A_K_M;
    A_K_M.funk(TheField,NNNN);
	Find_Clasters Clus = (*TheField).Find_Clastersarr[0];
	int I = 0;
	Point PC = Clus.Clasterarr[0].arr[Clus.Clasterarr[0].arr.size()-1];
	float D = dist(P, PC);
	for (int i = 0; i < NNNN; i++)
	{
		Point PC = Clus.Clasterarr[i].arr[Clus.Clasterarr[i].arr.size() - 1];
		if (D > dist(P, PC))
		{
			D = dist(P, PC);
			I = i;
		}
	}
	vector <Point> cluster = Clus.Clasterarr[I].arr;
	Delaunay_triangulation Triang = (*TheField).generate_delaunay_trinagulation(cluster);

	for (int i = 0; i < Triang.get_number_triangles(); i++)
	{
		if (Triang.get_triangle(i).Belongs_P(P) == 1)
		{
			float S, v;
			vector <float> p;
			S = v = 0;
			cout << "P=( " << P.x << ", " << P.y << ")" << " F(P)=" << sin((P.x * P.x + P.y * P.y)) << endl;
            Point PP(Triang.get_triangle(i).get_A());
			p.push_back((1 / (10 * dist(PP, P))) * (1 / (10 * dist(PP, P))));
			S += PP.get_function_value() * p[p.size() - 1];
			v += p[p.size() - 1];
			Point PPP = Triang.get_triangle(i).get_B();
			p.push_back((1 / (10 * dist(PPP, P))) * (1 / (10 * dist(PPP, P))));
			S += PPP.get_function_value() * p[p.size() - 1];
			v += p[p.size() - 1];
			Point PPPP = Triang.get_triangle(i).get_C();
			p.push_back((1 / (10 * dist(PPPP, P))) * (1 / (10 * dist(PPPP, P))));
			S += PPPP.get_function_value() * p[p.size() - 1];
			v += p[p.size() - 1];
			cout << "S=" << S << "v=" << v << endl;
			S = S / v;
			float xx = P.x * P.x;
			float yy = P.y * P.y;
			float s = (xx + yy);
			cout << "Oshibka = " << S - sin(s) << endl;

			ofstream g, h;
			g.open("Interpolation.txt");
			//g.open("C:\\My_Program\\INT\\Algoritms\\K_means_Ker.txt");
			h.open("Interpolation.plt");
			h << "set pm3d" << endl;
			h << "splot[-1.5:1.5][-1.5:1.5][-1.2:1.2] sin(x**2 + y**2), 'Interpolation.txt' " << endl;
			g << P.x << " " << P.y << " " << S << endl;
			Interpolation A_Interpolation;
			A_Interpolation.Evaluation(cluster);
			return 1;
		}
	}
	*/
	return -1;
}
