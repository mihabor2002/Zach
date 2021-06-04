#pragma once
#include <vector>
#include <fstream>
#include "field.h"
#include "sp_tree.h"
#include <cfloat>
using namespace std;
//211-Borovikov-Mikhail-2021
class create_sp
{
public:
	int rfind(vector <Point*> p1, Point p)
	{
		for (int i = 0; i < p1.size(); i++)
		{
			if ((*p1[i]).x == p.x && (*p1[i]).y==p.y)
			{
				return 0;
			}
		}
		return 1;
	}
	void funk(Field* TheField)//funlcia sozdaet gistogrammu rasstoyanii polya
	{
		int k, n, id_p, id_q;
		double min_dist;
		n = (*TheField).numb_Point;
		Sp_Tree tree = Sp_Tree(&((*TheField).Pointarr[0]));
		for (int i = 0; i <= n - 2; i++)
		{
			id_p = -1;
			id_q = -1;
			min_dist = DBL_MAX;
			for (int j = 0; j < tree.t_p.size(); j++)
				for (k = 0; k < n; k++)
				{
					if ((k != j) &&
						(rfind(tree.t_p,(*TheField).Pointarr[k]) == 1) &&
						((*TheField).MatrD.matrx[tree.t_p[j]->ID][k] < min_dist))
					{
						id_p = k;
						id_q = tree.t_p[j]->ID;
						min_dist = (*TheField).MatrD.matrx[tree.t_p[j]->ID][k];
					}
				}
			tree.addPoint(&((*TheField).Pointarr[id_p]), &((*TheField).Pointarr[id_q]));
			tree.t_p.push_back(&(*TheField).Pointarr[id_p]);
			tree.edges.push_back((*TheField).MatrD.matrx[id_p][id_q]);
		}
		std::ofstream tree_file;
		tree_file.open("SpaningTree.txt");
		tree.displayTree(tree_file);
		tree_file.close();
		tree_file.open("GistField.txt");
		for (int i = 0; i < tree.edges.size(); i++)
		{
			tree_file << tree.edges[i] << std::endl;
		}
		tree_file.close();
	}
};