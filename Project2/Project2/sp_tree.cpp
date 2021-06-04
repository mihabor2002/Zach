#include "sp_tree.h"
#include "field.h"
//211-Borovikov-Mikhail-2021
Sp_Tree::Sp_Tree()
{
	neighbors = {};
	t_p = {};
	edges = {};
}

Sp_Tree::Sp_Tree(Point *P)
{
	neighbors = {};
	t_p = {};
	edges = {};
	t_p.push_back(P);
}

void Sp_Tree::addPoint(Point* p, Point* q)
{
	if (point == q)
	{
		Sp_Tree* new_tree = new Sp_Tree(p);
		neighbors.push_back(new_tree);
	}
	else
		for (Sp_Tree* nb : neighbors)
			nb->addPoint(p, q);
}

void Sp_Tree::displayTree(std::ofstream& out)
{
	for (Sp_Tree* nb : neighbors)
	{
		point->printToFile(out, 0);
		nb->point->printToFile(out, 0);
		out << std::endl;
		nb->displayTree(out);
	}
}