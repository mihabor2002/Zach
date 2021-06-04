#pragma once
#include <vector>
#include "point.h"
#include <fstream>
//211-Borovikov-Mikhail-2021
class Sp_Tree
{
private:
	Point* point;
public:
	std::vector<Sp_Tree*> neighbors;
	std::vector<Point*> t_p;
	std::vector<double> edges;
	Sp_Tree();
	Sp_Tree(Point* p);
	void addPoint(Point* p, Point* q);
	void displayTree(std::ofstream& out);
};

