#include "Delaunay_triangulation.h"
//211-Borovikov-Mikhail-2021

Delaunay_triangulation::Delaunay_triangulation() { number_triangles = number_points = 0; }

Delaunay_triangulation::Delaunay_triangulation(const Delaunay_triangulation& Del_tr)
{
	triangles = Del_tr.triangles;
	number_triangles = Del_tr.number_triangles;
	number_points = Del_tr.number_points;
	triangle_indicators = Del_tr.triangle_indicators;
}

const Delaunay_triangulation& Delaunay_triangulation::operator=(const Delaunay_triangulation& Del_tr)
{
	triangles = Del_tr.triangles;
	number_triangles = Del_tr.number_triangles;
	number_points = Del_tr.number_points;
	triangle_indicators = Del_tr.triangle_indicators;
	return *this;
}

void Delaunay_triangulation::add_triangle(Triangle T)
{
	T.assign_id(number_triangles);
	triangles.push_back(T);
	number_triangles++;
	number_points = number_points + 3;
}

Triangle Delaunay_triangulation::get_triangle(int id) { return triangles[id]; }

int Delaunay_triangulation::get_number_triangles() { return number_triangles; }

int Delaunay_triangulation::find_triangle(Point p, Point q)
{
	for (int i = 0; i < number_triangles; i++)
	{
		Point a = triangles[i].get_A(), b = triangles[i].get_B(), c = triangles[i].get_C();
		if (a.ID == p.ID && b.ID == q.ID) return i;
		if (a.ID == p.ID && c.ID == q.ID) return i;
		if (b.ID == p.ID && a.ID == q.ID) return i;
		if (c.ID == p.ID && a.ID == q.ID) return i;
		if (b.ID == p.ID && c.ID == q.ID) return i;
		if (c.ID == p.ID && b.ID == q.ID) return i;
	}
	return -1;
}

void Delaunay_triangulation::create_triangle_indicators()
{
	triangle_indicators.resize(number_triangles);
	for (int i = 0; i < number_triangles; i++) triangle_indicators[i].resize(number_triangles);
	for (int i = 0; i < number_triangles; i++)
	{
		for (int j = 0; j < number_triangles; j++) triangle_indicators[i][j] = 0;
	}
	for (int i = 0; i < number_triangles; i++)
	{
		for (int j = 0; j < number_triangles; j++)
			if (triangles[i].get_neighbouring_triangle_id() == j) triangle_indicators[i][j] = triangle_indicators[j][i] = 1;
	}
}

vector <Point> Delaunay_triangulation::find_neighbouring_points(Point p)
{
	int i, j, k;
    vector <Point> neighbouring_points;

	for (i = 0; i < number_triangles; i++)
	{
		if (triangles[i].get_A().ID == p.ID)
		{
			k = i;
			if (triangles[k].get_A().ID != p.ID) neighbouring_points.push_back(triangles[k].get_A());
			if (triangles[k].get_B().ID != p.ID) neighbouring_points.push_back(triangles[k].get_B());
			if (triangles[k].get_C().ID != p.ID) neighbouring_points.push_back(triangles[k].get_C());
			/*for (j = 0; j < number_triangles; j++)
			{
				if (triangle_indicators[k][j] == 1 && j != k)
				{
					if (triangles[j].get_A().ID == p.ID || triangles[j].get_B().ID == p.ID ||
						triangles[j].get_C().ID == p.ID)
					{
						if (triangles[k].get_A().ID != p.ID) neighbouring_points.push_back(triangles[k].get_A());
						if (triangles[k].get_B().ID != p.ID) neighbouring_points.push_back(triangles[k].get_B());
						if (triangles[k].get_C().ID != p.ID) neighbouring_points.push_back(triangles[k].get_C());
					}
				}
			}*/
		}
		if (triangles[i].get_B().ID == p.ID)
		{
			k = i;
			if (triangles[k].get_A().ID != p.ID) neighbouring_points.push_back(triangles[k].get_A());
			if (triangles[k].get_B().ID != p.ID) neighbouring_points.push_back(triangles[k].get_B());
			if (triangles[k].get_C().ID != p.ID) neighbouring_points.push_back(triangles[k].get_C());
			for (j = 0; j < number_triangles; j++)
			{
				if (triangle_indicators[k][j] == 1 && j != k)
				{
					if (triangles[j].get_A().ID == p.ID || triangles[j].get_B().ID == p.ID ||
						triangles[j].get_C().ID == p.ID)
					{
						if (triangles[k].get_A().ID != p.ID) neighbouring_points.push_back(triangles[k].get_A());
						if (triangles[k].get_B().ID != p.ID) neighbouring_points.push_back(triangles[k].get_B());
						if (triangles[k].get_C().ID != p.ID) neighbouring_points.push_back(triangles[k].get_C());
					}
				}
			}
		}
		if (triangles[i].get_C().ID == p.ID)
		{
			k = i;
			if (triangles[k].get_A().ID != p.ID) neighbouring_points.push_back(triangles[k].get_A());
			if (triangles[k].get_B().ID != p.ID) neighbouring_points.push_back(triangles[k].get_B());
			if (triangles[k].get_C().ID != p.ID) neighbouring_points.push_back(triangles[k].get_C());
			/*for (j = 0; j < number_triangles; j++)
			{
				if (triangle_indicators[k][j] == 1 && j != k)
				{
					if (triangles[j].get_A().ID == p.ID || triangles[j].get_B().ID == p.ID ||
						triangles[j].get_C().ID == p.ID)
					{
						if (triangles[k].get_A().ID != p.ID) neighbouring_points.push_back(triangles[k].get_A());
						if (triangles[k].get_B().ID != p.ID) neighbouring_points.push_back(triangles[k].get_B());
						if (triangles[k].get_C().ID != p.ID) neighbouring_points.push_back(triangles[k].get_C());
					}
				}
			}*/
		}
	}
	for (i = 0; i < neighbouring_points.size(); i++)
	{
		for (j = i + 1; j < neighbouring_points.size(); j++)
		{
			if (neighbouring_points[i].ID == neighbouring_points[j].ID)
			{
				neighbouring_points.erase(neighbouring_points.begin() + j);
			}
		}
	}
	return neighbouring_points;
}

void Delaunay_triangulation::print_Delaunay_triangulation()
{
	ofstream f("Delaunay_triangulation.txt"),
		g("Delaunay_triangulation.plt");
	g << "set size square\n" << "plot 'Delaunay_triangulation.txt' with linespoints\n";
	for (int i = 0; i < number_triangles; i++)
	{
		f << triangles[i].get_A().x << " " << triangles[i].get_A().y << std::endl;
		f << triangles[i].get_B().x << " " << triangles[i].get_B().y << std::endl;
		f << triangles[i].get_C().x << " " << triangles[i].get_C().y << std::endl;
		f << triangles[i].get_A().x << " " << triangles[i].get_A().y << std::endl;
		f << std::endl << std::endl;
	}
	f.close();
	g.close();
}

void Delaunay_triangulation::print(int n)
{

}

void Delaunay_triangulation::print_gif(int n)
{

}