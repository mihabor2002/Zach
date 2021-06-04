#include "em_alg.h"
//211-Borovikov-Mikhail-2021
float opred(std::vector<std::vector<float>>& matrx)//funkcia shitaet opredelitel matritsi 2*2
{
	return matrx[0][0] * matrx[1][1] - matrx[0][1] * matrx[1][0];
}
std::vector<float> sobst_znach(std::vector<std::vector<float>>& matrx)//funkcia shitaet sobstvennie znacheniya matritsi 2*2
{
	float descr = (matrx[0][0] + matrx[1][1]) * (matrx[0][0] + matrx[1][1]) - 4 * opred(matrx);
	std::vector<float> korni;
	if (descr >= 0)
	{
		korni.push_back(((matrx[0][0] + matrx[1][1]) + sqrt(descr)) / 2);
		korni.push_back(((matrx[0][0] + matrx[1][1]) - sqrt(descr)) / 2);
	}
	return korni;
}
std::vector<std::vector<float>> sobstv_vector(std::vector<std::vector<float>>& matrx)//funkcia shitaet sobsvennie vektora matritsi 2*2
{
	std::vector<float> znach = sobst_znach(matrx);
	std::vector<std::vector<float>> vecs;
	if (znach.size() == 0)
	{
		return vecs;
	}
	vecs.push_back(std::vector<float>{matrx[0][0] - znach[0], matrx[1][0]});
	vecs.push_back(std::vector<float>{matrx[0][0] - znach[0], matrx[1][0]});
	vecs.push_back(znach);
	return vecs;
}
void em_alg::funk(Field* ThisField, int k)
{
	std::vector <Point> centr;//centralnie tochki 
	for (int i = 0; i < k; i++)
	{
		centr.push_back((*ThisField).Pointarr[i * ((*ThisField).numb_Point / k)]);
	}
	for (int i = 0; i < (*ThisField).numb_Point; i++)
	{
		for (int c = 0; c < k; c++)
		{
			(*ThisField).Pointarr[i].c_po_x.push_back(dist((*ThisField).Pointarr[i], centr[c]) / (dist((*ThisField).Pointarr[i], centr[0]) + dist((*ThisField).Pointarr[i], centr[1])));
			(*ThisField).Pointarr[i].x_po_c.push_back(dist((*ThisField).Pointarr[i], centr[c]) / (dist((*ThisField).Pointarr[i], centr[0]) + dist((*ThisField).Pointarr[i], centr[1])));
		}
	}
	std::ofstream outpu;
	outpu.open("output_for_em-1.txt");
	for (int i = 0; i < (*ThisField).numb_Point; i++)
	{
		outpu << (*ThisField).Pointarr[i].x << " " << (*ThisField).Pointarr[i].y << " ";
		float s = (*ThisField).Pointarr[i].c_po_x[0];
		int t = 0;
		for (int c = 0; c < k; c++)
		{
			if ((*ThisField).Pointarr[i].c_po_x[c] > s)
			{
				s = (*ThisField).Pointarr[i].c_po_x[c];
				t = c;
			}
		}
		outpu << t << std::endl;
	}
	outpu.close();
	for (int q = 0; q < 300; q++)
	{
		for (int i = 0; i < k; i++)
		{
			centr[i].P_c = 0;
			for (int j = 0; j < (*ThisField).numb_Point; j++)
			{
				centr[i].P_c += (*ThisField).Pointarr[j].c_po_x[i];
			}
			centr[i].P_c = float(centr[i].P_c) / (*ThisField).numb_Point;
		}
		std::vector <std::vector <float>> mu(k);
		for (int c = 0; c < k; c++)
		{
			mu[c].push_back(0);
			mu[c].push_back(0);
		}
		for (int c = 0; c < k; c++)
		{
			mu[c][0] = 0;
			mu[c][1] = 0;
			for (int i = 0; i < (*ThisField).numb_Point; i++)
			{
				mu[c][0] += (*ThisField).Pointarr[i].c_po_x[c] / float((*ThisField).numb_Point) / centr[c].P_c * (*ThisField).Pointarr[i].x;
				mu[c][1] += (*ThisField).Pointarr[i].c_po_x[c] / float((*ThisField).numb_Point) / centr[c].P_c * (*ThisField).Pointarr[i].y;
			}
		}
		std::vector <std::vector<std::vector <float>>> eps(k);
		for (int c = 0; c < k; c++)
		{
			std::vector <float> w;
			w.push_back(0);
			w.push_back(0);
			eps[c].push_back(w);
			eps[c].push_back(w);
		}
		for (int c = 0; c < k; c++)
		{
			eps[c][0][0] = 0;
			eps[c][0][1] = 0;
			eps[c][1][0] = 0;
			eps[c][1][1] = 0;
			for (int i = 0; i < (*ThisField).numb_Point; i++)
			{
				eps[c][0][0] += (*ThisField).Pointarr[i].c_po_x[c] / float((*ThisField).numb_Point) / centr[c].P_c * ((*ThisField).Pointarr[i].x - mu[c][0]) * ((*ThisField).Pointarr[i].x - mu[c][0]);
				eps[c][0][1] += (*ThisField).Pointarr[i].c_po_x[c] / float((*ThisField).numb_Point) / centr[c].P_c * ((*ThisField).Pointarr[i].x - mu[c][0]) * ((*ThisField).Pointarr[i].y - mu[c][1]);
				eps[c][1][1] += (*ThisField).Pointarr[i].c_po_x[c] / float((*ThisField).numb_Point) / centr[c].P_c * ((*ThisField).Pointarr[i].y - mu[c][1]) * ((*ThisField).Pointarr[i].y - mu[c][1]);
			}
			eps[c][1][0] = eps[c][0][1];
		}
		for (int i = 0; i < (*ThisField).numb_Point; i++)
		{
			float ch;
			for (int c = 0; c < k; c++)
			{
				ch = (((*ThisField).Pointarr[i].x - mu[c][0]) * ((*ThisField).Pointarr[i].x - mu[c][0]) * eps[c][1][1] - 2 * ((*ThisField).Pointarr[i].x - mu[c][0]) * ((*ThisField).Pointarr[i].y - mu[c][1]) * (eps[c][1][0]) + ((*ThisField).Pointarr[i].y - mu[c][1]) * ((*ThisField).Pointarr[i].y - mu[c][1]) * (eps[c][0][0])) / ((eps[c][0][0] * eps[c][1][1] - eps[c][1][0] * eps[c][1][0]));
				ch = exp(ch * (-0.5));
				ch = ch / sqrt(2) / sqrt(2 * asin(1)) / (eps[c][0][0] * eps[c][1][1] - eps[c][1][0] * eps[c][1][0]);
				(*ThisField).Pointarr[i].x_po_c[c] = ch;
				ch = 0;
			}
		}
		for (int i = 0; i < (*ThisField).numb_Point; i++)
		{
			for (int c = 0; c < k; c++)
			{
				float ch = 1;
				ch = (*ThisField).Pointarr[i].x_po_c[c] * centr[c].P_c;
				float ch1 = 0;
				for (int s = 0; s < k; s++)
				{
					ch1 += (*ThisField).Pointarr[i].x_po_c[s] * centr[s].P_c;
				}
				ch = ch / ch1;
				(*ThisField).Pointarr[i].c_po_x[c] = ch;
				ch = 1;
			}
		}
		std::ofstream outpu;
		outpu.open("output_for_em" + std::to_string(q + 2) + ".txt");
		for (int i = 0; i < (*ThisField).numb_Point; i++)
		{
			outpu << (*ThisField).Pointarr[i].x << " " << (*ThisField).Pointarr[i].y << " " << (*ThisField).Pointarr[i].c_po_x[0] << " ";
			float s = (*ThisField).Pointarr[i].c_po_x[0];
			int t = 0;
			for (int c = 0; c < k; c++)
			{
				if ((*ThisField).Pointarr[i].c_po_x[c] > s)
				{
					s = (*ThisField).Pointarr[i].c_po_x[c];
					t = c;
				}
			}
			outpu << t << std::endl;
		}
		outpu.close();
		outpu.open("ouput_for_em_ell" + std::to_string(q + 2) + ".txt");
		for (int c = 0; c < k; c++)
		{
			std::vector<std::vector<float>> vecs = sobstv_vector(eps[c]);
			int ind_max = (vecs[2][0] > vecs[2][1]) ? 0 : 1;
			if (vecs[ind_max][1] < 0)
			{
				vecs[ind_max][0] *= -1;
				vecs[ind_max][1] *= -1;
			}
			float angle = acos(vecs[ind_max][0] / sqrt(vecs[ind_max][0] * vecs[ind_max][0] + vecs[ind_max][1] * vecs[ind_max][1]));
			angle *= 180;
			angle /= asin(1) * 2;
			outpu << mu[c][0] << " " << mu[c][1] << " " << 20 * vecs[2][1 - ind_max] << " " << 20 * vecs[2][ind_max] << " " << angle << std::endl;
		}
		outpu.close();
		eps.clear();
		mu.clear();
	}
}