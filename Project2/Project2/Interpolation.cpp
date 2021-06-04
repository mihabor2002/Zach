#include "Interpolation.h"
#include "k_sredn.h"
float function(float x)
{
    return exp(- x);
}
void print_itog(Point Z, float y)
{
    ofstream f("interp_gnu.txt"), g("interp.txt");
    g << Z.x << " " << Z.y << " " << y;
    g.close();
    g << "splot (x*x-y*y)/10+5, 'Point.txt'\n";
    g << "set arrow from " << Z.x << "," << Z.y << "," << 0 << " to " << Z.x << "," << Z.y << "," << y << "\n";
}
float Interpolation::funk(Field* TheField, Point Z)
{
    k_sredn alg;
    vector <Point> points, cluster, neighbouring_points, points_1;
    int i, j, k = 0;
    float distant = 100, h = 0.1, sum_w = 0, sum_w_y = 0, r = 0, sum_y = 0, sum_eps = 0, mean_y = 0;
    Delaunay_triangulation T;
    vector <double> w, eps;

    for (i = 0; i < (*TheField).numb_Point; i++)
    {
        points.push_back((*TheField).Pointarr[i]);
    }
    alg.funk(TheField, 1);
    Find_Clasters F = (*TheField).Find_Clastersarr[(*TheField).Find_Clastersarr.size() - 1];
    for (i = 0; i < F.Clasterarr.size(); i++)
    {
        if (dist(Z, F.Clasterarr[i].arr[F.Clasterarr[i].arr.size() - 1]) < distant)
        {
            distant = dist(Z, F.Clasterarr[i].arr[F.Clasterarr[i].arr.size() - 1]);
            k = i;
        }
    }
    if (distant > 10)
    {
        cout << "Error. Points is too far away from the clusters.\n\n";
    }
    else
    {
        for (i = 0; i < F.Clasterarr[k].arr.size() - 1; i++)
        {
            cluster.push_back(F.Clasterarr[k].arr[i]);
        }
        cluster.push_back(Z);
        T = (*TheField).generate_delaunay_trinagulation(cluster);
        T.create_triangle_indicators();
        neighbouring_points = T.find_neighbouring_points(Z);
        cout << "nbp " << neighbouring_points.size();
        cout << function(6);
        for (i = 0; i < neighbouring_points.size(); i++)
        {
            cout <<" dist = " <<dist(Z, neighbouring_points[i]) / h;
            w.push_back(function(dist(Z, neighbouring_points[i])));
            cout << "w[j]" << w[i];
        }
        for (i = 0; i < neighbouring_points.size(); i++)
        {
            sum_w_y = sum_w_y + points[i].get_function_value() * w[i];
            sum_w = sum_w + w[i];
        }
        cout << "res = " << sum_w_y / sum_w;
        cluster.pop_back();
        for (i = 0; i < cluster.size(); i++)
        {
            mean_y = mean_y + points[i].get_function_value();
        }
        mean_y = mean_y * (1 / cluster.size());
        for (i = 0; i < cluster.size(); i++)
        {
            neighbouring_points.clear();
            w.clear();
            points_1.clear();
            for (j = 0; j < i; j++) points_1.push_back(cluster[i]);
            for (j = i + 1; j < points.size(); j++) points_1.push_back(cluster[i]);
            T = (*TheField).generate_delaunay_trinagulation(points_1);
            T.create_triangle_indicators();
            neighbouring_points = T.find_neighbouring_points(points[i]);
            for (j = 0; j < neighbouring_points.size(); j++) w.push_back(function(dist(cluster[i], neighbouring_points[i]) / h));
            for (j = 0; j < neighbouring_points.size(); j++)
            {
                sum_w_y = sum_w_y + (points_1[j].get_function_value()) * w[i];
                sum_w = sum_w + w[i];
            }
            eps.push_back(cluster[i].get_function_value() - sum_w_y / sum_w);
        }
        for (i = 0; i < cluster.size(); i++)
            sum_y = sum_y + ((cluster[i].get_function_value() - mean_y) * (cluster[i].get_function_value() - mean_y));
        for (i = 0; i < cluster.size(); i++) sum_eps = sum_eps + eps[i] * eps[i];
        r = 1 - sum_eps / sum_y;
        cout << "Real: " << Z.get_function_value() << ", sun_w"<<sum_w<<"  Forecast: " << sum_w_y / sum_w << ". r^2=" << r << "." << endl;
        print_itog(Z, sum_w_y / sum_w);
    }
    return sum_w_y / sum_w;
}
