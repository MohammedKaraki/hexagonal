#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>


#include "hexagonal.h"

using std::cout, std::endl;
using std::vector;

double mean(const vector<double> &v);
double standard_error(const vector<double> &v);


// easy-to-use timing functionality
std::vector<std::chrono::high_resolution_clock::time_point> time_points;

void push_timer()
{
    time_points.push_back(std::chrono::high_resolution_clock::now());
}

std::string last_timer()
{
    auto t2 = std::chrono::high_resolution_clock::now();
    auto t1 = time_points.back();
    std::ostringstream ss;
    ss << 
        std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
        << "ms";
    return ss.str();
}

std::string pop_timer() 
{
    std::string ret = last_timer();
    time_points.pop_back();
    return ret;
}

// file notation
// nx ny field_max beta_init

constexpr int nx = 10, ny = 10, nt = 50;
double beta_init = 2.;
constexpr double field_min = .1;
constexpr double field_max = 5.;
constexpr double field_int = .1;
// constexpr int n_field_vals = (field_max - field_min) / field_int;



int main()
{
    hexagonal<nx, ny, nt> lattice (beta_init, field_min);
    cluster_update_metadata metadata;


    auto thermalize = [&](int n1, int n2, int n3) { 
        int total_spatial = 0, total_temporal = 0, total_both = 0;
        int N = 0;
        while (total_spatial < n1 * (nx * ny) ||
                total_temporal < n2 * nt ||
                total_both < n3 * (nx * ny * nt)) {

            lattice.cluster_update(&metadata);

            N++;

            total_temporal += metadata.n_temporal_links;
            total_spatial += metadata.n_spatial_links;
            total_both += metadata.n_all_links;
        }

        return N;
    };

    auto measure = [&](int n1, int n2, int n3) { 
        int total_spatial = 0, total_temporal = 0, total_both = 0;
        double energy_sum = 0.;
        int N = 0;
        while (total_spatial < n1 * (nx * ny) ||
                total_temporal < n2 * nt ||
                total_both < n3 * (nx * ny * nt)) {

            lattice.cluster_update(&metadata);

            N++;
            energy_sum += lattice.energy();

            total_temporal += metadata.n_temporal_links;
            total_spatial += metadata.n_spatial_links;
            total_both += metadata.n_all_links;
        }

        return energy_sum / N;
    };

    cout << nx << "x" << ny << "x" << nt << endl;
    push_timer();
    for (double beta = beta_init; beta < 10. * beta_init; beta *= 1.5) {

        push_timer();
        thermalize(200, 100, 200);
        cout << "\nthermalization time: " << pop_timer() << endl;

        for (double field = field_min; field <= field_max; field += field_int)
        {
            push_timer();
            lattice.set_params(beta, field);

            thermalize(20, 20, 10);

            vector<double> data;
            vector<double> mag;
            constexpr int sample_size = 25;
            for (int i = 0; i < sample_size; ++i) {
                data.push_back(measure(2, 1, 2));
                mag.push_back(std::abs(lattice.magnetization()));
            }

            cout << beta << " \t" << field << " \t" << mean(data) << " \t" << 
                standard_error(data) << " \t" << mean(mag) << " \t" << 
                standard_error(mag) << " \t" << pop_timer() << " \t" <<
                last_timer() << endl;
        }
    }
    cout << "total time: " << pop_timer() << endl;

    return 0;
}


double mean(const vector<double> &v)
{
    return std::accumulate(v.begin(), v.end(), 0.) / v.size();
}

double standard_error(const vector<double> &v)
{
    double sample_mean = mean(v);
    double sq_sum = 0.;
    for (auto it = v.begin(); it != v.end(); ++it) {
        double val = *it;
        sq_sum += val*val;
    }
    double N = v.size();
    double sq_mean = sq_sum / N;
    double tmp1 = sq_mean - sample_mean * sample_mean;
    double tmp2 = tmp1 / (N-1);
    return std::sqrt(tmp2);
}

