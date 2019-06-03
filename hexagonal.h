#ifndef _HEXAGONAL_HEADER__
#define _HEXAGONAL_HEADER__

#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <random>

// coordinates of a site. sublattice={0, 1}
struct coords { int x, y, t, sublattice; };

template<int nx, int ny, int nt, int rand_seed = 1999>
class hexagonal // actually, honeycomb
{
private:
    std::mt19937 rand_gen{rand_seed};

private: 
    double beta, field;
    double delta, Gamma;

private:
    static constexpr auto n_cell3 = nx * ny * nt; // #unit cells in 2+1D
    static constexpr auto n_cell2 = nx * ny;// #unit cells in a 2D time-slice
    static constexpr auto n_site3 = 2 * n_cell3; // we have 2 atoms / unit cell
    static constexpr auto n_nbrs = 5; /* we have 5 nearest neighbors, i.e.,
                                        3 (spatial) + 2 (temporal axis) */
private:
    std::array<int, n_site3> state; // configuration of the classical system
    std::array<std::array<int, n_nbrs>, n_site3> nbr_indices; /* we will 
                                    calculate the indices for the neighbors 
                                    of each site once and for all */
private:
    // Return true with probability p
    bool prob_true(double p)
    {
        static const auto max = rand_gen.max();
        return (rand_gen() < p * max);
    }

public:
    int get_site_spin(int x, int y, int t, int sublattice) 
    {
        int index = coords_to_index(x, y, t, sublattice);
        return state[index];
    }

public:
    hexagonal(double beta_p, double field_p)
    {
        set_params(beta_p, field_p);
        generate_nbr_indices();
        generate_random_state();
    }

    void set_params(double beta_p, double field_p)
    {
        beta = beta_p;
        field = field_p;

        delta = beta / nt;

        double eta = -.5 * std::log(std::tanh(delta * field));
        Gamma = eta / delta;
    }

    void generate_nbr_indices()
    {
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int t = 0; t < nt; ++t) {
                    int index;
                    int sublattice;


                    sublattice = 0;
                    index = coords_to_index(x, y, t, sublattice);

                    nbr_indices[index][0] = 
                        coords_to_index(x, y, t+1, sublattice);
                    nbr_indices[index][1] = 
                        coords_to_index(x, y, t-1, sublattice);
                    nbr_indices[index][2] = 
                        coords_to_index(x, y-1, t, 1-sublattice);
                    nbr_indices[index][3] = 
                        coords_to_index(x+1, y-1, t, 1-sublattice);
                    nbr_indices[index][4] = 
                        coords_to_index(x, y, t, 1-sublattice);


                    sublattice = 1;
                    index = coords_to_index(x, y, t, sublattice);

                    nbr_indices[index][0] = 
                        coords_to_index(x, y, t+1, sublattice);
                    nbr_indices[index][1] = 
                        coords_to_index(x, y, t-1, sublattice);
                    nbr_indices[index][2] = 
                        coords_to_index(x, y+1, t, 1-sublattice);
                    nbr_indices[index][3] = 
                        coords_to_index(x+1, y+1, t, 1-sublattice);
                    nbr_indices[index][4] = 
                        coords_to_index(x, y, t, 1-sublattice);

                }
            }
        }
        // for (int i = 0; i < 5; ++i) {
            // print_coords(index_to_coords(nbr_indices[coords_to_index(0,0,0,1)][i]));
        // }
    }

    void generate_random_state()
    {
        for (auto it = state.begin(); it != state.end(); ++it) {
            *it = 1 - 2 * (rand_gen() % 2);
        }
    }


    int coords_to_index(int x, int y, int t, int sublattice) 
    {

        // take care of periodic boundary conditions without
        // (expensive) if statements
        x += nx; y += ny; t += nt;
        x %= nx; y %= ny; t %= nt;

        return sublattice * n_cell3 + t * n_cell2 + y * nx + x;
    }

    coords index_to_coords(int index)
    {
       coords ret;

       ret.sublattice = index / n_cell3;
       index %= n_cell3;

       ret.t = index / n_cell2;
       index %= n_cell2;

       ret.y = index / nx;
       index %= nx;

       ret.x = index;

       return ret;
    }


    int local_update()
    {
        int rand_index = rand_gen() % n_site3;
        int spin = state[rand_index];
        double flip_cost = 0.;

        int nbr_spin;

        // temporal axis bonds
        for (int i = 0; i < 2; ++i) {
            nbr_spin = state[nbr_indices[rand_index][i]];
            flip_cost += (+2.) * Gamma * spin * nbr_spin;
        }

        // spatial (always-ferromagnetic) bonds
        for (int i = 2; i < n_nbrs - 1; ++i) {
            nbr_spin = state[nbr_indices[rand_index][i]];
            flip_cost += (+2.) * spin * nbr_spin;
        }


        // last bond is anti-ferromagnetic exactly for even-x sites:
        nbr_spin = state[nbr_indices[rand_index][n_nbrs-1]];
        flip_cost += (rand_index % 2 == 0 ? (-1) : (+1)) *
            (+2.) * spin * nbr_spin;


        static int flips_succeeded = 0;

        double flip_prob = std::exp(-delta * flip_cost);
        if (prob_true(flip_prob)) {
            state[rand_index] *= -1;
                flips_succeeded++;
        }

        return flips_succeeded;
    }

    int cluster_update()
    {
        std::vector<int> cluster, old_set, new_set;
        int seed_index = rand_gen() % n_site3;// index of first site on cluster
        cluster.push_back(seed_index);
        old_set.push_back(seed_index);

        std::vector<bool> on_cluster; /* using which, it takes O(1) to check
                                   if a given site belongs to the cluster;
                                   no need to do a search every time. */
        on_cluster.resize(n_site3);
        for (auto it = on_cluster.begin(); it != on_cluster.end(); ++it) {
            *it = false;
        }
        on_cluster[seed_index] = true;


        const double p_xy = 1. - std::exp(-2. * delta);
        const double p_t = 1. - std::exp(-2. * delta * Gamma);

        while (!old_set.empty()) {
            new_set.clear();
            for (auto old_it = old_set.begin(); old_it != old_set.end(); 
                    ++old_it) {
                int my_spin_state = state[*old_it];

                // the upcoming for-loop will only involve ferro. bonds
                int n_ferro_bonds = n_nbrs;
                if (*old_it % 2 == 0) { // for even x, the last bond is anti-f
                    n_ferro_bonds--;
                    int last_nbr_index = nbr_indices[*old_it][n_nbrs-1];
                    if (!on_cluster[last_nbr_index]) {
                        int nbr_spin_state = state[last_nbr_index];
                        if (nbr_spin_state != my_spin_state) {
                            if (prob_true(p_xy)) {
                                new_set.push_back(last_nbr_index);
                                cluster.push_back(last_nbr_index);
                                on_cluster[last_nbr_index] = true;
                            }
                        }
                    }
                }

                for (auto nbr_it = 0; nbr_it < n_ferro_bonds; ++nbr_it) {
                    int nbr_index = nbr_indices[*old_it][nbr_it];
                    if (!on_cluster[nbr_index]) {
                        int nbr_spin_state = state[nbr_index];
                        if (nbr_spin_state == my_spin_state) {
                            double p = nbr_it < 2 ? p_t : p_xy; /* is it a 
                                        temporal or a spatial bond? */
                            if (prob_true(p)) {
                                new_set.push_back(nbr_index);
                                cluster.push_back(nbr_index);
                                on_cluster[nbr_index] = true;
                            }
                        }
                    }
                }
            }
            old_set = new_set;
        }

        for (auto it = cluster.begin(); it != cluster.end(); ++it) {
            state[*it] *= -1;
        }

        return cluster.size();
    }

    double magnetization()
    {
        // operator m is diagonal, so we evaluate its average using
        // the first time-slice only:
        double M = 0;
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int sublattice = 0; sublattice < 2; ++sublattice) {
                    const int t = 0;
                    M += state[coords_to_index(x, y, t, sublattice)];
                }
            }
        }
        return M / (nx * ny * 2);
    }


    std::string coords_to_string (coords c)
    {
        std::ostringstream ss;
        ss << "(" << c.x << ", " << c.y
            << ", " << c.t << ", " << c.sublattice << ")";
        return ss.str();
    }

    double get_Gamma() { return Gamma; }
    double get_delta() { return delta; }
};

#endif
