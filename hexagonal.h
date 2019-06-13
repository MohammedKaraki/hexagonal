#ifndef _HEXAGONAL_HEADER__
#define _HEXAGONAL_HEADER__

#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <random>

// coordinates of a site. sub_lattice={0, 1}
struct coords
{ 
    int x, y, t, sub_lattice; 
};

// feedback from the algorithm
struct cluster_update_metadata
{
    int n_all_links, n_spatial_links, n_temporal_links;
};

template<int nx, int ny, int nt, int rand_seed = 1999>
class hexagonal // actually, honeycomb
{
private:
    std::mt19937 rand_gen{rand_seed};

private: 
    double beta, field;
    double delta, eta, Gamma;

private:
    static constexpr auto n_cell3 = nx * ny * nt; // #unit cells in 2+1D
    static constexpr auto n_cell2 = nx * ny;// #unit cells in a 2D time-slice
    static constexpr auto n_site3 = 2 * n_cell3;// we have 2 atoms / unit cell
    static constexpr auto n_site2 = 2 * n_cell2;
    static constexpr auto n_nbrs = 5; /* we have 5 nearest neighbors, i.e.,
                                        3 (spatial) + 2 (temporal axis) */
private:
    std::array<int, n_site3> state; // configuration of the classical system
    std::array<std::array<int, n_nbrs>, n_site3> nbr_indices; /* we will 
                                    calculate the indices for the neighbors 
                                    of each site once and for all */

    // a list of magnetization history for calculating mag_susc.
    // Related functions: save_mag, clear_mag_history, mag_susc.
    std::vector<double> mag_history;

private:
    // return true with probability p
    bool prob_true(double p)
    {
        static const auto max = rand_gen.max();
        return (rand_gen() < p * max);
    }

public:
    int get_site_spin(int x, int y, int t, int sub_lattice) 
    {
        int index = coords_to_index(x, y, t, sub_lattice);
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

        eta = -.5 * std::log(std::tanh(delta * field));
        Gamma = eta / delta;
    }

    void generate_nbr_indices()
    {
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int t = 0; t < nt; ++t) {
                    int index;
                    int sub_lattice;

                    /* here is a diagram showing the links with nearest 
                       neighbors of a site on the 0 sub_lattice. The short
                       link means it's sometimes-ferromagnetic.
                            | /
                            |/
                            0
                            |     */
                    sub_lattice = 0;
                    index = coords_to_index(x, y, t, sub_lattice);

                    // temporal neighbors
                    nbr_indices[index][0] = 
                        coords_to_index(x, y, t+1, sub_lattice);
                    nbr_indices[index][1] = 
                        coords_to_index(x, y, t-1, sub_lattice);

                    // spatial always-ferromagnetic exchange
                    nbr_indices[index][2] = 
                        coords_to_index(x, y-1, t, 1-sub_lattice);
                    nbr_indices[index][3] = 
                        coords_to_index(x+1, y-1, t, 1-sub_lattice);

                    // spatial sometimes-ferromagnetic exchange
                    nbr_indices[index][4] = 
                        coords_to_index(x, y, t, 1-sub_lattice);



                    /* here is a diagram showing the links with nearest 
                       neighbors of a site on the 1 sub_lattice. The short
                       link means it's sometimes-ferromagnetic.
                            |
                            1
                           /|
                          / |     */
                    sub_lattice = 1;
                    index = coords_to_index(x, y, t, sub_lattice);

                    // temporal neighbors
                    nbr_indices[index][0] = 
                        coords_to_index(x, y, t+1, sub_lattice);
                    nbr_indices[index][1] = 
                        coords_to_index(x, y, t-1, sub_lattice);

                    // spatial always-ferromagnetic exchange
                    nbr_indices[index][2] = 
                        coords_to_index(x, y+1, t, 1-sub_lattice);
                    nbr_indices[index][3] = 
                        coords_to_index(x-1, y+1, t, 1-sub_lattice);

                    // spatial sometimes-ferromagnetic exchange
                    nbr_indices[index][4] = 
                        coords_to_index(x, y, t, 1-sub_lattice);

                }
            }
        }
    }

    void generate_random_state()
    {
        for (auto it = state.begin(); it != state.end(); ++it) {
            *it = 1 - 2 * (rand_gen() % 2);
        }
    }


    int coords_to_index(int x, int y, int t, int sub_lattice) 
    {

        // take care of periodic boundary conditions without
        // (expensive) if statements
        x += nx; y += ny; t += nt;
        x %= nx; y %= ny; t %= nt;

        return sub_lattice * n_cell3 + t * n_cell2 + y * nx + x;
    }

    coords index_to_coords(int index)
    {
       coords ret;

       ret.sub_lattice = index / n_cell3;
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

        // temporal axis interactions
        for (int i = 0; i < 2; ++i) {
            nbr_spin = state[nbr_indices[rand_index][i]];
            flip_cost += (+2.) * Gamma * spin * nbr_spin;
        }

        // spatial (always-ferromagnetic) interactions
        for (int i = 2; i < n_nbrs - 1; ++i) {
            nbr_spin = state[nbr_indices[rand_index][i]];
            flip_cost += (+2.) * spin * nbr_spin;
        }


        // last bond is antiferromagnetic exactly for even-x sites:
        nbr_spin = state[nbr_indices[rand_index][n_nbrs-1]];
#ifdef FRUSTRATION_OFF
        flip_cost += (+2.) * spin * nbr_spin;
#else
        flip_cost += (rand_index % 2 == 0 ? (-1) : (+1)) *
            (+2.) * spin * nbr_spin;
#endif


        static int flips_succeeded = 0;

        double flip_prob = std::exp(-delta * flip_cost);
        if (prob_true(flip_prob)) {
            state[rand_index] *= -1;
            flips_succeeded++;
        }

        return flips_succeeded;
    }



    int cluster_update(cluster_update_metadata *metadata)
    {
        // feedback variables (metadata)
        int n_temporal_links = 0, n_spatial_links = 0;


        std::vector<int> cluster, old_set, new_set;
        int seed_index = rand_gen() % n_site3; /* index of first site on 
                                                  cluster */
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

                // the upcoming for-loop will involve only ferromagnetic
                // exchange. Here we take care of the exchange with the last
                // neighbor, which can be either, depending on x parity.
                int n_ferro_bonds = n_nbrs;
#ifdef FRUSTRATION_OFF
                if (false) {
#else
                if (*old_it % 2 == 0) { /* for even x, the exchange is
                                           antiferromagnetic */
#endif
                    n_ferro_bonds--;
                    int last_nbr_index = nbr_indices[*old_it][n_nbrs-1];
                    if (!on_cluster[last_nbr_index]) {
                        int nbr_spin_state = state[last_nbr_index];
                        if (nbr_spin_state != my_spin_state) {
                            if (prob_true(p_xy)) {
                                new_set.push_back(last_nbr_index);
                                cluster.push_back(last_nbr_index);
                                on_cluster[last_nbr_index] = true;

                                n_spatial_links++;
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
                                        temporal or a spatial exchange? */
                            if (prob_true(p)) {
                                new_set.push_back(nbr_index);
                                cluster.push_back(nbr_index);
                                on_cluster[nbr_index] = true;

                                if (nbr_it < 2) {
                                    n_temporal_links++;
                                } else {
                                    n_spatial_links++;
                                }
                            }
                        }
                    }
                }
            }
            old_set = new_set;
        }

        /* 
         * is it correct to treat the cluster differently, applying 
         * conditional acceptance of the move, if its size is 1?
         *
         * Currently, this implementation always accepts the cluster, and I
         * don't understand why there should be an exception.
         *
         */

        for (auto it = cluster.begin(); it != cluster.end(); ++it) {
            state[*it] *= -1;
        }

        if (metadata != nullptr) {
            metadata->n_spatial_links = n_spatial_links;
            metadata->n_temporal_links = n_temporal_links;
            metadata->n_all_links = n_spatial_links + n_temporal_links;
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
                for (int sub_lattice = 0; sub_lattice < 2; ++sub_lattice) {
                    const int t = 0;
                    M += state[coords_to_index(x, y, t, sub_lattice)];
                }
            }
        }
        return M / (n_site2);
    }


    std::string coords_to_string (coords c)
    {
        std::ostringstream ss;
        ss << "(" << c.x << ", " << c.y << ", " << c.t << ", " 
            << c.sub_lattice << ")";

        return ss.str();
    }

    double get_Gamma() { return Gamma; }
    double get_delta() { return delta; }


    void save_mag()
    {
        mag_history.push_back(magnetization());
    }

    void clear_mag_history()
    {
        mag_history.clear();
    }


    // return magnetic susceptibility per spin
    // \Chi = \beta N (<m^2> - <m>^2).
    // Note that we only need evaluate the average over a single time-slice
    // since the operator m is diagonal.
    double mag_susc()
    {
        if (mag_history.size() <= 0) {
            return 0.;
        }

        double m_sum = 0., m2_sum = 0.;

        for (auto it = mag_history.begin(); it != mag_history.end(); ++it) {
            double m = *it;
            m_sum += m;
            m2_sum += m*m;
        }

        const int N = mag_history.size();
        double m_avg = m_sum / N;
        double m2_avg = m2_sum / N;

        return beta * n_site2 * (m2_avg - m_avg*m_avg);
    }

    /* e = (U + V)/n_site2;
     * E need to be evaluated using one slice only, be we'll average it over
     * nt slices.
     * For U, we'll loop over the sites on sub_lattice 0 only, each of which
     * has 3 spatial exchanges. This exhausts all the spatial exchanges.
     * For V, we need to consider the two sub_lattices to exhaust all the
     * temporal links.
     *
       */
    double energy()
    {
        // E = U + V
        double U, V;

        U = 0.;
        V = 0.;
        for (int t = 0; t < nt; ++t) { // average over nt slices

            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {

                    int my_index;
                    int my_spin;
                    for (int sub_lattice = 1; sub_lattice >= 0; --sub_lattice)
                    {

                        my_index = 
                            coords_to_index(x, y, t, sub_lattice);
                        my_spin = state[my_index];

                        V += std::exp(-2. * eta * my_spin * 
                                state[nbr_indices[my_index][0]]);
                    }

                    // at this point, my_index and my_spin correspond
                    // to sub_lattce = 0

                    // the always-ferromagnetic exchange contributions
                    U += (-1) * my_spin * (state[nbr_indices[my_index][2]]
                            + state[nbr_indices[my_index][3]]);

                    // the sometimes-ferromagnetic exchange contribution
                    int last_exchange_sign = 1 - 2*(x%2); /*
                          = +1 (anti-ferromagnetic) if x is even
                          = -1 (ferromagnetic) if x is odd */
#ifdef FRUSTRATION_OFF
                    last_exchange_sign = -1;
#endif
                    U += last_exchange_sign * my_spin *
                        state[nbr_indices[my_index][4]];
                }
            }

        }
        V = V * (-field) / nt;
        U /= nt;

        return (U + V) / n_site2;
    }

};

#endif

