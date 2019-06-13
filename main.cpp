// #define FRUSTRATION_OFF
#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <array>
#include <cmath>
#include <sstream>
#include <fstream>
#include <filesystem>

#include <ncurses.h>    // TUI
#include <unistd.h>     // usleep(microseconds)

#include "hexagonal.h"


// easy-to-use timing functionality
std::vector<std::chrono::high_resolution_clock::time_point> time_points;

void push_timer() { 
    time_points.push_back(std::chrono::high_resolution_clock::now());
}

std::string pop_timer() { 
    auto t2 = std::chrono::high_resolution_clock::now();
    auto t1 = time_points.back();
    time_points.pop_back();
    std::ostringstream ss;
    ss << 
        std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()
        << "μs";
    return ss.str();
}





int main()
{
    constexpr int nx = 20, ny = 20, nt = 50;
    double beta = 1., field = 1.;
    hexagonal<nx, ny, nt> lattice(beta, field);


    // information reported on screen
    std::string last_update_duration; // to report the time cost of updating
    double Gamma = lattice.get_Gamma(),
           delta = lattice.get_delta(),
           magnetization = 0., energy = 0.;
    int n_update = 0; // counter for number of updates performed;
    int t = 0; // time-slice to be drawn on screen
    int last_cluster_size = 0;
    int n_flips_succeeded = 0; // #spins local-updating succeeded in flipping
    bool save_mag = false; /* saving magnetization to a list to calculate
                              magnetic susceptibility later */


    // parameters window region
    const int param_win_L = 0, param_win_T = 0, param_win_W = 80, 
          param_win_H = 5;
    WINDOW *param_win;




    std::string data_file_path = "./data_file";

    // if an older file already exists, move the old file
    if (std::filesystem::exists(data_file_path)) {
        const int max_it = 1000;
        bool moved = false;
        for (int it = 0; it < max_it; ++it) {
            std::string moveto_file_path = data_file_path + "_old_" 
                + std::to_string(it);
            if (!std::filesystem::exists(moveto_file_path)) {
                std::filesystem::rename(data_file_path, moveto_file_path);
                moved = true;
                break;
            }
        }
        if (!moved) { 
            std::cerr << "couldn't move existing file" << std::endl;
            exit(1);
        }
    }

    std::ofstream data_file(data_file_path);
    if (!data_file) {
        std::cerr << "cannot open file for writing\nexiting ..." << std::endl;
        exit(1);
    }




    // initialize and configure ncurses
    initscr();
    curs_set(0);
    cbreak();
    noecho();
    start_color();
    init_pair(1, COLOR_BLUE, COLOR_BLACK);
    init_pair(2, COLOR_BLUE, COLOR_WHITE);
    param_win = newwin(param_win_H, param_win_W, param_win_L,
            param_win_T);
    nodelay(param_win, TRUE); /* don't let wgetch stop the execution 
                                 waiting for a keystroke */




    auto draw_param_win = [&]() {
        wmove(param_win, 1, 1);
        wprintw(param_win, "TmSlic:%d beta:%.3f h:%.3f size:%dx%dx%d "
                "Gamma:%.5f delta:%.3f ", t,
                beta, field, nx, ny, nt, last_cluster_size, n_flips_succeeded,
                Gamma, delta);
        wmove(param_win, 2, 1);
        wprintw(param_win, "mag'n:%.4f energy:%.4f ", magnetization, energy);

        if (save_mag) {
            wprintw(param_win, "suce:%.4f ", lattice.mag_susc());
        }
        wmove(param_win, 3, 1);
        wprintw(param_win, "update#:%d"
                " LstClusSiz:%d LstUpdDur:%s", n_update,
                last_cluster_size, last_update_duration.c_str());
        box(param_win, 0, 0);
        wrefresh(param_win);
    };



    auto draw_state = [&lattice, &t]() {
        for (int sublattice = 0; sublattice < 2; ++sublattice) {
            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {
                    attron(COLOR_PAIR(2-x%2));
                    move(2 * y + sublattice + param_win_T + param_win_H, 
                            6*x + 3*(y%2));
                    addch(' ');
                    addch(lattice.get_site_spin(x, y, t, sublattice) == +1 
                            ? '+' : '-');
                    addch(' ');
                }
            }
        }
        refresh();
    };



    int key;
    do {
        key = wgetch(param_win);
        usleep(1000);

        /* keys:
         * u/U: update the state
         * t/T: go to previous/next time-slice
         * h/H: decrease/increase the magnetic field by a factor of 1.1
         * b/B: decrease/increase the inverse-temperature by a factor of 1.1
         */
        switch (key) {
            case 'u': // update
            case 'U':
                push_timer();
                {
                    int total = 0;
                    while (total < 20 * nx * ny * nt) {
                        // it seems that the cluster evolution isn't ergodic?
                        // so, is this loop a valid solution for this?
                        for (int j = 0; j < 100; j++) {
                            n_flips_succeeded = lattice.local_update();
                        }

                        last_cluster_size = lattice.cluster_update();
                        total += last_cluster_size;
                        magnetization = lattice.magnetization();
                        energy = lattice.energy();
                        n_update++;
                        data_file << n_update << ' ' << magnetization << std::endl;
                    }
                }
                last_update_duration = pop_timer();

                if (save_mag) {
                    lattice.save_mag();
                }
                break;
            case 't': 
                t += nt - 1; t%=nt;
                break;
            case 'T': 
                t++; t%=nt;
                break;
            case 'h': 
                field /= 1.1;
                lattice.set_params(beta, field);
                Gamma = lattice.get_Gamma(); delta = lattice.get_delta();
                break;
            case 'H':
                field *= 1.1;
                lattice.set_params(beta, field);
                Gamma = lattice.get_Gamma(); delta = lattice.get_delta();
                break;
            case 'b':
                beta /= 1.1;
                lattice.set_params(beta, field);
                Gamma = lattice.get_Gamma(); delta = lattice.get_delta();
                break;
            case 'B':
                beta *= 1.1;
                lattice.set_params(beta, field);
                Gamma = lattice.get_Gamma(); delta = lattice.get_delta();
                break;
            case 'm':
            case 'M':
                if (save_mag) {
                    save_mag = false;
                    lattice.clear_mag_history();
                } else {
                    save_mag = true;
                }
                break;
            case 'r':
            case 'R':
                break;

        }
        draw_param_win();
        draw_state();
    } while (key != 'q');



    // finalize TUI
    delwin(param_win);
    endwin();

    return 0;
}

