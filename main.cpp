#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <array>
#include <cmath>
#include <sstream>

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
    constexpr int nx = 10, ny = 10, nt = 10;

    // initialize and configure ncurses
    initscr();
    curs_set(0);
    cbreak();
    noecho();
    start_color();
    init_pair(1, COLOR_BLUE, COLOR_BLACK);
    init_pair(2, COLOR_BLUE, COLOR_WHITE);


    double beta = 10., field = 1.;
    double Gamma, delta, magnetization = 0.;
    int n_update = 0; // counter for number of updates performed;
    std::string last_update_duration;
    int t = 0; // time-slice to be drawn on screen
    hexagonal<nx, ny, nt> lattice(beta, field);
    Gamma = lattice.get_Gamma(); delta = lattice.get_delta();


    const int param_win_L = 0, param_win_T = 0, param_win_W = 80, 
          param_win_H = 5;
    WINDOW *param_win = newwin(param_win_H, param_win_W, param_win_L,
            param_win_T);
    nodelay(param_win, TRUE); /* don't let wgetch stop the execution 
                                 waiting for a keystroke */


    int last_cluster_size = 0;
    int n_flips_succeeded = 0;

    auto draw_param_win = [&]() {
        wmove(param_win, 1, 1);
        wprintw(param_win, "TmSlic:%d beta:%.3f h:%.3f size:%dx%dx%d "
                "Gamma:%.5f delta:%.3f ", t,
                beta, field, nx, ny, nt, last_cluster_size, n_flips_succeeded,
                Gamma, delta);
        wmove(param_win, 2, 1);
        wprintw(param_win, "mag'n:%.4f", magnetization);
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

        switch (key) {
            case 'u':
            case 'U':
                push_timer();

                // it seems that the cluster evolution isn't ergodic?
                // so, is this loop a valid solution for this?
                for (int i = 0; i < 1000; i++) {
                    n_flips_succeeded = lattice.local_update();
                }

                last_cluster_size = lattice.cluster_update();
                magnetization = lattice.magnetization();
                n_update++;
                last_update_duration = pop_timer();
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
        }
        draw_param_win();
        draw_state();
    } while (key != 'q');

    // finalize TUI
    delwin(param_win);
    endwin();

    return 0;
}

