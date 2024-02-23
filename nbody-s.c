/**
 * Runs a simulation of the n-body problem in 3D.
 * 
 * To compile the program:
 *   gcc -Wall -O3 -march=native nbody-s.c matrix.c util.c -o nbody-s -lm
 * 
 * To run the program:
 *   ./nbody-s time-step total-time outputs-per-body input.npy output.npy
 * where:
 *   - time-step is the amount of time between steps (Δt, in seconds)
 *   - total-time is the total amount of time to simulate (in seconds)
 *   - outputs-per-body is the number of positions to output per body
 *   - input.npy is the file describing the initial state of the system (below)
 *   - output.npy is the output of the program (see below)
 * 
 * input.npy has a n-by-7 matrix with one row per body and the columns:
 *   - mass (in kg)
 *   - initial x, y, z position (in m)
 *   - initial x, y, z velocity (in m/s)
 * 
 * output.npy is generated and has a (outputs-per-body)-by-(3n) matrix with each
 * row containing the x, y, and z positions of each of the n bodies after a
 * given timestep.
 * 
 * See the PDF for implementation details and other requirements.
 * 
 * AUTHORS: Zachery Bingaman, Seth Coleman
 */

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "matrix.h"
#include "util.h"

// Gravitational Constant in N m^2 / kg^2 or m^3 / kg / s^2
#define G 6.6743015e-11

// Softening factor to reduce divide-by-near-zero effects
#define SOFTENING 1e-9


int main(int argc, const char* argv[]) {
    // parse arguments
    if (argc != 6 && argc != 7) { fprintf(stderr, "usage: %s time-step total-time outputs-per-body input.npy output.npy [num-threads]\n", argv[0]); return 1; }
    double time_step = atof(argv[1]), total_time = atof(argv[2]);
    if (time_step <= 0 || total_time <= 0 || time_step > total_time) { fprintf(stderr, "time-step and total-time must be positive with total-time > time-step\n"); return 1; }
    size_t num_outputs = atoi(argv[3]);
    if (num_outputs <= 0) { fprintf(stderr, "outputs-per-body must be positive\n"); return 1; }
    Matrix* input = matrix_from_npy_path(argv[4]);
    if (input == NULL) { perror("error reading input"); return 1; }
    if (input->cols != 7) { fprintf(stderr, "input.npy must have 7 columns\n"); return 1; }
    size_t n = input->rows;
    if (n == 0) { fprintf(stderr, "input.npy must have at least 1 row\n"); return 1; }
    size_t num_steps = (size_t)(total_time / time_step + 0.5);
    if (num_steps < num_outputs) { num_outputs = 1; }
    size_t output_steps = num_steps/num_outputs;
    num_outputs = (num_steps+output_steps-1)/output_steps;

    // variables available now:
    //   time_step    number of seconds between each time point
    //   total_time   total number of seconds in the simulation
    //   num_steps    number of time steps to simulate (more useful than total_time)
    //   num_outputs  number of times the position will be output for all bodies
    //   output_steps number of steps between each output of the position
    //   input        n-by-7 Matrix of input data
    //   n            number of bodies to simulate

    // start the clock
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // This all came from copilot, run it before we trust it next time
    Matrix* output = matrix_copy(input);
    // compute for all time steps
    for (size_t step = 0; step < num_steps; step++) {
        // compute the force on each body
        for (size_t i = 0; i < n; i++) {
            double force[3] = {0, 0, 0};
            for (size_t j = 0; j < n; j++) {
                if (i == j) { continue; }
                double r[3] = {
                    output->data[i*7+1] - output->data[j*7+1],
                    output->data[i*7+2] - output->data[j*7+2],
                    output->data[i*7+3] - output->data[j*7+3]
                };
                double dist = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
                double mag = G * output->data[i*7] * output->data[j*7] / (dist*dist + SOFTENING*SOFTENING);
                force[0] -= mag * r[0] / dist;
                force[1] -= mag * r[1] / dist;
                force[2] -= mag * r[2] / dist;
            }
            output->data[i*7+4] += force[0] / output->data[i*7];
            output->data[i*7+5] += force[1] / output->data[i*7];
            output->data[i*7+6] += force[2] / output->data[i*7];
        }
        // update the position of each body
        for (size_t i = 0; i < n; i++) {
            output->data[i*7+1] += time_step * output->data[i*7+4];
            output->data[i*7+2] += time_step * output->data[i*7+5];
            output->data[i*7+3] += time_step * output->data[i*7+6];
        }
    }

    // get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f secs\n", time);

    // save results
    //matrix_to_npy_path(argv[5], output);

    // cleanup


    return 0;
}
