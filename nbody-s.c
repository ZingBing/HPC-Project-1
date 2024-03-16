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
#include "helper_functions.h"

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

    //print input matrix to know initial body positions
    // printf("Initial Body Positions:\n");
    // for (size_t i = 0; i < n; i++) {
    //     printf("Body %zu: %f, %f, %f\n", i, input->data[i*input->cols + 1], input->data[i*input->cols + 2], input->data[i*input->cols + 3]);
    // }

    // start the clock
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // here's where we'll be writing code, copilot did a bad
    // allocate output matrix as num_outputs x 3*n 
    Matrix* output = matrix_create_raw(num_outputs, 3*n);
    if (output == NULL) { perror("error allocating output"); return 1; }
    
    // Create simple variable for number of columns in input matrix
    size_t icols = input->cols;
    
    // Create variables for position and velocity of each body
    double* position = malloc(3 * n * sizeof(double));
    double* velocity = malloc(3 * n * sizeof(double));
    // double* force = malloc(3 * n * sizeof(double));
    double* mass = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        mass[i] = input->data[i*icols];
        position[3*i] = input->data[i*icols + 1];
        position[3*i+1] = input->data[i*icols + 2];
        position[3*i+2] = input->data[i*icols + 3];
        velocity[3*i] = input->data[i*icols + 4];
        velocity[3*i+1] = input->data[i*icols + 5];
        velocity[3*i+2] = input->data[i*icols + 6];
    }

    // Save positions to row `0` of output
    save_position(output, position, 0, n);


    // print the initial positions of the bodies to see if they are correct
    // for (size_t i = 0; i < n; i++) {
    //     printf("Body %zu: %f, %f, %f\n", i, position[3*i], position[3*i+1], position[3*i+2]);
    // }

    // Run simulation for each time step 
    for (size_t t = 1; t < num_steps; t++) { 
        // compute time step...
        for (size_t i = 0; i < n; i++) {
            double x_force = 0;
            double y_force = 0;
            double z_force = 0;

            for (size_t j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                double distance = euclidean_distance_sans_sqrt(&position[i*3], &position[j*3]);
                double force = gravitation(mass[i], mass[j], &position[i*3], &position[j*3], distance);

                double accel = force / sqrt(distance);
                
                x_force += accel * (position[j*3+0] - position[i*3+0]);
                y_force += accel * (position[j*3+1] - position[i*3+1]);
                z_force += accel * (position[j*3+2] - position[i*3+2]);
            }

            double x_accel = get_acceleration(x_force, mass[i]);
            double y_accel = get_acceleration(y_force, mass[i]);
            double z_accel = get_acceleration(z_force, mass[i]);

            // Numerically integrate acceleration to get velocity
            velocity[3*i+0] += x_accel * time_step;
            velocity[3*i+1] += y_accel * time_step;
            velocity[3*i+2] += z_accel * time_step;

            // Numerically integrate velocity to get position
            position[3*i+0] += velocity[3*i+0] * time_step;
            position[3*i+1] += velocity[3*i+1] * time_step;
            position[3*i+2] += velocity[3*i+2] * time_step;
        }
    
        // Periodically copy the positions to the output data 
        if (t % output_steps == 0) { 
            // Save positions to row `t/output_steps` of output
            size_t output_row = t / output_steps;
            save_position(output, position, output_row, n);
        }
    } 
    
    // Save the final set of data if necessary 
    if (num_steps % output_steps != 0) { 
        // TODO: save positions to row `num_outputs-1` of output 
        save_position(output, position, num_outputs-1, n);
    }


    // get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f secs\n", time);

    // save results
    matrix_to_npy_path(argv[5], output);

    // cleanup
    matrix_free(output);
    matrix_free(input);
    free(position);
    free(velocity);
    free(mass);

    return 0;
}
