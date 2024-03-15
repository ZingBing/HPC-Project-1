/**
 * Declares all of the math helper functions (which are defined in math_helper.c).
 */

/**
* Functions used to calculate the gravitational force between two bodies
* Used in all four implementations of the n-body problem
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"

// Gravitational Constant in N m^2 / kg^2 or m^3 / kg / s^2
#define G 6.6743015e-11

// Softening factor to reduce divide-by-near-zero effects
#define SOFTENING 1e-9

/**
* Calculate the euclidean distance between two bodies in 3D space
*/
inline static double euclidean_distance_sans_sqrt(double* body_1_position, double* body_2_position) {
    double x = body_1_position[0] - body_2_position[0];
    double y = body_1_position[1] - body_2_position[1];
    double z = body_1_position[2] - body_2_position[2];
    return (x * x + y * y + z * z + SOFTENING);
}

/**
* Calculate the gravitational force between two bodies
*/
inline static double gravitation(double mass_body_1, double mass_body_2, double* body_1_position, double* body_2_position, double distance) {
    return G * ((mass_body_1 * mass_body_2) / distance);
}

/**
* Calculate the acceleration of a body due to a force
*/
inline static double get_acceleration(double force, double mass) {
    return force / mass;
}

inline static void save_position(Matrix* output, double* position, size_t output_row, size_t n) {
    for (size_t i = 0; i < n; i++) {
        output->data[output_row*3*n + 3*i] = position[3*i];
        output->data[output_row*3*n + 3*i+1] = position[3*i+1];
        output->data[output_row*3*n + 3*i+2] = position[3*i+2];
        }
}
