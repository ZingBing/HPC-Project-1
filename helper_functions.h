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
inline static double euclidean_distance(double* body_1_position, double* body_2_position) {
    double x = body_1_position[0] - body_2_position[0];
    double y = body_1_position[1] - body_2_position[1];
    double z = body_1_position[2] - body_2_position[2];
    return sqrt(x * x + y * y + z * z + SOFTENING);
}

/**
* Calculate the gravitational force between two bodies
*/
inline static double gravitation(double mass_body_1, double mass_body_2, double* body_1_position, double* body_2_position) {
    double top = mass_body_1 * mass_body_2;
    double bottom = euclidean_distance(body_1_position, body_2_position);
    printf("Top / Bottom: %f\n / %f\n", top, bottom);
    printf("bottom * bottom: %f\n", (bottom * bottom));
    printf("(top / bottom * bottom): %f\n", (top / (bottom * bottom)));
    printf("G * (top / bottom * bottom): %f\n", G * (top / (bottom * bottom)));
    return G * (top / (bottom * bottom));
}

/**
* Calculate the net force on a body due to all other bodies
*/
inline static double net_force(double* grav_forces, size_t number_of_grav_forces, double* body_1_position, double* body_2_position) {
    double sum = 0;
    for (int i = 0; i < number_of_grav_forces; i++) {
        sum += grav_forces[i] * ((body_2_position[0] - body_1_position[0]) / euclidean_distance(body_1_position, body_2_position));
    }

    return sum;
}

/**
* Calculate the acceleration of a body due to a force
*/
inline static double get_acceleration(double force, double mass) {
    return force / mass;
}

inline static void save_position(Matrix* output, double* position, size_t output_row, size_t n) {
    for (size_t i = 0; i < n; i++) {
        output->data[output_row + 3*i] = position[3*i];
        output->data[output_row + 3*i+1] = position[3*i+1];
        output->data[output_row + 3*i+2] = position[3*i+2];
        }
}
