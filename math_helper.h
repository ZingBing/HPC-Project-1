/**
 * Declares all of the math helper functions (which are defined in math_helper.c).
 */

/**
* Calculate the euclidean distance between two bodies in 3D space
*/
double euclidean_distance(double* body_1_position, double* body_2_position);

/**
* Calculate the gravitational force between two bodies
*/
double gravitation(double mass_body_1, double mass_body_2, double* body_1_position, double* body_2_position);