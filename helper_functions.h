
inline static void save_position(Matrix* output, double* position, size_t output_row, size_t n) {
    for (size_t i = 0; i < n; i++) {
        output->data[output_row + 3*i] = position[3*i];
        output->data[output_row + 3*i+1] = position[3*i+1];
        output->data[output_row + 3*i+2] = position[3*i+2];
        }
}
