/*
 *  author: Suhas Vittal
 *  date:   7 August 2023
 *
 *  Syndromes should be in the format:
 *
 *  <bit0><bit1><bit2>[...].<obs>
 *  i.e.
 *          00000010010101.0
 *          00011101010111.1
 * */

#include <stdio.h>

#include "uf_common.h"

int main(int argc, char* argv[]) {
    if (argc != 1) {
        fprintf(stderr, "Please provide a file containing syndromes and observables.\n");
        return 1;
    }

    const char* syndrome_file = argv[1];

    FILE* fin = fopen(syndrome_file, "r");
    int awaiting_obs = 0;

    uint64_t syndrome[PACKED_INPUT_DATA_ARRAY_SIZE];
    uint64_t offset = 0;

    uint64_t corr = 0;

    uint64_t logical_errors = 0;
    uint64_t number_of_syndromes = 0;

    while (true) {
        char x = fgetc(fin);
        if (x == '.') {
            awaiting_obs = 1;
        } else if (x == '\n') {
            awaiting_obs = 0;
            offset = 0;
            /* DECODE HERE */
            decode(syndrome, &corr);
            if (corr != 0)  logical_errors++;
            number_of_syndromes++;
        } else if (x == '0' || x == '1') {
            if (awaiting_obs) {
                corr = x - '0';
            } else {
                uint64_t major = offset >> 6;
                uint64_t minor = offset & 0x3f;
                syndrome[major][minor] = x - '0';
            }
        }
    }
    double logical_error_rate = ((double) logical_errors) / ((double) number_of_syndromes);

    printf("Logical error rate: %.3e\n", logical_error_rate);

    fclose(fin);
}

