#pragma once

#include <stddef.h>

typedef struct {
    int count;
    float min;
    float max;
    void* bins;
    size_t bin_count;
    size_t max_bin_count;
} histo_t;

void histo_init(histo_t*, size_t max_bin_count);
void histo_destroy(histo_t*);

void histo_print(histo_t*);

void histo_insert(histo_t*, float key);
float histo_query(histo_t*, float quantile);
