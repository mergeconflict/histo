#include "histo.h"

#include <stdio.h>

void histo_lol(histo_t* this) {
    printf(
        "min: %f, p50: %f, p90: %f, p95: %f, p99: %f, p995: %f, p999: %f, p9995: %f, p9999: %f, max: %f\n",
        histo_query(this, 0),
        histo_query(this, 0.5),
        histo_query(this, 0.9),
        histo_query(this, 0.95),
        histo_query(this, 0.99),
        histo_query(this, 0.995),
        histo_query(this, 0.999),
        histo_query(this, 0.9995),
        histo_query(this, 0.9999),
        histo_query(this, 1));
}

int main(int argc, char** argv) {
    histo_t histo;
    histo_init(&histo, 40);
    while (1) {
        float key;
        scanf("%f", &key);
        histo_insert(&histo, key);
        histo_print(&histo);
        histo_lol(&histo);
    }
    return 0;
}
