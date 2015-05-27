#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const size_t bin_size = sizeof(float) + sizeof(unsigned int);

typedef struct {
    int count;
    float min;
    float max;
    void* bins;
    size_t bin_count;
    size_t max_bin_count;
} histo_t;

inline static float* k(const void* bin) {
    return (float*) bin;
}

inline static unsigned int* v(const void* bin) {
    return (unsigned int*) (bin + sizeof(float));
}

void histo_init(histo_t* this, size_t max_bin_count) {
    this->count         = 0;
    this->min           = 0.0/0.0;
    this->max           = 0.0/0.0;
    this->bins          = malloc((max_bin_count + 1) * bin_size);
    this->bin_count     = 0;
    this->max_bin_count = max_bin_count;
}

void histo_destroy(histo_t* this) {
    free(this->bins);
}

void histo_print(histo_t* this) {
    printf(
        "histo_t( count: %d, min: %f, max: %f, bins: { ",
        this->count,
        this->min,
        this->max);
    for (void* bin = this->bins;
         bin < this->bins + (this->bin_count * bin_size);
         bin += bin_size)
    {
        if (bin != this->bins) printf(", ");
        printf("%f -> %d", *k(bin), *v(bin));
    }
    printf(" } )\n");
}

void histo_insert(histo_t* this, float key) {
    // if the histogram is empty, just insert the new key as a single bin.
    if (this->bin_count == 0) {
        this->count = 1;
        this->min = key;
        this->max = key;
        *k(this->bins) = key;
        *v(this->bins) = 1;
        this->bin_count = 1;
        return;
    }

    // increment the population count, and update the min and/or max as needed.
    this->count += 1;
    if (key < this->min) this->min = key;
    if (key > this->max) this->max = key;

    // find the first bin whose key ≥ the new key.
    void* bin;
    void* overflow_bin = this->bins + this->bin_count * bin_size;
    for (bin = this->bins;
         *k(bin) < key && bin < overflow_bin;
         bin += bin_size);

    // if the current bin's key equals the new key, just increment the bin's
    // value, rather than inserting a new bin.
    if (bin != overflow_bin && *k(bin) == key) {
        *v(bin) += 1;
        return;
    }

    // TODO: only use memmove if we're inserting a new bin (that is, if it
    // doesn't overflow). reimplement the merge case below WITHOUT memmove,
    // merging in place. this also means we can allocate one less bin.

    // at this point, the current bin's key must be greater than the new key.
    // shift this and all subsequent bins to the right, and insert a new bin.
    memmove(bin + bin_size, bin, overflow_bin - bin);
    *k(bin) = key;
    *v(bin) = 1;

    // if insertion didn't overflow, just update the bin_count rather than
    // merging.
    if (this->bin_count < this->max_bin_count) {
        this->bin_count += 1;
        return;
    }

    // find the pair of bins with the closest keys, to merge.
    bin = this->bins;
    float delta = *k(bin + bin_size) - *k(bin);
    for (void* cur_bin = bin + bin_size;
         cur_bin < overflow_bin;
         cur_bin += bin_size)
    {
        float cur_delta = *k(cur_bin + bin_size) - *k(cur_bin);
        if (cur_delta < delta) {
            bin = cur_bin;
            delta = cur_delta;
        }
    }

    // merge these two bins such that the new key is the weighted average of the
    // old keys, and the new value is the sum of the old values.
    float k1 = *k(bin), k2 = *k(bin + bin_size);
    int   v1 = *v(bin), v2 = *v(bin + bin_size);
    *k(bin) = (k1*v1 + k2*v2) / (v1 + v2);
    *v(bin) = v1 + v2;

    // shift all subsequent bins back to the left.
    bin += bin_size;
    memmove(bin, bin + bin_size, overflow_bin - bin);
}

float histo_query(histo_t* this, float quantile) {
    // TODO: when there are exact counts, we should get exact results
    float lhs_count = 0, rhs_count = 0;
    void* bin = this->bins;
    void* overflow_bin = this->bins + this->bin_count * bin_size;
    while (1) {
        float lhs_key, rhs_key;
        if (bin == this->bins) {
            lhs_key = this->min;
            lhs_count = 0;
            rhs_key = *k(bin);
            rhs_count = *v(bin) / 2;
        } else if (bin == overflow_bin) {
            lhs_key = *k(bin - bin_size);
            lhs_count = rhs_count;
            rhs_key = this->max;
            rhs_count = this->count;
        } else {
            lhs_key = *k(bin - bin_size);
            lhs_count = rhs_count;
            rhs_key = *k(bin);
            rhs_count = lhs_count + (float) (*v(bin - bin_size) + *v(bin)) / 2;
        }

        if (rhs_count < quantile * this->count) {
            bin += bin_size;
        } else {
            // TODO: trapezoid, not rectangle
            return lhs_key + (rhs_key - lhs_key) * (quantile * this->count - lhs_count) / (rhs_count - lhs_count);
        }
    }
}

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
