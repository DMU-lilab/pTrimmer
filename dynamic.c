#include "dynamic.h"

dynamic_t global_align(const char *s1, int m, const char *s2, int n, int error)
{
    int flags = SEMIGLOBAL;
    int degenerate = 0;
    dynamic_t result = {0};

    /*
    DP Matrix:
               s2 (j)
            ----------> n
           |
    s1 (i) |
           |
           V
           m
    */

    Entry* column; /* only a single column of the DP matrix is stored */
    column = (Entry*)malloc((m+1)*sizeof(Entry));
    if (column == NULL)
        return result;

    int i, j, best_i, best_j, best_cost, best_matches, best_origin;

    //initialize first column
    for (i = 0; i <= m; ++i) {
        column[i].matches = 0;
        column[i].cost = (flags & START_WITHIN_SEQ1) ? 0 : i * DELETION_COST;
        column[i].origin = (flags & START_WITHIN_SEQ1) ? -i : 0;
    }

    best_i = m;
    best_j = 0;
    best_cost = column[m].cost;
    best_matches = 0;
    best_origin = column[m].origin;

    // maximum no. of errors
    int k = error;
    int last = k + 1;
    if (flags & START_WITHIN_SEQ1) {
        last = m;
    }
    // iterate over columns
    for (j = 1; j <= n; ++j) {
        // remember first entry
        Entry tmp_entry = column[0];

        // fill in first entry in this column TODO move out of loop
        if (flags & START_WITHIN_SEQ2) {
            column[0].cost = 0;
            column[0].origin = j;
            column[0].matches = 0;
        } else {
            column[0].cost = j * INSERTION_COST;
            column[0].origin = 0;
            column[0].matches = 0;
        }
        for (i = 1; i <= last; ++i) {
            int match = (s1[i-1] == s2[j-1])
                        || ((degenerate & ALLOW_WILDCARD_SEQ1) && (s1[i-1] == 'N'))
                        || ((degenerate & ALLOW_WILDCARD_SEQ2) && (s2[j-1] == 'N'));
            int cost_diag = tmp_entry.cost + (match ? MATCH_COST : MISMATCH_COST);
            int cost_deletion = column[i].cost + DELETION_COST;
            int cost_insertion = column[i-1].cost + INSERTION_COST;

            int origin, cost, matches;
            if (cost_diag <= cost_deletion && cost_diag <= cost_insertion) {
                // MATCH or MISMATCH
                cost = cost_diag;
                origin = tmp_entry.origin;
                matches = tmp_entry.matches + match;
            } else if (cost_insertion <= cost_deletion) {
                // INSERTION
                cost = cost_insertion;
                origin = column[i-1].origin;
                matches = column[i-1].matches;
            } else {
                // DELETION
                cost = cost_deletion;
                origin = column[i].origin;
                matches = column[i].matches;
            }

            // remember current cell for next iteration
            tmp_entry = column[i];

            column[i].cost = cost;
            column[i].origin = origin;
            column[i].matches = matches;

        }
        while (column[last].cost > k) {
            last--;
        }
        if (last < m) {
            last++;
        } else {
            // found
            // if requested, find best match in last row
            if (flags & STOP_WITHIN_SEQ2) {
                // length of the aligned part of string1
                int length = m + min(column[m].origin, 0);
                int cost = column[m].cost;
                int matches = column[m].matches;
                if (cost <= error && (matches > best_matches || (matches == best_matches && cost < best_cost))) {
                    // update
                    best_matches = matches;
                    best_cost = cost;
                    best_origin = column[m].origin;
                    best_i = m;
                    best_j = j;
                }
            }

        }
        // column finished
    }

    if (flags & STOP_WITHIN_SEQ1) {
        // search in last column // TODO last?
        for (i = 0; i <= m; ++i) {
            int length = i + min(column[i].origin, 0);
            int cost = column[i].cost;
            int matches = column[i].matches;
            if (cost <= error && (matches > best_matches || (matches == best_matches && cost < best_cost))) {
                // update best
                best_matches = matches;
                best_cost = cost;
                best_origin = column[i].origin;
                best_i = i;
                best_j = n;
            }
        }
    } free(column);

    if (best_origin >= 0) {
        result.start1 = 0;
        result.start2 = best_origin;
    } else {
        result.start1 = -best_origin;
        result.start2 = 0;
    }
    result.stop1 = best_i; result.stop2 =best_j;
    result.matches = best_matches; 
    result.errors = best_cost;

#ifdef _DYNAMIC_MAIN
    fprintf(stderr, "start1:%d\tstop1:%d\nstart2:%d\tstop2:%d\nmatches:%d\nerrors:%d\n", result.start1, result.stop1, result.start2, result.stop2, result.matches, result.errors);
#endif
    return result;
}

#ifdef _DYNAMIC_MAIN
int main(int argc, char **argv)
{
    global_align(argv[1], strlen(argv[1]), argv[2], strlen(argv[2]), 3);

    return 0;

}
#endif

