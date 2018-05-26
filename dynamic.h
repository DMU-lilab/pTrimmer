/*
Copyright (c) 2010,2011 Marcel Martin <marcel.martin@tu-dortmund.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum FLAG {
    START_WITHIN_SEQ1 = 1,
    START_WITHIN_SEQ2 = 2,
    STOP_WITHIN_SEQ1 = 4,
    STOP_WITHIN_SEQ2 = 8,
    SEMIGLOBAL = 15,
    ALLOW_WILDCARD_SEQ1 = 1,
    ALLOW_WILDCARD_SEQ2 = 2,
};

#define GAPCHAR '\0'

#define DELETION_COST 1
#define INSERTION_COST 1
#define MISMATCH_COST 1
#define MATCH_COST 0

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

typedef struct {
    int start1, stop1;
    int start2, stop2;
    int matches;
    int errors;
} dynamic_t;


/*! @typedef amp_t
 @abstract structure for the DP matrix
 @field cost           the insert length for the paired primer
 @field matches        number of matches in this alignment
 @field fwdprim        where the alignment originated: negative for positions 
                       withinseq1, positive for position within seq2
*/
typedef struct {
    int cost;
    int matches;
    int origin;
} Entry;

dynamic_t global_align(const char *, int, const char *, int, int);
