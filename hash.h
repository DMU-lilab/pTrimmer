/* The hash function is based on separate chaining */


#ifndef HASH_H
#define HASH_H

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#define KEYLEN 128 // maxmum length for hash key
#define LOCNUM 8   // initial number for location index


/*! @typedef loc_t
 @abstract structure for the location index
 @field ploc         primer index in the primelist
 @field frlo         forward or revrse index [0=>forward and 1=>reverse]
 @field sloc         the index of key in the primer sequence
*/
typedef struct __loc_t {
    int ploc;
    int frloc;
    int sloc;
} loc_t;


/*! @typedef node_t
 @abstract structure for the node in hash table
 @field num          the number in the loc
 @field loc          the list for the collided primer who have same key
 @field key          the key for the hash table
 @field next         next node for key who have same hash value
*/
typedef struct __node_t node_t;
struct __node_t {
    int num;
    char key[KEYLEN];
    loc_t *loc;
    node_t *next;
};


/*! @typedef hash_t
 @abstract structure for the hash table
 @field size         hash table size
 @field table        root pointer for the hash table
*/
typedef struct __hash_t {
    uint32_t size;
    node_t *table;
} hash_t;


/*! @typedef status_t
 @abstract structure for the status return by search function
 @field size         find[0=>not found  and 1=>found]
 @field table        the node adress for the serch
 @NOTE: 
    if (find == 0):
        node => where could insert the new node
    else (find == 1):
        node => where the loc_t should be appended on
*/
typedef struct __status {
    int find;
    node_t *node;
} status;


/* prototype function */
int Hash(const char *, int); // hash function
float LoadFactor(hash_t *); // calculate the loadfactor for the hash table
hash_t *InitHash(int); // hash table initialization
void Search(char *, hash_t *, status *); // search the key loaction in the hash table
void Insert(char *, loc_t *, hash_t *); // insert a new key-value pair into hash table
void UpdateLoc(node_t *, loc_t *); // update the loc of primer index

#endif
