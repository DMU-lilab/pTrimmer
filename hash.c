#include "hash.h"

#define ISHEAD(_node) (!_node->next && !_node->num)

/*! @funciton: get next prime number
 *   @parmeters:
 *   @    size       the prime number next to the given size
 *   @return:
 *   @    tril       a proper prime number
*/
static int NextPrime(int size)
{
    int flag = 1;
    int tril = size;

    while (1) {
        tril++; flag = 1;
        for (int i=2; i < tril; ++i) {
            if (tril % i == 0) {
                flag = 0; break;
            }
        }
        if (flag) 
            return tril;
    }
}


/*! @funciton: calculate the load factor of the hash table
 *   @parmeters:
 *   @    T          the hash table [hast_t *]
 *   @return:
 *   @               the load factor [float]
*/
float LoadFactor(hash_t *T)
{
    int num = 0;
    node_t *node;

    for (int i=0; i < T->size; i++) {
        node = &T->table[i];
        if (ISHEAD(node)) continue;
        do {
            num++;
            node = node->next;
        } while (node);
    }
    return (num / (float)T->size);
}

/*! @funciton: hash function
 *   @parmeters:
 *   @    key        the key for the hash table [const char *]
 *   @    size       the table size [int]
 *   @return:
 *   @               the index in the hash table
*/
int Hash(const char *key, int size)
{
    uint32_t val = 0;

    while (*key != '\0')
        val = (val << 5) -val + *key++;

    return val % size;
}

/*! @funciton: hash table initialization
 *   @parmeters:
 *   @    size       the hash table size [int]
 *   @return:
 *   @    H          the pointer to the hash table [hash_t *]
*/
hash_t *InitHash( int size )
{
    hash_t *H;

    H = (hash_t *)malloc(sizeof(hash_t));
    if (!H)  
        goto _memerror;

    H->size = NextPrime(size);
    H->table = (node_t *)calloc(H->size, sizeof(node_t));
    if (!H->table)
        goto _memerror;

    return H;

  _memerror:
      fprintf(stderr, "[Err::%s] \
          Failed to allocate memory\n", __func__); exit(-1);
}

/*! @funciton: creat a new node
 *   @parmeters:
 *   @    void       
 *   @return:
 *   @    node       the pointer to the new node [node_t *]
*/
static node_t *NewNode(void)
{
    node_t *node;

    node = (node_t *)calloc(1, sizeof(node_t));
    if (!node)
        goto _memerror;
    node->loc = (loc_t *)calloc(1, sizeof(loc_t));
    if (!node->loc)
        goto _memerror;

    return node;

  _memerror:
      fprintf(stderr, "[Err::%s] \
          Failed to allocate memory\n", __func__); exit(-1);

}

/*! @funciton: CORE function for the hash search
 *   @parmeters:
 *   @    key       the key for the hash table [char *]
 *   @    H         the pointer to the hash table [hash_t *]
 *   @    s         the search status [status *]
 *   @return:
 *   @    void
*/
void Search(char *key, hash_t *H, status *s)
{
    node_t *p, *t;

    p = &H->table[Hash(key, H->size)];

    if (ISHEAD(p)) {
        /* the head node is empty */
        s->find = 0; s->node = p;
        return ;
    }
    do {
        t = p;
        if (!strcmp(key, p->key)) {
            s->find = 1; s->node = p;
            return ;
        }
        p = p->next;
    } while (p);

    /* can't find the key in the hash table */
    s->find = 0; s->node = t;
}

/*! @funciton: hash table initialization
 *   @parmeters:
 *   @    size       the hash table size [int]
 *   @return:
 *   @    H          the pointer to the hash table [hash_t *]
*/
void Insert(char *key, loc_t *loc, hash_t *H)
{
    status S;

    Search(key, H, &S);

    if (!S.find && ISHEAD(S.node)) { 
        /* insert the loc info to the head node */
        strcpy(S.node->key, key);
        UpdateLoc(S.node, loc); return ;
    }
    if (!S.find) { // insert the loc info to the new node
        S.node->next = NewNode();
        strcpy(S.node->next->key, key);
        UpdateLoc(S.node->next, loc);
        return ;
    }
    if (S.find) { // just update the loc info
        UpdateLoc(S.node, loc); return ;
    }
}

/*! @funciton: hash table initialization
 *   @parmeters:
 *   @    node      the pointer to the node [node_t *]
 *   @    loc       the kmer who has collided key with others [int]
 *   @return:
 *   @    void
*/
void UpdateLoc(node_t *node, loc_t *loc)
{
    int locnum;

    locnum = node->num;
    if (locnum % LOCNUM == 0) {
        int newsize = locnum + LOCNUM;
        loc_t *tem = (loc_t *)realloc(node->loc, newsize*sizeof(loc_t));
        if (!tem) 
            goto _memerror;
        node->loc = tem;
    }
    node->loc[locnum].ploc = loc->ploc;
    node->loc[locnum].frloc = loc->frloc;
    node->loc[locnum].sloc = loc->sloc;

    node->num++;

    return ;

  _memerror:
      fprintf(stderr, "[Err:%s] Failed to alloc memory\n", __func__);
}


