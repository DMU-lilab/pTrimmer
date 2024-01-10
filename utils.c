/*************************************************************************
    > File Name: utils.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2024年01月10日 星期三 17时39分16秒
 ************************************************************************/


#include "utils.h"


kstring_t *k_strcpy(kstring_t *kdest, const char *src)
{
    size_t src_len;
    src_len = strlen(src) + 1;

    if (kdest->m < src_len) {
        kdest->m = src_len; kroundup32(kdest->m);
        err_realloc(kdest->s, kdest->m, char);
    }
    kdest->l = src_len-1;
    strcpy(kdest->s, src);
    return kdest;
}


kstring_t *k_strncpy(kstring_t *kdest, const char *src, size_t n)
{
    size_t src_len;
    src_len = n + 1;

    if (kdest->m < src_len) {
        kdest->m = src_len; kroundup32(kdest->m);
        err_realloc(kdest->s, kdest->m, char);
    }
    kdest->l = n;
    memcpy(kdest->s, src, kdest->l*sizeof(char));
    kdest->s[kdest->l] = '\0';
    return kdest;
}


kstring_t *k_strcat(kstring_t *kdest, const char *src)
{
    size_t target_len;
    target_len = strlen(src) + kdest->l + 1;

    if (kdest->m < target_len) {
        kdest->m = target_len; kroundup32(kdest->m);
        err_realloc(kdest->s, kdest->m, char);
    }
    strcpy(kdest->s+kdest->l, src);
    kdest->l = target_len - 1;
    return kdest;
}