#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"
#include "thresholds.h"

#define MAXDEPTH 600
#define RA  0 
#define RC  1
#define RG  2
#define RT  3
#define AA  0 
#define AC  4
#define AG  8
#define AT 12
#define SO 16
#define BI 32
#define MO 64
#define REF_MASK   3
#define ALT_MASK  12
#define SOMA_MASK 16
#define BIAL_MASK 32
#define MONO_MASK 64

typedef struct scan_t {
    int beg, end, chr;
    samfile_t *in;
    ucharArray *site;
    ucharArray *info;
} scan_t;
extern int luaopen_readscan (lua_State *L);
