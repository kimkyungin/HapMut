#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <ctype.h>
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#define MAX_UCHAR 256
#define MAXLEN 256
#define checkarray(L) (ucharArray *)luaL_checkudata(L, 1, "uchar.array")

typedef struct ucharArray {
    int size;
    unsigned char values[1];
} ucharArray;
extern int luaopen_ucharray (lua_State *L);


