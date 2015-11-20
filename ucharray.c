#include "ucharray.h"

static int newarray (lua_State *L) {
    int i;
    size_t nbytes;
    ucharArray *a;
    const char *fname, *chrtag;
    char line[MAXLEN];

    int n = luaL_checkint(L, 1);
    luaL_argcheck(L, n >= 1, 1, "invalid size");
    nbytes = sizeof(ucharArray) + (n - 1)*sizeof(unsigned char);
    a = (ucharArray *)lua_newuserdata(L, nbytes);
    a->size = n;

    /* luaL_checktype(L, 2, LUA_TSTRING); /\* genome name *\/ */
    /* luaL_checktype(L, 3, LUA_TSTRING); /\* chromosome starting point *\/ */
    
    fname = lua_tostring(L, 2);
    chrtag = lua_tostring(L, 3);
    if (fname == NULL) {
	for (i = 0; i < n; i++)
	    a->values[i] = 0;
    } else {
    	FILE *fp = fopen(fname, "r");
    	if (fp == NULL)
    	    luaL_error(L, "can't open '%s'\n", fname);
    	else {
	    while (fgets(line, MAXLEN, fp) != NULL)
		if (strstr(line, chrtag))
		    break;
	    char c;
	    int k = 0;
	    while ((c = getc(fp)) != '>') 
		if (isalpha(c))
		    a->values[k++] = c;

	    assert(k == n);
    	    /* for (i = 0; i < n; i++) */
    	    /* 	a->values[i] = getc(fp); */
	}
    	fclose(fp);
    }

    luaL_getmetatable(L, "uchar.array");
    lua_setmetatable(L, -2);
    return 1;
}
static int setarray (lua_State *L) {
    ucharArray *a = checkarray(L);
    int index = luaL_checkint(L, 2) - 1;
    int value = luaL_checkint(L, 3);
    luaL_argcheck(L, 0 <= index && index < a->size, 2, "index out of range");
    luaL_argcheck(L, 0 <= value && value < MAX_UCHAR, 3, "invalid value, too large or too small for unsigned char");
    a->values[index] = value;
    return 0;
}
static int getarray (lua_State *L) {
    ucharArray *a = checkarray(L);
    int index = luaL_checkint(L, 2) - 1;
    luaL_argcheck(L, 0 <= index && index < a->size, 2, "index out of range");
    lua_pushinteger(L, a->values[index]);
    return 1;
}
static int getsize (lua_State *L) {
    ucharArray *a = checkarray(L);
    luaL_argcheck(L, a != NULL, 1, "'array' expected");
    lua_pushinteger(L, a->size);
    return 1;
}
static int array2string (lua_State *L) {
    ucharArray *a = checkarray(L);
    lua_pushfstring(L, "array(%d)", a->size);
    return 1;    
}
static const struct luaL_Reg arraylib_m [] = {
    {"__newindex", setarray},
    {"__index", getarray},
    {"__len", getsize},
    {"__tostring", array2string},
    {NULL, NULL}
};
static const struct luaL_Reg arraylib_f [] = {
    {"new", newarray},
    {NULL, NULL}
};
int luaopen_ucharray (lua_State *L) {
    luaL_newmetatable(L, "uchar.array");
    luaL_setfuncs(L, arraylib_m, 0);
    luaL_newlib(L, arraylib_f);
    return 1;
}


