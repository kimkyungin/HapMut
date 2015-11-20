typedef struct extract_t {
    int beg, end, chr;
    samfile_t *in;
    ucharArray *info;
    int (*find) (unsigned char);
    lua_State *L;
    int pos_count;
} extract_t;
extern int luaopen_snpextract (lua_State *L);
