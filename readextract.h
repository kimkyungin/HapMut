typedef struct read_extract_t {
    int beg, end, chr;
    samfile_t *in;
    unsigned *list;
    int nlist;
    unsigned char *ref;
    lua_State *L;
    int pos_count;
} read_extract_t;
extern int luaopen_readextract (lua_State *L);
