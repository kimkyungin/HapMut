#include "ucharray.h"
#include "readscan.h"
#include "readextract.h"

static int index_found = 0;
static int in_list (unsigned pos, unsigned *list, int nlist);
static int distance (unsigned char ref, unsigned char all);
static int read_extract_c (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
static int read_extract (lua_State *L) {
    const char *bam = lua_tostring(L, 1);
    int chr = lua_tointeger(L, 2) - 1;
    luaL_argcheck(L, 0 <= chr && chr < 24, 2, "invalid chromosome number");
    luaL_checktype(L, 3, LUA_TTABLE); /* list should be sorted as passed in */
    int tlist = 3;
    ucharArray *a = (ucharArray *)luaL_checkudata(L, 4, "uchar.array");
    unsigned char *ref = a->values;

    int len = luaL_len(L, tlist); 
    unsigned *list = (unsigned *)malloc(len*sizeof(unsigned));
    int i;
    for (i = 0; i < len; i++) {
    	lua_pushinteger(L, i + 1);
    	lua_rawget(L, tlist);
    	list[i] = lua_tointeger(L, -1);
    	lua_pop(L, 1);
    }

    read_extract_t read;
    read.beg = 0;
    /* read.end = 0x7fffffff; */
    read.end = a->size - 1;
    read.in = samopen(bam, "rb", 0);
    read.chr = chr;
    read.list = list;
    read.nlist = len;
    read.ref = ref;
    read.L = L;
    read.pos_count = 1;

    lua_newtable(L); /* out_pos */
    lua_newtable(L); /* out */
    sampileup(read.in, -1, read_extract_c, &read); 
    lua_pushvalue(L, -2); /* exchange out_pos and out */
    
    free(list);
    samclose(read.in);
    return 2; /*return table*/
}
static int read_extract_c (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
    read_extract_t *tmp = (read_extract_t *)data;
    if ((int)pos >= tmp->beg && (int)pos <= tmp->end && (int)tid == tmp->chr) {
	if (in_list((int)pos, tmp->list, tmp->nlist)) {
	    unsigned char r = tmp->ref[(int)pos];
	    bam1_t *ba;
	    bam_pileup1_t plk;
	    unsigned char *seq, *qual;
	    char allele, *qname;
	    int k, qpos;
	    double bcerr;

	    lua_pushinteger(tmp->L, pos);
	    lua_newtable(tmp->L);
	    for (k = 0; k < n; k++) {
		plk = pl[k];
		qpos = plk.qpos;
		ba = plk.b;
		qname = bam1_qname(ba);
		seq = bam1_seq(ba);
		qual = bam1_qual(ba);
		bcerr = phred2error[(int)qual[qpos]];
		allele = bam_nt16_rev_table[bam1_seqi(seq, qpos)];

		lua_pushinteger(tmp->L, k + 1);
		lua_newtable(tmp->L);
		
		lua_pushstring(tmp->L, "name");
		lua_pushstring(tmp->L, qname);
		lua_settable(tmp->L, -3);

		lua_pushstring(tmp->L, "eps");
		lua_pushnumber(tmp->L, bcerr);
		lua_settable(tmp->L, -3);

		lua_pushstring(tmp->L, "dist");
		lua_pushinteger(tmp->L, distance(r, allele));
		lua_settable(tmp->L, -3);

		lua_pushstring(tmp->L, "dna");
		lua_pushfstring(tmp->L, "%c", allele);
		lua_settable(tmp->L, -3);

		lua_pushstring(tmp->L, "hp");
		lua_pushnumber(tmp->L, 0.5);
		lua_settable(tmp->L, -3);

		lua_settable(tmp->L, -3);		    
	    }
	    lua_settable(tmp->L, -3);  /* out[pos] */
	    lua_pushinteger(tmp->L, tmp->pos_count++);
	    lua_pushinteger(tmp->L, pos);
	    lua_settable(tmp->L, -4); /* out_pos[pos_count] = pos */
	}
    }
    return 0;
}
static int in_list(unsigned pos, unsigned *list, int n) {
    int i;
    for (i = index_found; i < n; i++)
	if (pos == list[i]) { /*zero-based*/
	    index_found = i;
	    return 1;
	}
    return 0;
}
static int distance (unsigned char ref, unsigned char all) {
    if (ref == 'A')
	switch (all) {
	case 'A': return 0;
	case 'C': return 1;
	case 'G': return 2;
	case 'T': return 3;
	default: return 4;
	}
    else if (ref == 'C')
	switch (all) {
	case 'A': return 3;
	case 'C': return 0;
	case 'G': return 1;
	case 'T': return 2;
	default: return 4;
	}
    else if (ref == 'G')
	switch (all) {
	case 'A': return 2;
	case 'C': return 3;
	case 'G': return 0;
	case 'T': return 1;
	default: return 4;
	}
    else if (ref == 'T')
	switch (all) {
	case 'A': return 1;
	case 'C': return 2;
	case 'G': return 3;
	case 'T': return 0;
	default: return 4;
	}
    else
	return 4;
}
static const struct luaL_Reg readextract [] = {
    {"extract", read_extract},
    {NULL, NULL}
};
int luaopen_readextract (lua_State *L) {
    luaL_newlib(L, readextract);
    return 1;
}

