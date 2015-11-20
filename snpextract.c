#include "ucharray.h"
#include "readscan.h"
#include "snpextract.h"

static int find_somatic (unsigned char x) { return x & SOMA_MASK; }
static int find_biallelic (unsigned char x) { return x & BIAL_MASK; }
static int find_monoallelic (unsigned char x) { return x & MONO_MASK; }
static unsigned char decode_ref (unsigned char x);
static unsigned char decode_alt (unsigned char x);
static int snp_extract_c (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
static int snp_extract_depth_c (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
static int snp_extract (lua_State *L) {
    const char *bam = lua_tostring(L, 1);
    int chr = lua_tointeger(L, 2) - 1;
    luaL_argcheck(L, 0 <= chr && chr < 24, 2, "invalid chromosome number");
    ucharArray *info = (ucharArray *)lua_touserdata(L, 3);
    const char *key = lua_tostring(L, 4);

    extract_t read;
    read.beg = 0;
    /* read.end = 0x7fffffff; */
    read.end = info->size - 1;
    read.in = samopen(bam, "rb", 0);
    read.chr = chr;
    read.info = info;
    if (strcmp(key, "somatic") == 0)
	read.find = find_somatic;
    else if (strcmp(key, "biallelic") == 0)
	read.find = find_biallelic;
    else if (strcmp(key, "monoallelic") == 0)
	read.find = find_monoallelic;
    read.L = L;
    read.pos_count = 1;

    lua_newtable(L); /*snp*/
    lua_newtable(L); /*snp_pos*/
    sampileup(read.in, -1, snp_extract_c, &read); 
    
    samclose(read.in);
    return 2; /*return table*/
}
static int snp_extract_c (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
    extract_t *tmp = (extract_t *)data;
    if ((int)pos >= tmp->beg && (int)pos <= tmp->end && (int)tid == tmp->chr) {
	unsigned char value = tmp->info->values[(int)pos];
	if (value > 0 && tmp->find(value)) {
	    bam1_t *ba;
	    bam_pileup1_t plk;
	    unsigned char *seq, *qual;
	    char allele, *qname;
	    int k, qpos, distance;
	    double bcerr;

	    lua_pushinteger(tmp->L, tmp->pos_count); 
	    lua_pushinteger(tmp->L, pos);
	    lua_settable(tmp->L, -3); /* snp_pos[pos_count] = pos */

	    lua_pushvalue(tmp->L, -2); /* switch to snp */

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
		if (allele == decode_ref(value))      
		    distance = 0;
		else if (allele == decode_alt(value)) 
		    distance = 1;
		else                    
		    distance = 2;

		lua_pushstring(tmp->L, qname);
		lua_newtable(tmp->L);
		lua_pushinteger(tmp->L, 1);
		lua_pushinteger(tmp->L, distance);
		lua_settable(tmp->L, -3);
		lua_pushinteger(tmp->L, 2);
		lua_pushnumber(tmp->L, bcerr);
		lua_settable(tmp->L, -3);
		lua_settable(tmp->L, -3);
	    }
	    lua_settable(tmp->L, -3);
	    lua_settop(tmp->L, -2); /*switch to snp_pos*/
	    tmp->pos_count++;
	}
    }
    return 0;
}
static int snp_extract_depth (lua_State *L) {
    const char *bam = lua_tostring(L, 1);
    int chr = lua_tointeger(L, 2) - 1;
    luaL_argcheck(L, 0 <= chr && chr < 24, 2, "invalid chromosome number");
    ucharArray *info = (ucharArray *)lua_touserdata(L, 3);
    const char *key = lua_tostring(L, 4);

    extract_t read;
    read.beg = 0;
    read.end = 0x7fffffff;
    read.in = samopen(bam, "rb", 0);
    read.chr = chr;
    read.info = info;
    if (strcmp(key, "somatic") == 0)
	read.find = find_somatic;
    else if (strcmp(key, "biallelic") == 0)
	read.find = find_biallelic;
    else if (strcmp(key, "monoallelic") == 0)
	read.find = find_monoallelic;
    read.L = L;
    read.pos_count = 1;

    lua_newtable(L); /*snp*/
    lua_newtable(L); /*snp_pos*/
    sampileup(read.in, -1, snp_extract_depth_c, &read); 
    
    samclose(read.in);
    return 2; /*return table*/
}
static int snp_extract_depth_c (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
    extract_t *tmp = (extract_t *)data;
    if ((int)pos >= tmp->beg && (int)pos <= tmp->end && (int)tid == tmp->chr) {
	unsigned char value = tmp->info->values[(int)pos];
	if (value > 0 && tmp->find(value)) {
	    (void) pl;
	    lua_pushinteger(tmp->L, tmp->pos_count); 
	    lua_pushinteger(tmp->L, pos);
	    lua_settable(tmp->L, -3); /* snp_pos[pos_count] = pos */

	    lua_pushvalue(tmp->L, -2); /* switch to snp */

	    lua_pushinteger(tmp->L, pos);
	    lua_pushinteger(tmp->L, n); /* push the read depth instead of true */
	    lua_settable(tmp->L, -3);
	    lua_settop(tmp->L, -2); /*switch to snp_pos*/
	    tmp->pos_count++;
	}
    }
    return 0;
}
static unsigned char decode_ref (unsigned char x) {
    switch (x & REF_MASK) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    }
    return EXIT_FAILURE;
}
static unsigned char decode_alt (unsigned char x) {
    switch (x & ALT_MASK) {
    case 0: return 'A';
    case 4: return 'C';
    case 8: return 'G';
    case 12: return 'T';
    }
    return EXIT_FAILURE;
}
static const struct luaL_Reg snpextract [] = {
    {"extract_depth", snp_extract_depth},
    {"extract", snp_extract},
    {NULL, NULL}
};
int luaopen_snpextract (lua_State *L) {
    luaL_newlib(L, snpextract);
    return 1;
}
