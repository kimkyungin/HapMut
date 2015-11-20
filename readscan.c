#include "ucharray.h"
#include "readscan.h"

static int max3 (int x, int y, int z, int *pos);
static int encode_site (int na, int nc, int ng, int nt, char ref, int *conf);
static int decode_site (lua_State *L);
static int read_scan_c (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
static int read_scan (lua_State *L) {
    const char *bam = lua_tostring(L, 1);
    int chr = lua_tointeger(L, 2) - 1;
    luaL_argcheck(L, 0 <= chr && chr < 24, 2, "invalid chromosome number");
    ucharArray *site = (ucharArray *)lua_touserdata(L, 3);
    ucharArray *info = (ucharArray *)lua_touserdata(L, 4);
	
    scan_t read;
    read.beg = 0;
    /* read.end = 0x7fffffff; */
    read.end = info->size - 1;
    read.in = samopen(bam, "rb", 0);
    read.chr = chr;
    read.site = site;
    read.info = info;
    sampileup(read.in, -1, read_scan_c, &read);
    samclose(read.in);
    return 0;
}
static int read_scan_c (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
    scan_t *tmp = (scan_t *)data;
    if ((int)pos >= tmp->beg && (int)pos <= tmp->end && (int)tid == tmp->chr) {
 	char ref = tmp->site->values[(int)pos];
	if (ref == 'A' || ref == 'C' || ref == 'G' || ref == 'T') {
	    bam1_t *ba;
	    bam_pileup1_t plk;
	    unsigned char *seq, *qual;
	    char allele;
	    int k, qpos, qsum = 0;
	    int n_alt = 0, na = 0, nc = 0, ng = 0, nt = 0, nn = 0, nm = 0, nr = 0;
	    for (k = 0; k < n; k++) {
		plk = pl[k];
		ba = plk.b;
		qpos = plk.qpos;
		seq = bam1_seq(ba);
		qual = bam1_qual(ba);
		qsum += (int)qual[qpos];
		allele = bam_nt16_rev_table[bam1_seqi(seq, qpos)];
		switch (allele) {
		case 'A': {na++; break;}
		case 'C': {nc++; break;}
		case 'G': {ng++; break;}
		case 'T': {nt++; break;}
		case 'N': {nn++; break;}
		case 'M': {nm++; break;}
		case 'R': {nr++; break;}
		default: {
		    fprintf(stderr, "Unknown DNA type: '%c'\n", allele);
		    return EXIT_FAILURE;
		}
		}
	    }
	    int allele_conf = 0;
	    n_alt = encode_site(na, nc, ng, nt, ref, &allele_conf);
	    
	    /* if nn > 0 and alignment is correct then the site should have many N */
	    if (9 < n && n < MAXDEPTH && nn == 0 && nm == 0 && nr == 0) {
		int cind = n - 10;
		qsum /= n;
		/* if qsum < 10 or qsum > 60, then discard the site */
		if (10 <= qsum && qsum <= 60) {
		    /* if qsum == 20, then 11th row is the position, which is rind = 10 */
		    int rind = qsum - 10; 
		    
		    if (germ_err[rind][cind] <= n_alt && n_alt <= bial_err[rind][cind])
			tmp->info->values[pos] += BI;
		    else if (n_alt > bial_err[rind][cind])
			tmp->info->values[pos] += MO;
		    if (seq_err[rind][cind] < n_alt)
		        tmp->info->values[pos] += SO;
		    /* else if (seq_err[rind][cind] < n_alt)  */
		    /* 	tmp->info->values[pos] += SO; */
		    tmp->info->values[pos] += allele_conf;
		}
	    }
	}
    }
    return 0;
}
static int max3 (int x, int y, int z, int *pos) {
    int max = x, imax = 1;
    if (y > max) {max = y; imax = 2;}
    if (z > max) {max = z; imax = 3;}
    *pos = imax;
    return max;
}
static int encode_site (int na, int nc, int ng, int nt, char ref, int *conf) {
    int nalt = 0, pos = 0;
    switch (ref) {
    case 'A': {
	nalt = max3(nc, ng, nt, &pos);
	if      (pos == 1) *conf = RA + AC;
	else if (pos == 2) *conf = RA + AG;
	else if (pos == 3) *conf = RA + AT;
	else
	    fprintf(stderr, "Position error\n");
	return nalt;
    }
    case 'C': {
	nalt = max3(na, ng, nt, &pos);
	if      (pos == 1) *conf = RC + AA;
	else if (pos == 2) *conf = RC + AG;
	else if (pos == 3) *conf = RC + AT;
	else
	    fprintf(stderr, "Position error\n");
	return nalt;
    }
    case 'G': {
	nalt = max3(na, nc, nt, &pos);
	if      (pos == 1) *conf = RG + AA;
	else if (pos == 2) *conf = RG + AC;
	else if (pos == 3) *conf = RG + AT;
	else
	    fprintf(stderr, "Position error\n");
	return nalt;
    }
    case 'T': {
	nalt = max3(na, nc, ng, &pos);
	if      (pos == 1) *conf = RT + AA;
	else if (pos == 2) *conf = RT + AC;
	else if (pos == 3) *conf = RT + AG;
	else
	    fprintf(stderr, "Position error\n");
	return nalt;
    }
    default: {
	return EXIT_FAILURE;
    }
    }
}
static int decode_site (lua_State *L) { 
    ucharArray *info = (ucharArray *)lua_touserdata(L, 1);
    int n = info->size;
    unsigned char *value = info->values;

    lua_newtable(L);
    int ts = lua_gettop(L);
    lua_newtable(L);
    int tb = lua_gettop(L);
    lua_newtable(L);
    int tm = lua_gettop(L);
	
    int i;
    for (i = 0; i < n; i++) {
	unsigned char v = value[i];
	if ((SOMA_MASK & v) > 0) {
	    if ((BIAL_MASK & v) > 0) { 
		lua_pushinteger(L, i);
		lua_pushboolean(L, 1);
		lua_settable(L, ts);
		lua_pushinteger(L, i);
		lua_pushboolean(L, 1);
		lua_settable(L, tb);
	    } else if ((MONO_MASK & v) > 0) {
		lua_pushinteger(L, i);
		lua_pushboolean(L, 1);
		lua_settable(L, ts);
		lua_pushinteger(L, i);
		lua_pushboolean(L, 1);
		lua_settable(L, tm);
	    } else {
		lua_pushinteger(L, i);
		lua_pushboolean(L, 1);
		lua_settable(L, ts);
	    }
	} else {
	    if ((BIAL_MASK & v) > 0) { 
		lua_pushinteger(L, i);
		lua_pushboolean(L, 1);
		lua_settable(L, tb);
	    } else if ((MONO_MASK & v) > 0) {
		lua_pushinteger(L, i);
		lua_pushboolean(L, 1);
		lua_settable(L, tm);
	    } 
	}
    }
    return 3;
}
static const struct luaL_Reg readscan [] = {
    {"scan", read_scan},
    {"decode", decode_site},
    {NULL, NULL}
};
int luaopen_readscan (lua_State *L) {
    luaL_newlib(L, readscan);
    return 1;
}

