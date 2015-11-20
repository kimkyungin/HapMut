#include "blocklib.h"

static double alpha_cutoff = 0.01;
enum { REF_HAP = 0, ALT_HAP = 1 };
enum { SAME = 0, DIFFERENT = 1, WRONG = 2 };
static double parallel_prob (lua_State *L, int tsnp, int tsnp_pos, int k) {
    lua_pushinteger(L, k);
    lua_gettable(L, tsnp_pos); /* snp_pos[k] */
    lua_gettable(L, tsnp); /* snp[snp_pos[k]] == a */
    int ta = lua_gettop(L);

    lua_pushinteger(L, k + 1);
    lua_gettable(L, tsnp_pos); /* snp_pos[k+1] */
    lua_gettable(L, tsnp); /* snp[snp_pos[k+1]] == b */
    int tb = lua_gettop(L);
    
    double hpr = 1.0;

    const char *key;
    lua_pushnil(L);
    while (lua_next(L, ta) != 0) {
	key = lua_tostring(L, -2); /* S: key, a[k] */
	int tak = lua_gettop(L);

	lua_pushstring(L, key);
	lua_gettable(L, tb);
	int tbk = lua_gettop(L); /* S: key, a[k], b[k] */

	if (!lua_isnil(L, tbk)) {
	    double parallel = 0.0, crossed = 0.0;
	    lua_pushinteger(L, 1);
	    lua_gettable(L, tak);
	    int v = lua_tointeger(L, -1); /* a[k][1] */
	    lua_pushinteger(L, 2);
	    lua_gettable(L, tak);
	    double x = lua_tonumber(L, -1); /* a[k][2] */
	    lua_pushinteger(L, 1);
	    lua_gettable(L, tbk);
	    int w = lua_tointeger(L, -1); /* b[k][1] */
	    lua_pushinteger(L, 2);
	    lua_gettable(L, tbk);
	    double y = lua_tonumber(L, -1); /* b[k][2] */
	    lua_pop(L, 4);

	    if (v == 0) {
		if (w == 0) {
		    parallel = (1 - x)*(1 - y) + x/3*y/3;
		    crossed  = (1 - x)*y/3     + x/3*(1 - y);
		} else if (w == 1) {
		    parallel = (1 - x)*y/3     + x/3*(1 - y);
		    crossed  = (1 - x)*(1 - y) + x/3*y/3;
		} else if (w == 2) {
		    parallel = (1 - x)*y*2/3   + x/3*y*2/3;
		    crossed  = (1 - x)*y*2/3   + x/3*y*2/3;
		}
	    } else if (v == 1) {
		if (w == 0) {
		    parallel = x/3*(1 - y)   + (1 - x)*y/3;
		    crossed  = x/3*y/3       + (1 - x)*(1 - y);
		} else if (w == 1) {
		    parallel = x/3*y/3       + (1 - x)*(1 - y);
		    crossed  = x/3*(1 - y)   + (1 - x)*y/3;
		} else if (w == 2) {
		    parallel = x/3*y*2/3     + (1 - x)*y*2/3;
		    crossed  = x/3*y*2/3     + (1 - x)*y*2/3;
		}
	    } else if (v == 2) {
		if (w == 0) {
		    parallel = x*2/3*(1 - y) + x*2/3*y/3;
		    crossed  = x*2/3*y/3     + x*2/3*(1 - y);
		} else if (w == 1) {
		    parallel = x*2/3*y/3     + x*2/3*(1 - y);
		    crossed  = x*2/3*(1 - y) + x*2/3*y/3;
		} else if (w == 2) {
		    parallel = x*2/3*y*2/3   + x*2/3*y*2/3;
		    crossed  = x*2/3*y*2/3   + x*2/3*y*2/3;
		}
	    }
	    hpr *= crossed/parallel;
	}
	lua_pop(L, 2); /* S: key */
    }
    lua_pop(L, 2); /* pop a, b */
    return 1./(1. + hpr);
}
static int make_hapblock (lua_State *L) {
    luaL_checktype(L, 1, LUA_TTABLE); /* snp */
    luaL_checktype(L, 2, LUA_TTABLE); /* snp_pos */
    lua_newtable(L); /* hapblock */
    int tsnp = 1, tsnp_pos = 2, thapblock = 3;
    int nsnp = luaL_len(L, tsnp_pos);
    int nhapblock = 0;

    int k = 1;
    while (k < nsnp) {
	lua_newtable(L); /* u */
	lua_pushinteger(L, 1);
	lua_pushinteger(L, k);
	lua_gettable(L, tsnp_pos);
	lua_settable(L, -3);
	int tu = lua_gettop(L);
	int nu = luaL_len(L, tu);

	lua_newtable(L); /* v */
	lua_pushinteger(L, 1);
	lua_pushinteger(L, REF_HAP);
	lua_settable(L, -3);
	int tv = lua_gettop(L);
	int nv = luaL_len(L, tv);
	
	while (k < nsnp) {
	    double p = parallel_prob(L, tsnp, tsnp_pos, k);
	    int parallel = p > 1 - alpha_cutoff;
	    int crossed = p < alpha_cutoff;
	    if (!(parallel || crossed)) 
		break;

	    lua_pushinteger(L, ++nu);
	    lua_pushinteger(L, k + 1);
	    lua_gettable(L, tsnp_pos);
	    lua_settable(L, tu);

	    lua_pushinteger(L, nv);
	    lua_gettable(L, tv);
	    int vv = lua_tointeger(L, -1);
	    lua_pop(L, 1);
	    if ((parallel && (vv == REF_HAP)) || (crossed && (vv == ALT_HAP))) {
		lua_pushinteger(L, ++nv);
		lua_pushinteger(L, REF_HAP);
		lua_settable(L, tv);
	    } else if ((parallel && (vv == ALT_HAP)) || (crossed && (vv == REF_HAP))) {
		lua_pushinteger(L, ++nv);
		lua_pushinteger(L, ALT_HAP);
		lua_settable(L, tv);
	    }
	    k++;
	}
	lua_pushinteger(L, ++nhapblock);
	lua_newtable(L);
	lua_pushstring(L, "pos");
	lua_pushvalue(L, tu);
	lua_settable(L, -3);
	lua_pushstring(L, "hap1");
	lua_pushvalue(L, tv);
	lua_settable(L, -3);
	lua_settable(L, thapblock);

	lua_pop(L, 2); /* pop u, v */
	k++;
    }
    return 1;
}
static int make_blockread (lua_State *L) {
    luaL_checktype(L, 1, LUA_TTABLE); /* snp */
    luaL_checktype(L, 2, LUA_TTABLE); /* hapblock */
    int tsnp = 1, thapblock = 2;
    int nhapblock = luaL_len(L, thapblock);

    lua_newtable(L);
    int tblockread = lua_gettop(L);
    
    int i;
    for (i = 1; i <= nhapblock; i++) {
	lua_pushinteger(L, i);
	lua_gettable(L, thapblock); /* hapblock[i] */
	int ihapblock = lua_gettop(L);
	lua_pushstring(L, "pos");
	lua_gettable(L, ihapblock); /* hapblock[i].pos */
	int ipos = lua_gettop(L);
	int npos = luaL_len(L, ipos);
	lua_pushstring(L, "hap1");
	lua_gettable(L, ihapblock); /* hapblock[i].hap1 */
	int ihap1 = lua_gettop(L);
	
	lua_newtable(L); /* read */
	int tread = lua_gettop(L); 

	int j;
	for (j = 1; j <= npos; j++) {
	    lua_pushinteger(L, j);
	    lua_gettable(L, ipos);
	    int x = (int) lua_tonumber(L, -1);

	    lua_gettable(L, tsnp); /* snp[x], #S = 7 */
	    int snp_x = lua_gettop(L);

	    const char *name;
	    int val, hap1j;
	    lua_pushnil(L);
	    while (lua_next(L, snp_x) != 0) {
		name = lua_tostring(L, -2);
		lua_pushinteger(L, 1);
		lua_gettable(L, -2); /* snp_x[name][1] */
		val = (int) lua_tonumber(L, -1);
		lua_pop(L, 2); /* pop snp_x[name] snp_x[name][1] */
		
		lua_pushinteger(L, j);
		lua_gettable(L, ihap1);
		hap1j = (int) lua_tonumber(L, -1);

		lua_pushstring(L, name);
		lua_gettable(L, tread); /* read[name] */
		if (lua_isnil(L, -1)) {
		    lua_pop(L, 1);
		    lua_pushstring(L, name);
		    lua_newtable(L);
		    lua_settable(L, tread);
		    lua_pushstring(L, name);
		    lua_gettable(L, tread);
		} /* S: name hap1j read[name] */

		lua_pushinteger(L, x); /* S: name hap1j read[name] x */
		if (val == 2) {
		    lua_pushinteger(L, WRONG);
		    lua_settable(L, -3);
		} else if (val != hap1j) {
		    lua_pushinteger(L, DIFFERENT);
		    lua_settable(L, -3);
		} else if (val == hap1j) {
		    lua_pushinteger(L, SAME);
		    lua_settable(L, -3);
		} /* S: name hap1j read[name] */
		lua_pop(L, 2); /* S: ... name */
	    }
	    lua_pop(L, 1); /* pop snp[x] */
	}

	lua_pushnil(L);
	while (lua_next(L, tread) != 0) {
	    const char *y = lua_tostring(L, -2);
	    double ratio = 1.0;
	    lua_pushnil(L);
	    while (lua_next(L, -2) != 0) {
		int x = (int) lua_tonumber(L, -2);
		int w = (int) lua_tonumber(L, -1);
		lua_pushinteger(L, x);
		lua_gettable(L, tsnp); /* snp[x] */
		lua_pushstring(L, y);
		lua_gettable(L, -2); /* snp[x] snp[x][y] */
		lua_pushinteger(L, 2);
		lua_gettable(L, -2); /* snp[x] snp[x][y] snp[x][y][2] */
		double eps = lua_tonumber(L, -1);
		if (w == SAME)
		    ratio *= eps/(1 - eps);
		else if (w == DIFFERENT)
		    ratio *= (1 - eps)/eps;
		else if (w == WRONG)
		    ratio *= 1.;
		lua_pop(L, 4); /* remain x */		
	    }
	    lua_pushstring(L, "prob1");
	    lua_pushnumber(L, 1./(1. + ratio));
	    lua_settable(L, -3);
	    lua_pop(L, 1);	    
	}
	lua_pushinteger(L, i);
	lua_newtable(L);

	lua_pushstring(L, "read");
	lua_pushvalue(L, tread);
	lua_settable(L, -3);

	lua_pushstring(L, "begin");
	lua_pushinteger(L, 1);
	lua_gettable(L, ipos);
	lua_settable(L, -3);

	lua_pushstring(L, "einde");
	lua_pushinteger(L, npos);
	lua_gettable(L, ipos);
	lua_settable(L, -3);
	
	lua_settable(L, tblockread); /* blockread[i] */
	
	lua_pop(L, 4); /* pop hapblock[i] hapblock[i].pos hapblock[i].hap1 read */
    }
    return 1;
}
static int select_block (lua_State *L, int tb, int k, int t) {
    lua_pushinteger(L, k);
    lua_gettable(L, tb); /* block[k] */
    lua_pushstring(L, "read");
    lua_gettable(L, -2); /* block[k].read */
    int left = lua_gettop(L);
    lua_pushinteger(L, k + 1);
    lua_gettable(L, tb); /* block[k + 1] */
    lua_pushstring(L, "read");
    lua_gettable(L, -2); /* block[k + 1].read */
    int right = lua_gettop(L);

    /* S: block[k] block[k].read block[k+1] block[k+1].read */

    int lcount = 0, rcount = 0;
    int i;
    const char *name;
    for (i = 1; i <= luaL_len(L, t); i++) {
	lua_pushinteger(L, i);
	lua_gettable(L, t); /* S: t[i] */
	lua_pushstring(L, "name");
	lua_gettable(L, -2); /* S: t[i] t[i].name */
	name = lua_tostring(L, -1);

	lua_pushstring(L, name);
	lua_gettable(L, left);
	if (lua_isnil(L, -1)) {
	    lua_pop(L, 1);
	} else {
	    lcount++;
	    lua_pop(L, 1);
	}

	lua_pushstring(L, name);
	lua_gettable(L, right);
	if (lua_isnil(L, -1)) {
	    lua_pop(L, 1);
	} else {
	    rcount++;
	    lua_pop(L, 1);
	}
	lua_pop(L, 2);
    }
    lua_pop(L, 4);
    if (rcount > lcount)
	return k + 1;
    else
	return k;
}
static void assign (lua_State *L, int t, int r) {
    const char *name;
    int i;
    for (i = 1; i <= luaL_len(L, t) ; i++) {
	lua_pushinteger(L, i);
	lua_gettable(L, t); /* S: t[i] */
	int ti = lua_gettop(L);
	lua_pushstring(L, "name");
	lua_gettable(L, ti); /* S: t[i] t[i].name */
	name = lua_tostring(L, -1);
	lua_pop(L, 1); /* S: t[i] */

	lua_pushstring(L, name);
	lua_gettable(L, r); /* S: t[i] r[name] */
	if (lua_isnil(L, -1)) 
	    lua_pop(L, 2); /* pop t[i] + nil */
	else {
	    lua_pushstring(L, "prob1");
	    lua_gettable(L, -2); /* S: t[i] r[name] r[name].prob1 */
	    double p = lua_tonumber(L, -1);
	    lua_pushstring(L, "hp");
	    lua_pushnumber(L, p);
	    lua_settable(L, ti);
	    lua_pop(L, 3); /* pop t[i] + r[name] + r[name].prob1 */
	}
    }
}
static int getread (lua_State *L, int x, int tb) {
    lua_pushinteger(L, x);
    lua_gettable(L, tb); /* block[x] */
    lua_pushstring(L, "read");
    lua_gettable(L, -2); /* block[x] block[x].read */
    int r = lua_gettop(L); 
    return r;
}
static int assign_prob (lua_State *L) {
    luaL_checktype(L, 1, LUA_TTABLE); /* out */
    luaL_checktype(L, 2, LUA_TTABLE); /* out_pos */
    luaL_checktype(L, 3, LUA_TTABLE); /* block */
    int tout = 1, toutpos = 2, tblock = 3;
    int nout = luaL_len(L, toutpos);
    int nblock = luaL_len(L, tblock);

    int *begin = (int *)malloc((nblock + 1)*sizeof(int));
    int *einde = (int *)malloc((nblock + 1)*sizeof(int));
    begin[0] = -1; /* unused entry, index start from 1 */
    einde[0] = -1; /* unused entry, index start from 1 */
    
    int i;
    for (i = 1; i <= nblock; i++) { /* S: out block */
    	lua_pushinteger(L, i);
    	lua_gettable(L, tblock); /* block[i] */
    	lua_pushstring(L, "begin");
    	lua_gettable(L, -2); /* block[i] block[i].begin */
    	begin[i] = (int) lua_tonumber(L, -1);
    	lua_pop(L, 1); /* S[-1] = block[i] */
	
    	lua_pushstring(L, "einde");
    	lua_gettable(L, -2); /* block[i] block[i].einde */
    	einde[i] = (int) lua_tonumber(L, -1);
    	lua_pop(L, 2);
    }

    /* S: out out_pos block */
    assert(lua_gettop(L) == 3);

    int blockpos = 0; 
    for (i = 1; i <= nout; i++) {
	lua_pushinteger(L, i);
	lua_gettable(L, toutpos); 
	int pos = (int) lua_tonumber(L, -1);
	lua_pushinteger(L, pos);
	lua_gettable(L, tout);
	int t = lua_gettop(L);
	/* int nt = luaL_len(L, t); */
	
	/* S: out out_pos block out_pos[i] out[out_pos[i]] */
	assert(t == 5); 

	while (blockpos <= nblock) {
	    if (0 < blockpos && blockpos < nblock) {
		if (begin[blockpos] < pos && pos < begin[blockpos + 1])
		    break;
	    } else if (blockpos == 0) {
		if (pos < begin[1])
		    break;
	    } else if (blockpos == nblock) {
		if (pos > begin[nblock])
		    break;
	    }
	    blockpos++;
	}

	if (0 < blockpos && blockpos < nblock) {
	    if (begin[blockpos] < pos && pos < einde[blockpos]) {
		assign(L, t, getread(L, blockpos, tblock));
		lua_pop(L, 2); /* getread increas two stacks */
	    } else if (einde[blockpos] < pos && pos < begin[blockpos + 1]) {
		int x = select_block(L, tblock, blockpos, t);
		assign(L, t, getread(L, x, tblock));
		lua_pop(L, 2);
	    }
	} else if (blockpos == 0) {
	    assign(L, t, getread(L, 1, tblock));
	    lua_pop(L, 2);
	} else if (blockpos == nblock) {
	    assign(L, t, getread(L, nblock, tblock));
	    lua_pop(L, 2);
	}

	lua_pop(L, 2);
    }
    free(begin);
    free(einde);
    return 0;
}
static const struct luaL_Reg blocklib [] = {
    {"mkhapblock", make_hapblock},
    {"mkblockread", make_blockread},
    {"assignpr", assign_prob},
    {NULL, NULL}
};
int luaopen_blocklib (lua_State *L) {
    luaL_newlib(L, blocklib);
    return 1;
}


