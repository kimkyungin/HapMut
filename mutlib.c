#include "mutlib.h"

static int max16 (double *x);
static void baf (lua_State *L, int t, double *x1, double *x2) {
    int n = luaL_len(L, t), nalt = 0;
    double h1 = 0.0, a1 = 0.0;

    int i;
    for (i = 1; i <= n; i++) {
	lua_pushinteger(L, i);
	lua_gettable(L, t); /* t[i] */
	lua_pushstring(L, "dist");
	lua_gettable(L, -2);
	int dist = (int) lua_tonumber(L, -1);
	lua_pop(L, 1);
	lua_pushstring(L, "hp");
	lua_gettable(L, -2);
	double hp = lua_tonumber(L, -1);
	lua_pop(L, 1);
	if (dist != 0) {
	    nalt++;
	    a1 += hp;
	}
	h1 += hp;
	lua_pop(L, 1);
    }

    *x1 = a1/h1;
    *x2 = (nalt - a1)/(n - h1);
    return;	
}
static int posterior (lua_State *L) {
    luaL_checktype(L, 1, LUA_TTABLE); /* t */
    int tt = 1; /* stack position of t */
    int n = luaL_len(L, tt);
    lua_newtable(L); /* s */
    int ts = 2;

    double M = 6e-6; /* mutation rate */
    double G[3] = { (1 - M)*(1 - M), M*(1 - M)/3, M*M/9 };
    
    double x1, x2;
    baf(L, tt, &x1, &x2);

    double s[16], t[16], r[16];
    int i;
    for (i = 0; i < 16; i++) {
	s[i] = 0.0;
	t[i] = 1.0;
	r[i] = 1.0;
    }

    for (i = 1; i <= n; i++) {
	lua_pushinteger(L, i);
	lua_gettable(L, tt);
	lua_pushstring(L, "hp");
	lua_gettable(L, -2);
	double h1 = lua_tonumber(L, -1);
	double h2 = 1. - h1;
	lua_pop(L, 1);
	lua_pushstring(L, "eps");
	lua_gettable(L, -2);
	double q = lua_tonumber(L, -1);
	lua_pop(L, 1);
	lua_pushstring(L, "dist");
	lua_gettable(L, -2);
	int dist = (int) lua_tonumber(L, -1);
	lua_pop(L, 2);

	if (dist == 0) {
	    double y   = 1. - q;
	    double y1  = x1*q/3. + (1. - x1)*(1. - q);
	    double y2  = x2*q/3. + (1. - x2)*(1. - q);
	    double y01 = y *h1 + y2*h2;
	    double y10 = y1*h1 + y *h2;
	    double y11 = y1*h1 + y2*h2;
	    s[0] = y;
	    s[1] = s[2] = s[3]  = y01;
	    s[4] = s[8] = s[12] = y10;
	    s[5] = s[6] = s[7]  = s[9] = s[10] = s[11] = s[13] = s[14] = s[15] = y11; 
	} else if (dist == 1) {
	    double y   = q/3.;
	    double y1  = x1*(1. - q) + (1. - x1)*q/3.;
	    double y2  = x2*(1. - q) + (1. - x2)*q/3.;
	    double y11 = y1*h1 + y2*h2;
	    double y10 = y1*h1 + y *h2;
	    double y01 = y *h1 + y2*h2;
	    s[5] = y11;
	    s[4] = s[6] = s[7]  = y10;
	    s[1] = s[9] = s[13] = y01;
	    s[0] = s[2] = s[3]  = s[8] = s[10] = s[11] = s[12] = s[14] = s[15] = y;
	} else if (dist == 2) {
	    double y   = q/3.;
	    double y1  = x1*(1. - q) + (1. - x1)*q/3.;
	    double y2  = x2*(1. - q) + (1. - x2)*q/3.;
	    double y22 = y1*h1 + y2*h2;
	    double y20 = y1*h1 + y *h2;
	    double y02 = y *h1 + y2*h2;
	    s[10] = y22;
	    s[8] = s[9] = s[11] = y20;
	    s[2] = s[6] = s[14] = y02;
	    s[0] = s[1] = s[3]  = s[4] = s[5] = s[7] = s[12] = s[13] = s[15] = y;
	} else if (dist == 3) {
	    double y   = q/3.;
	    double y1  = x1*(1. - q) + (1. - x1)*q/3.;
	    double y2  = x2*(1. - q) + (1. - x2)*q/3.;
	    double y33 = y1*h1 + y2*h2;
	    double y30 = y1*h1 + y *h2;
	    double y03 = y *h1 + y2*h2;
	    s[15] = y33;
	    s[12] = s[13] = s[14] = y30;
	    s[3]  = s[7]  = s[11] = y03;
	    s[0]  = s[1]  = s[2]  = s[4] = s[5] = s[6] = s[8] = s[9] = s[10] = y;
	} else {
	    luaL_error(L, "Unknown DNA distance: %d\n", dist);
	}
	int j;
	for (j = 0; j < 16; j++)
	    t[j] *= s[j];
    }
    r[0]  = t[0] *G[0]; r[1]  = t[1] *G[1]; r[2]  = t[2] *G[1]; r[3]  = t[3] *G[1];
    r[4]  = t[4] *G[1]; r[5]  = t[5] *G[2]; r[6]  = t[6] *G[2]; r[7]  = t[7] *G[2];
    r[8]  = t[8] *G[1]; r[9]  = t[9] *G[2]; r[10] = t[10]*G[2]; r[11] = t[11]*G[2];
    r[12] = t[12]*G[1]; r[13] = t[13]*G[2]; r[14] = t[14]*G[2]; r[15] = t[15]*G[2];

    int k = max16(r);
    if (k == 1 || k == 2 || k == 3 || k == 4 || k == 8 || k == 12) {    
	double sum = 0.0;
	for (i = 0; i < 16; i++)
	    sum += r[i];
	
	lua_pushstring(L, "prob");
	lua_newtable(L);
	for (i = 0; i < 16; i++) {
	    lua_pushinteger(L, i + 1);
	    /* lua_pushnumber(L, r[i]/sum); */
	    lua_pushnumber(L, round(1000*r[i]/sum)/1000);
	    lua_settable(L, -3);
	}
	lua_settable(L, ts);
	
	lua_pushstring(L, "imax");
	if (k == 1 || k == 4) 
	    lua_pushnumber(L, 1);
	else if (k == 2 || k == 8)
	    lua_pushinteger(L, 2);
	else
	    lua_pushinteger(L, 3);
	lua_settable(L, ts);
    } else {
	lua_pushboolean(L, 0);
    }
    return 1;
}
static int max16 (double *x) {
    double max = 0.0;
    int i, imax = -1;
    for (i = 0; i < 16; i++)
	if (x[i] > max) {
	    max = x[i];
	    imax = i;
	}
    return imax;
}
static const struct luaL_Reg mutlib [] = {
    {"posterior", posterior},
    {NULL, NULL}
};
int luaopen_mutlib (lua_State *L) {
    luaL_newlib(L, mutlib);
    return 1;
}



