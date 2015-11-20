OS= $(shell uname)

CC= gcc

# Please modify the following directory information
# samtools-0.1.19 and lua 5.2 header and library directory
ISAMTOOLS= ${HOME}/Downloads/samtools-0.1.19
LSAMTOOLS= ${HOME}/Downloads/samtools-0.1.19
ILUA= /usr/local/include 
LLUA= /usr/local/lib

ifeq ($(OS), Darwin)
CFLAGS= $(IDIR) -Wall -Wextra -O2 
LDFLAGS= $(LDIR) -bundle -undefined dynamic_lookup 
LIBS= -lbam -lz
endif

ifeq ($(OS), Linux)
CFLAGS= $(IDIR) -Wall -Wextra -O2 -fPIC 
LDFLAGS= $(LDIR) -shared
LIBS= -lbam -lz -lpthread
endif

IDIR= -I$(ILUA) -I$(ISAMTOOLS)
LDIR= -L$(LLUA) -L$(LSAMTOOLS)


SO= ucharray.so readscan.so snpextract.so readextract.so blocklib.so mutlib.so
OBJ= ucharray.o readscan.o snpextract.o readextract.o blocklib.o mutlib.o
HDR= thresholds.h ucharray.h readscan.h snpextract.h readextract.h blocklib.h mutlib.h

.SUFFIXES: .so .o .c

all: $(SO) thresholds.h

$(OBJ): $(HDR)

thresholds.h: ./generate-thresholds.R
	Rscript ./generate-thresholds.R 0.5

.o.so: 
	gcc $(LDFLAGS) -o $*.so $*.o $(LIBS)

.c.o:
	gcc $(CFLAGS) -c $*.c -o $*.o

.PHONY: clean
clean:
	rm -f *.o *~

