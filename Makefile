CC= gcc
CFLAGS= -g -Wall
CLIB= -lz -lm

OBJ02= seldup.o bamutil.o util.o bgzf.o xorshift1024.o

OBJ06= rmdupsam.o util.o

seldup: $(OBJ02)
	$(CC) $(CFLAGS) -o $@ $(OBJ02) $(CLIB)

rmdupsam: $(OBJ06)
	$(CC) $(CFLAGS) -o $@ $(OBJ06) $(CLIB)

SRC= $(OBJ01:.o=.c)
INC= bamutil.h  bgzf.h  dcpm.h  gffData.h  khash.h  util.h  xorshift1024.h

VERSION=20180815

dist:
	mkdir seldup-dist.$(VERSION)
	cp *.[ch] Makefile LICENSE.md README.md seldup_license.txt seldup-dist.$(VERSION)
	tar czvf seldup-dist.$(VERSION).tgz seldup-dist.$(VERSION)

clean:
	rm -f $(OBJ02) seldup \
	rm -f $(OBJ06) rmdupsam

