IDIR =../include
CC=clang
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

BINDIR=../bin
APPNAME=mom

LIBS=-lm

_DEPS=libmath.h libtriangles.h
DEPS= $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ= mom.o libmath.o libtriangles.o 
OBJ= $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
		$(CC) -c -g -o $@ $< $(CFLAGS)

$(APPNAME): $(OBJ)
		$(CC) -o $(BINDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
		rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(BINDIR)/$(APPNAME) 
