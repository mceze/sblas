TOPDIR = /Users/mdebarro/Desktop/sblas

INCDIR = $(TOPDIR)/include
SRCDIR = $(TOPDIR)/src
LIBDIR = $(TOPDIR)/lib
CC = ~/homebrew/bin/gcc-4.9
LD = ~/homebrew/bin/gcc-4.9
CFLAGS = -g -O3 -fPIC -arch i386
LDFLAGS = -shared -arch i386
LIBS = -fopenmp

OBJLIST = $(SRCDIR)/sblas_io.o \
          $(SRCDIR)/sblas_linalg.o \
          $(SRCDIR)/sblas_utils.o \
          $(SRCDIR)/sblas_linsolv.o \
          $(SRCDIR)/sblas_aux.o

DEPS = $(INCDIR)/%.h

all: sblasobj libsblas

sblasobj:
	@(cd $(SRCDIR); $(CC) -c $(SRCDIR)/*.c $< $(CFLAGS) $(LIBS) -I$(INCDIR))


libsblas: $(OBJLIST)
	$(LD) $(LDFLAGS) -o $@.dylib $^ $(LIBS)
	mv $@.dylib $(LIBDIR)/

.PHONY: clean

clean:
	rm -f $(SRCDIR)/*.o *~ core $(INCDIR)/*~ $(LIBDIR)/*
