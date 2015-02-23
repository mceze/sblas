TOPDIR = $(PWD)

INCDIR = $(TOPDIR)/include
SRCDIR = $(TOPDIR)/src
LIBDIR = $(TOPDIR)/lib
CC = clang
LD = clang
CFLAGS = -g -O3 -fPIC
LDFLAGS = -dynamiclib
LIBS = -lm

OBJLIST = $(SRCDIR)/sblas_io.o \
          $(SRCDIR)/sblas_linalg.o \
          $(SRCDIR)/sblas_utils.o \
          $(SRCDIR)/sblas_linsolv.o \
          $(SRCDIR)/sblas_aux.o

DEPS = $(INCDIR)/%.h

all: sblasobj libsblas.dylib

sblasobj:
	@(cd $(SRCDIR); $(CC) -c $(SRCDIR)/*.c $< $(CFLAGS) $(LIBS) -I$(INCDIR))


libsblas.dylib: $(OBJLIST)
	$(LD) $(LDFLAGS) -o $@ $(CFLAGS) $^ $(LIBS)
	mv $@ $(LIBDIR)/

.PHONY: clean

clean:
	rm -f $(SRCDIR)/*.o *~ core $(INCDIR)/*~ $(LIBDIR)/*
