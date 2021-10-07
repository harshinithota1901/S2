#
# GNU makefile, application portion.
#
# "make" or "make all" to build executable.
# "make clean" to delete object code.
#
# some users have reported a "high performance arithmetic trap" on DEC Alphas
# and on LINUX.  if you get this, try adding "-DTOLERANCE_ZERO" to CFLAGS --
# it apparently works on Alphas; don't know about Linux.
#

KERNELDIR = ../../kernel
CC = gcc
CFLAGS = -O2
TARGET = gp

uobjects = function.o app.o
uheaders = appdef.h app.h function.h

include $(KERNELDIR)/GNUmakefile.kernel
