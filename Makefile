CC		= gcc
CPP		= g++
LINK		= g++

NTLPREFIX	= /usr/local
NTLINCLUDE      = -I$(NTLPREFIX)/include
NTLLIB          = -L$(NTLPREFIX)/lib -lntl -lgmp

PREFFLAGS       =
CPPFLAGS        = $(PREFFLAGS) -O2 -Wall $(NTLINCLUDE) -I.
LINKFLAGS       = $(NTLLIB)

ALL_PROGS	= ecm-test
COMMON_OBJS	= EC_p.o pair_ZZ_long.o ZZFactoring.o
COMMON_HEADERS	=
ALL_DIRS	=


.SUFFIXES: .cpp .o

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $*.cpp


all:	$(ALL_DIRS) $(ALL_PROGS)

install:	$(ALL_PROGS)
	cp $(ALL_PROGS) ../../bin

ecm-test:	ecm-test.o $(COMMON_OBJS)
	$(LINK) -o $@ $^ $(LINKFLAGS)

touch:
	touch $(ALL_SOURCE)

clean:
	rm -f *% *~ *.o core a.out $(ALL_PROGS)
