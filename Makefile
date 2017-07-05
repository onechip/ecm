CC		= gcc
CXX		= g++
LINK		= g++

NTLPREFIX	= /usr/local
NTLINCLUDE      = -I$(NTLPREFIX)/include
NTLLIB          = -L$(NTLPREFIX)/lib -lntl -lgmp

PREFFLAGS       = -O2 -Wall -I.
CXXFLAGS        = $(PREFFLAGS) $(NTLINCLUDE)
LINKFLAGS       = $(NTLLIB)

ALL_PROGS	= ecm-test factor
COMMON_OBJS	= EC_p.o ZZFactoring.o
COMMON_HEADERS	=
ALL_DIRS	=


all:	$(ALL_DIRS) $(ALL_PROGS)

ecm-test:	$(COMMON_OBJS) ecm-test.o
	$(LINK) -o $@ $^ $(LINKFLAGS)

factor:	$(COMMON_OBJS) factor.o
	$(LINK) -o $@ $^ $(LINKFLAGS)

clean:
	rm -f *% *~ *.o core a.out $(ALL_PROGS)
