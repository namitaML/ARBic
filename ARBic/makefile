VER=1.0
DIST=ARBic$(VER)
PROGS=ARBic
SRCS=struct.cpp read_array.cpp  make_graph.cpp get_options.cpp fib.cpp write_block.cpp cluster.cpp colbase.cpp main.cpp expand.cpp dag_longest_path.cpp 
OBJS=$(SRCS:.cpp=.o) 
CC=g++


LDFLAGS= -lm -fopenmp
CPPFLAGS= -std=c++0x -g -fopenmp -Wall -I.  -DVER=$(VER)
OMPFLAGS=-fopenmp

all: $(PROGS)

${PROGS}: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

.o:
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f $(PROGS)
	rm -f *.o
	rm -f *.rules
	rm -f *.chars
	rm -f *.blocks
	rm -f *.expansion

dist:
	$(MAKE) clean
	cd .. && tar czvf $(DIST).tar.gz $(DIST)/

Ecoli.chars:
	./${PROGS} -i Ecoli

CRP.blocks:
	./${PROGS} -i CRP 

test: Ecoli.chars CRP.blocks
	./${PROGS} -i Ecoli.chars -b CRP.blocks -s
