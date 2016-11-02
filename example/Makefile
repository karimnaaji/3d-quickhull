
all: qh.out quickhull_c99.out random_generation.out

clean:
	rm -f *.out hull_*

quickhull_c99.out: quickhull.c ../quickhull.h
	gcc -std=c99 $< -Wall -Wno-null-dereference -o $@ -I..

qh.out: command_line.cpp ../quickhull.h
	g++ $< -O3 -o $@ -I.. -I../tests -std=c++14

random_generation.out: random_generation.cpp ../quickhull.h ../quickhull_debug.h
	g++ $< -O3 -o $@ -I.. -Wall -Wno-null-dereference -Wno-unused-function -std=c++14

run-bench: qh.out
	./bench.sh
	cat log.txt

