program:
	gcc -o main main.c -O3 -Wall -lm

vec:
	gcc -o main main.c -O3 -Wall -lm -fopt-info-vec -fopt-info-vec-missed

debug:
	gcc -o main main.c -O3 -Wall -g -Wall -lm

python:
	gcc -fPIC -shared -o ../python_code/supernova.so orbit.h -O3 -lm -fopt-info-vec