program:
	gcc -o main main.c physics.c solvers.c orbit.c -O3 -Wall -lm

vec:
	gcc -o main main.c physics.c solvers.c orbit.c -O3 -Wall -lm -fopt-info-vec

debug:
	gcc -o main main.c physics.c solvers.c orbit.c -O3 -Wall -g -Wall -lm

python:
	gcc -fPIC -shared -o ../python_code/supernova.so physics.c solvers.c orbit.c -O3 -lm -fopt-info-vec