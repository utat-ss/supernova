program:
	gcc -o main main.c physics.c butcher.c solvers.c orbit.c -O3 -Wall -lm

python:
	gcc -fPIC -shared -o ../python_code/supernova.so physics.c butcher.c solvers.c orbit.c -O3 -lm