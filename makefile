all:
	g++ -std=c++11 -O2 -Wall -fsanitize=address -fsanitize=leak -fsanitize=undefined -ggdb 3d.cpp Automata3D.cpp Colisionador.cpp Cuerpo.cpp LB3D.cpp Random.cpp Vector.cpp
