energias.pdf tiempos.pdf: datos.dat tiempos.dat
	python tiempos.py
	python grafica.py
datos.dat tiempos.dat: compila
	./a.out 1 
	./a.out 2
	./a.out 4
compila: FPUT.c
	gcc -fopenmp FPUT.c -lm

clean:
	rm -f datos.dat tiempos.dat energia.pdf tiempos.pdf a.out
