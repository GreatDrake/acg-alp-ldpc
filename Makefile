bin:
	@mkdir -p bin

CMD = g++ -std=c++17 -pthread -L . -lglpk -O3 -I.

main.o: bin
	${CMD} -c main.cpp -o bin/main.o

main: main.o
	${CMD} -o bin/main bin/main.o

run: main
	./bin/main 2> log.txt

optimize_H.o: bin
	${CMD} -c optimize_H.cpp -o bin/optimize_H.o

optimize_H: optimize_H.o
	${CMD} -o bin/optimize_H bin/optimize_H.o

optimize: optimize_H
	./bin/optimize_H 2> log.txt

clean:
	rm -rf bin/
