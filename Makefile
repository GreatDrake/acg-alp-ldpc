bin:
	@mkdir -p bin

CMD = g++ -std=c++17 -pthread -L . -lglpk -O3 -I.

main.o: bin
	${CMD} -c main.cpp -o bin/main.o

main: main.o
	${CMD} -o bin/main bin/main.o

clean:
	rm -rf bin/

run: main
	./bin/main 2> log.txt
