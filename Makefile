bin:
	@mkdir -p bin

CMD = g++ -std=c++17 -pthread -L . -lglpk -O3 -I.
CMD_DEBUG = g++ -std=c++17 -pthread -g -fsanitize=address,undefined
CMD_NOGLPK = g++ -std=c++17 -pthread -O3 -I.

run: bin
	${CMD_NOGLPK} -c main.cpp -o bin/main.o
	${CMD} -o bin/main bin/main.o
	./bin/main 2> log.txt

optimize: bin
	${CMD_NOGLPK} -c optimize_H.cpp -o bin/optimize_H.o
	${CMD_NOGLPK} -o bin/optimize_H bin/optimize_H.o
	./bin/optimize_H 2> log.txt

clean:
	rm -rf bin/

run_qpadmm_params: bin
	${CMD_NOGLPK} -c qpadmm_params.cpp -o bin/qpadmm_params.o
	${CMD_NOGLPK} -o bin/qpadmm_params bin/qpadmm_params.o
	./bin/qpadmm_params 2> log.txt