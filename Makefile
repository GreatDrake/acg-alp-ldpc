compile_bp:
	g++ -std=c++17 -O3 bp.cpp -pthread -o bp

run_bp: compile_bp
	./bp

compile_acgalp:
	g++ -std=c++17 -O3 acgalp.cpp -lglpk -L . -o acgalp

run_acgalp: compile_acgalp
	./acgalp
