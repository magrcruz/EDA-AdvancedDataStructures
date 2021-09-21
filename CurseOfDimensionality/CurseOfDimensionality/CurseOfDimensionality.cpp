#include "CurseOfDimensionality.h"
#include <chrono>

using namespace std;

int main() {
	const int experimentos = 6;
	int dimensiones[experimentos] = { 2,5,10,15,20,25 };
	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	for (int i = 0; i < experimentos; i++) {
		cout << endl<<"Experimento con " << dimensiones[i] << " dimensiones" << endl;
		start = std::chrono::high_resolution_clock::now();
		Experimento exp(dimensiones[i]);
		exp.comenzarExp();
		end = std::chrono::high_resolution_clock::now();
		int64_t duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		cout << "Tiempo empleado en el experimento: " << duration << "ms" << endl;
	}

	//El tiempo incluye instrucciones como new

	return 0;
}
