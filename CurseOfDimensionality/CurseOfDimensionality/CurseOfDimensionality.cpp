#include "CurseOfDimensionality.h"

using namespace std;

int main() {
	const int experimentos = 6;
	int dimensiones[experimentos] = { 2,5,10,15,20,25 };
	for (int i = 0; i < experimentos; i++) {
		cout << endl << dimensiones[i] << " dimensiones" << endl;
		Experimento exp(dimensiones[i]);
		exp.comenzarExp();
	}

	return 0;
}
