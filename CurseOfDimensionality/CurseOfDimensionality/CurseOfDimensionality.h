#pragma once

#include <iostream>
#include <math.h>
#include <random>
#include <iomanip>
#define ld long double
#define _datos 20000

using namespace std;

struct Point {
	int dimensiones = 0;
	int* coord = 0;
	Point() {};
	~Point();
	void setDimension(int d);
};

ostream& operator <<(ostream &o, const Point& p) {
	o << "(" << setw(4) << p.coord[0];
	for (int i = 1; i < p.dimensiones; i++)
		o << "," << setw(4) << p.coord[i];
	o << ")";
	return o;
}

struct Experimento {
	int const datos = _datos;
	int dimensiones;
	int tabla[11] = {0,0,0,0,0,0,0,0,0,0,0};
	ld *distancias = 0;
	Point* puntos = 0;

	Experimento(int d);
	~Experimento();
	void comenzarExp();
	void generarPts();
	void generarDis();
	void genTabla();

	void printTabla();
	void printPts();
	void printPtsDistance();
};

Point::~Point() {
	if(coord) delete[] coord;
}

void Point::setDimension(int d){
	dimensiones = d;
	coord = new int[d];
}

ld distancia(const Point& p1, const Point& p2) {
	int D = p1.dimensiones;
	ld salida = 0, diferencia;
	for (int i = 0; i < D; i++) {
		diferencia = p1.coord[i] - p2.coord[i];
		salida += diferencia * diferencia;
	}
	return sqrt(salida);
}

void findMinMax(ld* &arr, int n, ld &pmin, ld &pmax) {
	pmax = pmin = arr[0];
	for (int i = 1; i < n; i++) {
		if (arr[i] < pmin) pmin = arr[i];
		else if (arr[i] > pmax) pmax = arr[i];
	}
}

Experimento::Experimento(int d) {
	dimensiones = d;
	distancias = new ld[datos - 1];
}

Experimento::~Experimento(){
	if (puntos) delete[] puntos;
	if (distancias) delete[] distancias;
}

void Experimento::comenzarExp(){
	generarPts();
	generarDis();
	//printPtsDistance();
	genTabla();
	printTabla();
}

void Experimento::generarPts() {
	puntos = new Point[datos];

	random_device rd;  //obtain a seed
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	uniform_int_distribution<> distrib(1, 1000);

	for (int i = 0; i < datos; i++) {
		puntos[i].setDimension(dimensiones);
		int* coord = new int[dimensiones];
		for (int j = 0; j < dimensiones; j++)
			coord[j] = distrib(gen);
		puntos[i].coord = coord;
	}
	//Random number generation based on https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution
}

void Experimento::generarDis(){
	for (int i = 0; i < datos - 1; i++) {
		distancias[i] = distancia(puntos[i], puntos[datos - 1]);
	}
}

void Experimento::genTabla(){
	ld pmin, pmax;
	findMinMax(distancias, datos-1, pmin, pmax);
	for (int i = 0; i < datos-1; i++) {
		ld aux = (distancias[i] - pmin) / (pmax - pmin) *10;//ratio*10
		tabla[int(aux)]++;
	}
}

void Experimento::printTabla(){
	int total = tabla[10];
	for (int i = 0; i < 10; i++) {
		cout << "0." << i << " -> " << tabla[i] << endl;
		total += tabla[i];
	}
	cout << "1.0 -> " << tabla[10] << endl ;
	cout << "Total: " << total << endl;
}

void Experimento::printPts() {
	for (int i = 0; i < datos; i++)
		cout << puntos[i] << endl;
	cout << endl << endl;
}

void Experimento::printPtsDistance() {
	cout << puntos[datos - 1] << endl;
	for (int i = 0; i < datos -1 ; i++)
		cout << puntos[i] << " -> " << distancias[i] << endl;
	cout << endl;
}
