#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#define PI 3.14159265358979323846

using Complex = std::complex<double>;

std::vector<Complex> read(std::string filename);
void write(std::vector<Complex> x,std::string filename);
double error(std::vector<Complex> x,std::vector<Complex> y);
void DFT(std::vector<Complex>& x, bool flag);
void FFT(std::vector<Complex>& x, bool flag);
std::vector<Complex> Converse(std::vector<Complex> x, std::vector<Complex> y);
std::vector<Complex> FastConverse(std::vector<Complex> x, std::vector<Complex> y);