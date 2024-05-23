#include "functions.h"

std::vector<Complex> read(std::string filename) {
    std::vector<Complex> y;
    std::ifstream in(filename); 
	filename = "../pythoncode/" + filename;
    if (in.is_open()){
        double a, b;
        while (in >> a >> b)
            y.push_back((a + b * Complex(0.0, 1.0)));
    }
    else{std::cout<<"Error while reading"<<std::endl;}
    in.close();
    return y;
}

void write(std::vector<Complex> x,std::string filename){
    std::ofstream out;
    std::string name = filename + ".txt";
	out.open(name);
	if (out.is_open()) {
		for (int i = 0; i < x.size(); i++)
			out << std::setprecision(20) << x[i].real() << " " << std::setprecision(20) << x[i].imag() << std::endl;
	}
	else {std::cout<<"Error while writing"<<std::endl;}
	out.close();
}

double error(std::vector<Complex> x,std::vector<Complex> y){
    double err=0;
    for(int i=0;i<x.size();i++)
        err+=abs(pow(x[i]-y[i],2));
    return err;
}

void DFT(std::vector<Complex> &x, bool flag) {
	int N = x.size();
	std::vector<Complex> y(N);
	Complex sum(0, 0);

	for (int k = 0; k < N; ++k) {
		sum = 0;
		for (int j = 0; j < N; ++j) {
			sum += x[j] * Complex(cos(2 * M_PI * k * j / N), pow(-1, flag) * sin(2 * M_PI * k * j / N));
		}
		y[k] = sum / sqrt(N);
	}
	x.swap(y);
}

void FFT(std::vector<Complex> &x, bool flag) {
	int N = x.size();
	int n = log2(N);
	Complex h;
	Complex y1;
	int ind; //индексы для перестановок
	int new_ind;	

	//перестановка
	for (int i = 1; i < N; ++i) {
		ind = i;
		new_ind = 0;
		for (int j = 0; j < n; ++j)
		{
			new_ind = (new_ind << 1) | (ind & 1);
			ind >>= 1;
		}
		if(new_ind > i) swap(x[i], x[new_ind]);
	}
		
	//сам бпф
	for (int k = 1; k < n + 1; ++k) {
		for (int j = 0; j < pow(2, n - k); ++j) {
			for (int l = 0; l < pow(2, k - 1); ++l) {
				h = Complex(cos(2 * M_PI * l / pow(2, k)), pow(-1, flag) * sin(2 * M_PI * l / pow(2, k))) * x[j * pow(2, k) + l + pow(2, k - 1)]; //w_n^k*y1(k)
				y1 = x[j * pow(2, k) + l] + h;
				x[j * pow(2, k) + l + pow(2, k - 1)] = (x[j * pow(2, k) + l] - h) / sqrt(2);
				x[j * pow(2, k) + l] = y1 / sqrt(2);
			}
		}
	}
}

std::vector<Complex> Converse(std::vector<Complex> x, std::vector<Complex> y) {
	if (x.size() > y.size()) std::swap(x, y);
	int n = x.size();
	int m = y.size();
	std::vector<Complex> result(m + n);
	Complex sum;
	for (int i = 0; i < m; i++) {
		sum = 0;
		for (int j = 0; j < i + 1; j++) {
			sum += x[j] * y[i - j];
		}
		result[i] = sum;
	}
	for (int i = m; i < result.size(); i++) {
		sum = 0;
		for (int j = 0; j < m; j++) {
			if (i - j < m) sum += x[j] * y[i - j];
		}
		result[i] = sum;
	}
	return result;
}

std::vector<Complex> FastConverse(std::vector<Complex> x, std::vector<Complex> y) {
	Complex zero = 0;
	int m = x.size();
	int l = y.size();
	int n = m > l ? 2 * m : 2 * l;
	int N = pow(2, ceil(log2(n)));

	std::vector<Complex> c(N);
	x.resize(N, zero);
	y.resize(N, zero);
	FFT(x, 0);
	FFT(y, 0);
	double s = sqrt(N);
	for (int i = 0; i < N; i++) {
		c[i] = s * x[i] * y[i];
	}
	FFT(c, 1);
	return c;
}