#include "functions.h"
#include <map>
#include <chrono>

void ex3();
void ex4();
void ex7();
void ex8();

int main()
{
	ex3();
	ex4();
	ex7();
	ex8();
	return 0;
}

void ex3() {
    std::vector<Complex> test = read("sig1.txt");
    std::vector<Complex> test_dpf(test);
    std::vector<Complex> test_fpf(test);

    std::cout<<"A)"<<std::endl;
    DFT(test_dpf,1);
	write(test_dpf,"DFT");
    DFT(test_dpf,0);
	write(test_dpf,"IDFT");
    std::cout<<"MSE between original and RDFT(DFT(X)) = "<<error(test,test_dpf)<<std::endl;

    std::cout<<"\nB)"<<std::endl;
    FFT(test_fpf,1);
	write(test_fpf,"FFT");
    FFT(test_fpf,0);
	write(test_fpf,"IFFT");
    std::cout<<"MSE between original and RFFT(FFT(X)) = "<<error(test,test_fpf)<<std::endl;

    std::cout<<"\nC)"<<std::endl;
}

void ex4() {
    int N=10;
	std::vector<Complex> time_dft(N);
	std::vector<Complex> time_fft(N);
	std::vector<Complex> dft;
	std::vector<Complex> fft;

	for (int k = 0; k < N; k++) {
		dft.resize(pow(2, k), 1);
        fft.resize(pow(2, k), 1);

		auto start = std::chrono::steady_clock::now();
		DFT(dft, 0);
		auto end = std::chrono::steady_clock::now();
		time_dft[k] = Complex(k,std::chrono::duration<double, std::nano>(end - start).count());

		start = std::chrono::steady_clock::now();
		FFT(fft, 0);
		end = std::chrono::steady_clock::now();
		time_fft[k] =Complex(k,std::chrono::duration<double, std::nano>(end - start).count());
	}
    write(time_dft,"time_dft");
    write(time_fft,"time_fft");
}


void ex7() {
	std::vector<Complex> conv;
	std::vector<Complex> fconv;
	std::vector<Complex> a(16);
	std::vector<Complex> b(16);
	a=read("sig1.txt");
	b=read("sig2.txt");

	conv = Converse(a, b);

    write(conv,"Conv");
	fconv = FastConverse(a, b);
    write(fconv,"FastConv");
}

void ex8() {
    int N = 14;
    std::vector<Complex> time_conv1(N);
    std::vector<Complex> time_fconv1(N);
    std::vector<Complex> time_conv2(N);
    std::vector<Complex> time_fconv2(N);
    std::vector<Complex> a;
    std::vector<Complex> b;
    std::vector<Complex> c1;
    a.resize(256, 1);

    for (int i = 0; i < N; i++) {
        b.resize(pow(2, i), 1);

        auto start = std::chrono::steady_clock::now();
        c1 = Converse(a, b);
        auto end = std::chrono::steady_clock::now();
        time_conv1[i] = Complex(i, std::chrono::duration<double, std::nano>(end - start).count());

        start = std::chrono::steady_clock::now();
        c1 = FastConverse(a, b);
        end = std::chrono::steady_clock::now();
        time_fconv1[i] = Complex(i, std::chrono::duration<double, std::nano>(end - start).count());
    }

    for (int i = 0; i < N; i++) {
        a.resize(pow(2, i), 1);
        b.resize(pow(2, i), 1);

        auto start = std::chrono::steady_clock::now();
        c1 = Converse(a, b);
        auto end = std::chrono::steady_clock::now();
        time_conv2[i] = Complex(i, std::chrono::duration<double, std::nano>(end - start).count());

        start = std::chrono::steady_clock::now();
        c1 = FastConverse(a, b);
        end = std::chrono::steady_clock::now();
        time_fconv2[i] = Complex(i, std::chrono::duration<double, std::nano>(end - start).count());
    }

    write(time_conv1, "time_conv1");
    write(time_conv2, "time_conv2");
    write(time_fconv1, "time_fastconv1");
    write(time_fconv2, "time_fastconv2");

}