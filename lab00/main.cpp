#include <iostream>
#include <cmath>

template <typename T>
T x_n ( int n ) {
    return n <= 0 ? std::log(1.2) : 1 / ( (T) n ) - 5 * x_n<T>(n - 1);
}



int main () {
    int n = 20;
    std::cout << "Wzor rekurencyjny:\n";
    std::cout << n << "\t" << "wynik\n";
    
    for (int i = 0; i <= n; i++)
        std::cout << i 
            << "\t"<< x_n<float>(i) 
            << "\t"<< x_n<double>(i)
            << "\t"<< x_n<long double>(i) 
            << std::endl;

    

    // std::cout << x_n<int>(n) << std::endl;
    // std::cout << x_n<float>(n) << std::endl;
    // std::cout << x_n<double>(n) << std::endl;
    // std::cout << x_n<long double>(n) << std::endl;
    // return 0;z
}