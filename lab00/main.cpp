#include <iostream>
#include <cmath>

//// PROTOTYPES ////

template <typename Type>
Type x_n (int n);

template <typename Type>
Type get_x_20 (int precision);

template <typename Type>
Type reverse_x_n (int precision, int n = 0);

//// MAIN ////

int main () {
    int n = 20;

    //// KROK 1 ////

    std::cout << "Wzor rekurencyjny:\n";
    std::cout << "i\tfloat\tdouble\tldouble\n";

    for (int i = 0; i <= n; i++)
        std::cout << i 
            << "\t"<< x_n<float>(i) 
            << "\t"<< x_n<double>(i)
            << "\t"<< x_n<long double>(i) 
            << std::endl;

    std::cout << std::endl;
    std::cout << "x20 : " << x_n<long double>(20) << std::endl;
    std::cout << "x0  : "  << x_n<long double>(0) << std::endl;
    std::cout << std::endl;

    //// KROK 2 i 3 ////

    std::cout << "Odwrotna rekurencja:\n";
    std::cout << "i\tfloat\tdouble\tldouble\n";

    int precision = 20;

    for (int i = 0; i <= n; i++)
        std::cout << i 
            << "\t"<< reverse_x_n<float>(precision, i) 
            << "\t"<< reverse_x_n<double>(precision , i)
            << "\t"<< reverse_x_n<long double>(precision , i) 
            << std::endl;

    std::cout << std::endl;
    std::cout << "x20 : " << get_x_20<long double>(1000) << std::endl;
    std::cout << "x0  : "  << reverse_x_n<long double>(precision) << std::endl;
    
    return 0;
}


//// CODE ////

template <typename Type>
Type x_n (int n) {
    return n > 0 
        ? 1 / ((Type) n) - 5 * x_n<Type>(n - 1)
        : std::log(1.2); 
}

template <typename Type>
Type get_x_20 (int precision) {
    Type result = 0;
    int multipler = 5;
    for (int n = 0; n < precision; n++) {
        result +=  1 / ((Type) multipler * ((Type) n + 20));
        multipler *= - 5;
    }
    return result;
}

template <typename Type>
Type reverse_x_n (int precision, int n)  {
    return n < 20 
        ? 1 / ((Type) n + 1) - reverse_x_n<Type>(precision, n + 1)
        : get_x_20<Type>(precision); 
}