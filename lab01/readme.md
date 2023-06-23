## Laboratorium 1: Arytmetyka komputerowa 

Niech

```math
x_n=\int_0^1 t^n(t+5)^{-1}dt
```

Całka powyższa może być obliczona za pomocą wzoru rekurencyjnego:

```math
x_n=\frac{1}{n}-5x_{n-1}
```

przy czym

```math
x_0 = ln6 - ln5 = ln1.2
```

Stosując ten wzór rekurencyjny obliczyć $x_0, x_1, ..., x_{20}$. Zwrócić uwagę, czy i dla
jakich $n$ obliczony wynik jest ujemny. Czy tak powinno być? 

Obliczyć analitycznie $x_{20}$ wykorzystując do funkcji podcałkowej szereg Taylora, a
następnie na podstawie tego wzoru wyliczyć w komputerze wartość $x_{20}$. Po
wyznaczeniu $x_{20}$ zastosować rekurencję wstecz wykorzystując odpowiednio
przekształcony wór (*), tj. obliczyć kolejno $x_{19}, x_{18}, … , x_0$. Czy otrzymane $x_0$ jest
poprawne? Wykonać obliczenia dla różnej precyzji zmiennych (float, double, long
double). Skomentować i spróbować objaśnić otrzymane wyniki. 