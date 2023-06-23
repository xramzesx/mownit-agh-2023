## Laboratorium 8: Rozwiązywanie równań nieliniowych

### Polecenie

Stosując metodę Newtona oraz metodę siecznych wyznacz pierwiastki równania
$f(x)=0$ w zadanym przedziale $[a, b]$. Dla metody Newtona wybierz punkty startowe rozpoczynając
od wartości końców przedziału, zmniejszając je o $0.1$ w kolejnych eksperymentach
numerycznych. Odpowiednio dla metody siecznej jeden z końców przedziału stanowić powinna
wartość punktu startowego dla metody Newtona, a drugi – początek, a następnie koniec przedziału
$[a, b]$.

Porównaj liczbę iteracji dla obu tych metod (dla różnych dokładności $\rho$), stosując jako kryterium
stopu: 

- kryterium przyrostowe:

```math
\left | x^{(i+1)}-x^{(i)} \right | < \rho
```

- kryterium rezydualne:

```math
\left |f(x^{(i)}) \right | < \rho
```

### Badana funkcja

```math
f(x) = x^2 - 20 sin^{12}(x)
```
```math
\frac{df}{dx}=2x - 240sin^{11}(x)cos(x)
```