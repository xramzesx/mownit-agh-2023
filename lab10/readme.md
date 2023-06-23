## Laboratorium 10: Rozwiązywanie układów równań liniowych metodami iteracyjnymi 

Dany jest układ równań liniowych *Ax=b*.
Elementy macierzy A są zadane wzorem:

```math
\begin{cases}
a_{1j} = 8 \\
a_{ij} = = \frac{1}{\left|i-j\right| + 4} \ dla\  i \ne 1
\end{cases} \\
i,j = 1, ..., n
```


Przyjmij wektor *x* jako dowolną $n$-elementową permutację ze zbioru $\left\\{ -1, 1 \right\\}$ i oblicz wektor *b*. 


### Polecenie nr 1

Metodą Jacobiego rozwiąż układ równań liniowych *Ax=b* (przyjmując jako niewiadomą wektor *x*),
przyjmując kolejno kryterium stopu: 

- kryterium przyrostowe:

```math
\left | x^{(i+1)}-x^{(i)} \right | < \rho
```

- kryterium rezydualne:

```math
\left |Ax^{(i)} - b \right | < \rho
```

Obliczenia wykonaj dla różnych rozmiarów układu n, dla różnych wektorów początkowych, a także
różnych wartości $\rho$ w kryteriach stopu. _(Podaj, jak liczono normę.)_ Wyznacz liczbę iteracji oraz
sprawdź różnicę w czasie obliczeń dla obu kryteriów stopu.
Sprawdź dokładność obliczeń. 

### Polecenie nr 2

Dowolną metodą znajdź promień spektralny *macierzy iteracji* (dla różnych rozmiarów układu –
takich, dla których znajdowane były rozwiązania układu). Sprawdź, czy spełnione są założenia o
zbieżności metody dla zadanego układu.

Opisz metodę znajdowania promienia spektralnego. 