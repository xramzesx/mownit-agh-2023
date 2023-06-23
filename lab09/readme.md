## Laboratorium 9: Rozwiązywanie układów równań liniowych metodami bezpośrednimi

Dany jest układ równań liniowych *Ax=b*.

### Polecenie nr 1

Elementy macierzy *A* o wymiarze $n$x$n$ są określone wzorem:

```math
\begin{cases}
a_{1j} = 1 \\
a_{ij} = = \frac{1}{i+j-1} \ dla\  i \ne 1
\end{cases} \\
i,j = 1, ..., n
```

Przyjmij wektor x jako dowolną n−elementową permutację ze zbioru { 1, -1 } i oblicz wektor b.
Następnie metodą eliminacji Gaussa rozwiąż układ równań liniowych Ax=b (przyjmując jako
niewiadomą wektor x). Przyjmij różną precyzję dla znanych wartości macierzy A i wektora b.
Sprawdź, jak błędy zaokrągleń zaburzają rozwiązanie dla różnych rozmiarów układu (porównaj –
zgodnie z wybraną normą – wektory x obliczony z x zadany). Przeprowadź eksperymenty dla
różnych rozmiarów układu.

### Polecenie nr 2

Powtórz eksperyment dla macierzy zadanej wzorem:
```math
\begin{cases}
a_{ij} = \frac{2i}{j} \ dla j >= i \\
a_{ij} = a_{ji} \ dla\ j < i
\end{cases} \\
i,j = 1, ..., n
```

Porównaj wyniki z tym, co otrzymano w przypadku układu z punktu 1). Spróbuj uzasadnić, skąd
biorą się różnice w wynikach. Sprawdż uwarunkowanie obu układów.

### Polecenie nr 3

Powtórz eksperyment dla jednej z macierzy zadanej wzorem poniżej (macierz i parametry podane
w zadaniu indywidualnym). Następnie rozwiąż układ metodą przeznaczoną do rozwiązywania
układów z macierzą trójdiagonalną. Porównaj wyniki otrzymane dwoma metodami (czas,
dokładność obliczeń i zajętość pamięci) dla różnych rozmiarów układu. Przy porównywaniu
czasów należy pominąć czas tworzenia układu. Opisz, jak w metodzie dla układów z macierzą
trójdiagonalną przechowywano i wykorzystywano macierz A. 

```math
\begin{cases}
a_{i,i} = 5 \\
a_{i,i+1} = \frac{1}{i+m} \\
a_{i,i-1} = = \frac{k}{i+m+1} \\
a_{i,j} = 0 \ dla \ j < - 1 \ oraz \ j > i + 1
\end{cases} \\
i,j = 1, ..., n
```