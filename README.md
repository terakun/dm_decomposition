# dm_decomposition
DM分解のC++実装
### How to Use

data.txt
```
7 7
0 1 0 1 2 0 0
0 0 2 1 0 5 0
0 0 0 4 0 0 1
1 0 0 1 0 0 0
3 1 0 0 1 1 0
0 0 0 3 0 0 1
0 0 1 0 0 3 4
```

```
$ ./main data.txt
original matrix:
  0  1  0  1  2  0  0
  0  0  2  1  0  5  0
  0  0  0  4  0  0  1
  1  0  0  1  0  0  0
  3  1  0  0  1  1  0
  0  0  0  3  0  0  1
  0  0  1  0  0  3  4
strongly connected components:
      2  0  1  3  0  1  3
    ---------------------
  0:  0  1  0  1  2  0  0
  1:  0  0  2  1  0  5  0
  3:  0  0  0  4  0  0  1
  2:  1  0  0  1  0  0  0
  0:  3  1  0  0  1  1  0
  3:  0  0  0  3  0  0  1
  1:  0  0  1  0  0  3  4
upper triangular block matrix: 
  1  2  0  0  0  1  0
  1  1  0  1  3  0  0
  0  0  2  5  0  1  0
  0  0  1  3  0  0  4
  0  0  0  0  1  1  0
  0  0  0  0  0  4  1
  0  0  0  0  0  3  1
```
