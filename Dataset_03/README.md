**Example 3**  
The dataset can be found [here](https://drive.google.com/file/d/1JrNm7KxWd2Mdmphj9nDOJUsZhrDaPGzh/view?usp=sharing)  
PANCAN with no missing values  
10884 samples  
20531 variables
  

| Comp | Iter | Eig | Var(%) | VarCum(%) | ConvValue| Orthogonality |
| --- | ---- | --- | ----- | ------ | ----- | ----- |
|   1  |  72 |  1349.9 | 12.40 | 12.40 | 1.5e-24 |       1 |
|   2  | 252 |   799.0 |  7.34 | 19.74 |   4e-24 | 2.6e-08 |
|   3  |  89 |   695.0 |  6.39 | 26.13 | 2.8e-24 | 2.3e-08 |
|   4 |  402 |   476.6 |  4.38 | 30.51 | 2.4e-24 | 1.9e-08 |
Elapsed time is 146 seconds   4.22 GB  

  
  
```matlab
tic,[coeff] = pca(X,'Algorithm','svd','NumComponents',4);toc
```
Elapsed time is 479.761866 seconds. 12.73 GB

```matlab
tic,[coeff] = pca(X,'Algorithm','eig','NumComponents',4);toc
```
Elapsed time is 638.926541 seconds. 18.85 GB

```matlab
tic,[coeff] = pca(X,'Algorithm','als','NumComponents',4);toc
```
Elapsed time is 1073.187737 seconds. 12.73 GB