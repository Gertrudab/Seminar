# Seminar
Emission Imaging

The goal of this assignment is implement code supporting Jeffrey A. Fessler's book's "Image Reconstruction: Algorithms and Analysis" 8th chapter: Emission Imaging. 

Implementation contains 4 fully implemented functions:
```console
P_y (µ_n, µ_t, t1, t2, y, n_d, r)
```
This function calculates the joint distribution of the recorded events Y

```console
lambda_x (x, µ_n, µ_t, t1, t2)
```
This function calculates the average emission rate density function over the time interval t_2-t_1
```console
calculate_y (i, µ_n, µ_t, t1, t2, r)
```
This function calculates the mean of Y_i (number of events recorded by the i_th detector)
```console
calculate_y (i, µ_n, µ_t, t1, t2, r)
```
This function calculates the log-likelihood for the average emission rate density function in binned-mode