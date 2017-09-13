# hypercontractivity

Estimator of the hypercontractivity coefficient from samples

This is a python package for estimating hypercontractivity coefficient from samples. The algorithm is proposed in Section 3 in http://arxiv.org/abs/1709.03999/

Sample Usage: 

    import HC_estimator

    s = HC_estimator.HC(x,y)
       
Input: 
  
    x -- 2D array of size n by d_x(1D array of size n if d_x = 1), where n is sample size and d_x is dimension of x

    y -- 2D array of size n by d_y(1D array of size n if d_y = 1), where n is sample size and d_y is dimension of y
       
Output: 

    A scalar s which is an estimate of s(X;Y)

Parameters: See comments in the code for more details

See demo.py for a simple example.


**************** Edited from here **************** 
See WHO_experiments.R for the R code for the WHO data experiment in Section 4.2. 

WHO dataset and MIC estimator:

    WHO.csv, MINE.jar, MINE.r  downloaded from http://www.exploredata.net/
**************** to here **************** 


Contact wgao9@illinois.edu
