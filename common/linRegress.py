#150507 Mary K. Thompson
#got this from here: Greg Pinero, http://www.answermysearches.com/how-to-do-a-simple-linear-regression-in-python/124/
#modified such that it will return the individual residuals
from math import sqrt

def linreg(X, Y):
    """
    Summary
        Linear regression of y = ax + b
    Usage
        real, real, real = linreg(list, list)
    Returns coefficients to the regression line "y=ax+b" from x[] and y[], and R^2 Value
    """
    if len(X) != len(Y):  raise ValueError, 'unequal length'
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in map(None, X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
    meanerror = residual = 0.0
    resis=[]
    for x, y in map(None, X, Y):
        meanerror = meanerror + (y - Sy/N)**2
        residual = residual + (y - a * x - b)**2 #this is the sum of the squared resdiuals, but I want it to just return them
        resis.append(y-a*x-b)
    RR = 1 - residual/meanerror
    ss = residual / (N-2)
    Var_a, Var_b = ss * N / det, ss * Sxx / det
   
    return a, b, RR, resis

if __name__=='__main__':
    #testing
    X=[1,2,3,4]
    Y=[357.14,53.57,48.78,10.48]
    
    a, b, RR, res= linreg(X,Y)
    print 'a', a #slope
    print 'b', b #y-intercept
    print 'R2', RR #R-squared value
    print 'residuals', res #residuals
    
    '''
    #should be:
    a -104.477
    b 378.685
    R2 0.702499063713
    residuals [82.93199999999996, -116.16099999999994, -16.473999999999933, 49.703000000000145]
    '''