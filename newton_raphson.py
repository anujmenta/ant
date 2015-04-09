def newton_raphson(x_, f, df, epsilon=0.001):
    n = 1
    while True:
        x = x_ - f(x_)/df(x_)
        if abs(x - x_) < epsilon:
            return x
        x_ = x
        print "Value of x after %s iterations is: %s" % (n, x_)
        n = n+1
        if n > 50:
            print "Divergent too many iterations"
