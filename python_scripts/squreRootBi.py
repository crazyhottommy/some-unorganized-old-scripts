def squreRootBi(x, epsilon):
    '''assume x>=0 and epsilon > 0
    Return y s.t y*y is within epsilon of x'''
    assert x>= 0, 'x must be non-negative, not' + str(x) 
    assert epsilon > 0 , 'epsilon must be positive, not' + str(epsilon)
    low = 0
    high= max (x, 1.0) # incase x is smaller than 1, 0.5**2=0.25
    guess=(low+high)/2.0
    count=1
    while abs(guess **2 -x ) > epsilon and count <= 100:
        print 'low:', low, 'high:',high, 'guess:', guess
        if guess **2 < x:
            low=guess
        else:
            high=guess
        guess=(low+high)/2.0
        count+=1
    assert count <= 100, 'iteration count exceeded'
    print 'Bi method, num.interations:', count, 'Estimate:', guess
    return guess
