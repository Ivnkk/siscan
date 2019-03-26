function L = Laplace(I,x,y)

L = I(x+1,y) + I(x-1,y) + I(x,y-1) +I(x,y+1) - 4*I(x,y);