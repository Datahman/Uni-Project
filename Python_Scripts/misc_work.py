# DDA algorithm..
# Take two points, (X_1,Y_1), (X_2,Y_1) then calculate element differences, dx,dy.
import matplotlib.pyplot as plt 
def DDA_algorithm(x_1,y_1,x_2,y_2):
    dx = x_2 - x_1
    dy = y_2- y_1
    if abs(dx) > abs(dy):
        steps = abs(dx)
        Xincrement = dx / steps
    else:
        steps = abs(dy)
    
        Yincrement = dy / steps 
        print Yincrement
    #print steps
        x = 0
        y = 0
        x_,y_ = [], []
        for i in range(0,steps,++1):
            x = x + Xincrement
            y = y + Yincrement
            
            x_.append(x)
            y_.append(y)

        print y_

DDA_algorithm(20,40,100,80)
