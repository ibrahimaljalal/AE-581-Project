import numpy as np
import sympy as sp
from matplotlib import pyplot as plt
from sympy import *
import math
import time as time2


def bicycleController(wsString,alphaString,t0,tf,dt,x0=0,y0=0,phi0=0,idealCase=True):
    
    #the t variable will be treated as a symbolic variable
    t=sp.symbols("t")

    #the time is an array from t0 to tf incremented by dt ([t0 dt 2dt ...tf])
    time=np.arange(t0,tf+dt,dt)

    #the radius of the wheels in cm
    r=2
    #The distance between the centers of the front and back wheels in cm
    d=15

#####################################################################
    #The code between the #s will convert the input strings
    #into symbolic expresins than will subtitute the values 
    #of the time in the expresins and output an array of 
    #ws and alpha with the same elements as the time
    #len(time) =len(ws)=len(alpha)

    #Some parts of this sectin  are written to handle some possible errors 
    wsSymbolic=sp.sympify(wsString)
    alphaSymbolic=sp.sympify(alphaString)
    fws=sp.lambdify(t,wsSymbolic,modules=["numpy"])
    falpha=sp.lambdify(t,alphaSymbolic,modules=["numpy"])
    ws=fws(time)
    alpha=falpha(time)

    try:
        len(ws)
    except:
        ws=ws*np.ones(len(time))

    try:
        len(alpha)
    except:
        alpha=alpha*np.ones(len(time))


    
    if (idealCase==False):
        
        alpha=np.array(list(map(round,alpha)))

#####################################################################


#####################################################################
    #This section of the program is simply applying the
    #the mathematical models for the project example

    vs=r*ws
    v=vs*np.cos(np.radians(alpha))
    w=(vs/d)*np.sin(np.radians(alpha))

    phi=np.zeros(len(time))
    x=np.zeros(len(time))
    y=np.zeros(len(time))

    phi[0]=phi0
    x[0]=x0
    y[0]=y0


    for i in range(len(time)-1):
        x[i+1]=x[i]+v[i]*np.cos(phi[i])*dt
        y[i+1]=y[i]+v[i]*np.sin(phi[i])*dt
        phi[i+1]=phi[i]+w[i]*dt



    return (x,y,phi,time,ws,alpha)

#####################################################################



#####################################################################

#This function was removed from the report and presintation
#that is because the reasoning may not be completely right
#however the code works fine :) !!
def testForSharpEdge(x,y,dy2OverdxThreshhold,dx2OverdyThreshold):

    sharp=False

    dyOverdx=np.zeros(len(y))



    for i in range(len(y)-1):
        dyOverdx[i]=(y[i+1]-y[i])/(x[i+1]-x[i])

    dyOverdx[len(y)-1]=(dyOverdx[len(y)-2]-dyOverdx[len(y)-3])/(x[len(y)-2]-x[len(y)-3])*(x[len(y)-1]-x[len(y)-2])+dyOverdx[len(y)-2]


    
    dy2Overdx=np.zeros(len(y))

    for i in range(len(y)-2):

        dy2Overdx[i]=(dyOverdx[i+1]-dyOverdx[i])/(x[i+1]-x[i])

    
    dy2Overdx[len(y)-2]=(dy2Overdx[len(y)-3]-dy2Overdx[len(y)-4])/(x[len(y)-3]-x[len(y)-4])*(x[len(y)-2]-x[len(y)-3])+dy2Overdx[len(y)-3]
    dy2Overdx[len(y)-1]=(dy2Overdx[len(y)-2]-dy2Overdx[len(y)-3])/(x[len(y)-2]-x[len(y)-3])*(x[len(y)-1]-x[len(y)-2])+dy2Overdx[len(y)-2]
    
    
    
    
    dxOverdy=np.zeros(len(y))


    for i in range(len(y)-1):
        dxOverdy[i]=((x[i+1]-x[i]))/(y[i+1]-y[i])

    
    dxOverdy[len(y)-1]=(dxOverdy[len(y)-2]-dxOverdy[len(y)-3])/(y[len(y)-2]-y[len(y)-3])*(y[len(y)-1]-y[len(y)-2])+dxOverdy[len(y)-2]


    dx2Overdy=np.zeros(len(y))

    for i in range(len(y)-2):
        dx2Overdy[i]=(dxOverdy[i+1]-dxOverdy[i])/(y[i+1]-y[i])

    dx2Overdy[len(y)-2]=(dx2Overdy[len(y)-3]-dx2Overdy[len(y)-4])/(y[len(y)-3]-y[len(y)-4])*(y[len(y)-2]-y[len(y)-3])+dx2Overdy[len(y)-3]
    dx2Overdy[len(y)-1]=(dx2Overdy[len(y)-2]-dx2Overdy[len(y)-3])/(y[len(y)-2]-y[len(y)-3])*(y[len(y)-1]-y[len(y)-2])+dx2Overdy[len(y)-2]
    

    index=-1

    for i in range(len(y)):
        if (dy2Overdx[i]>dy2OverdxThreshhold and dx2Overdy[i]>dx2OverdyThreshold):
            sharp=True
            index=i
            break



    return dyOverdx,dy2Overdx,dxOverdy,dx2Overdy,sharp,index


#####################################################################


#####################################################################
#This part will give us the first and second figures

def plotMainResults(x,y,time,ws,alpha):
    fig=plt.figure()

    ax1=fig.add_subplot(221)
    ax1.set_title("ws vs t")
    ax1.set_xlabel("time in seconds")
    ax1.set_ylabel("angular speed in rads/s")

    ax2=fig.add_subplot(222)
    ax2.set_title("alpha vs t")
    ax2.set_xlabel("time in seconds")
    ax2.set_ylabel("steering angle in degrees")

    ax3=fig.add_subplot(212)
    ax3.set_title("y vs x")
    ax3.set_xlabel("x position in cm")
    ax3.set_ylabel("y position in cm")

    ax1.plot(time,ws,color=(1,0,0))
    ax2.plot(time,alpha,color=(0,1,0))
    ax3.plot(x,y,color=(0,0,1))

    plt.show()
#####################################################################


#####################################################################
#This function will simply plot the actual and ideal cases in one plot.
#Note that the limits for x and y are from 0 to 100 cm
def plotIdealAndActual(xi,yi,xa,ya):
    plt.xlim((0,100))
    plt.ylim((0,100))
    plt.xlabel("x position in cm")
    plt.ylabel("y position in cm")
    plt.plot(xi,yi,label="Ideal Case")
    plt.plot(xa,ya,label="Real Case")
    plt.legend()
    plt.show()
#####################################################################



#####################################################################

#The actual program starts from here and will use the previous functions
if __name__=="__main__":

    
    print("""
    This program will show three figures according to the user input

    1-First figure is for the ideal case
    2-Second figure is for the actual case (the angle is a multiple of 1 degree)
    3-Third is for ideal vs real

    Notes:
    1-Write ** insted of ^ if you want to use power in your expression
    2-The expresions are functions of t 
    3-Some expresins my require to decrease dt to get clear graphs
    4-The last figure limits are from 0 100 cm for both x and y 
    so make sure that you put your first two expresins such that
    the output x and y will be within this region  
    
    """)

    #This part is for the AE581.exe file
    # wsString=str(input("ws(t) in rads/s = "))
    # alphaString=str(input("alpha(t) in degrees = "))
    # t0=float(input("t0 in seconds = "))
    # tf=float(input("tf in seconds = "))
    # dt=float(input("dt in seconds = "))
    # x0=float(input("x0 in cm = "))
    # y0=float(input("y0 in cm = "))
    # phi0=float(input("phi0 in radians = "))
    # sleep=2
    
    
    #This part is for the AE581.py file
    wsString="5**3"
    alphaString="sin(t)"
    t0=0
    tf=7
    dt=0.0001
    x0=0
    y0=0
    phi0=0
    sleep=0.1


    (xi,yi,phii,time,wsi,alphai)=bicycleController(wsString,alphaString,t0,tf,dt,x0=0,y0=0,phi0=0,idealCase=True)
    (xa,ya,phia,time,wsa,alphaa)=bicycleController(wsString,alphaString,t0,tf,dt,x0=0,y0=0,phi0=0,idealCase=False)

    plotMainResults(xi,yi,time,wsi,alphai)
    plotMainResults(xa,ya,time,wsa,alphaa)
    plotIdealAndActual(xi,yi,xa,ya)

    differenceInX=xi[len(xi)-1]-xa[len(xa)-1]
    differenceInY=yi[len(yi)-1]-ya[len(ya)-1]

    distanceDifference=math.sqrt(differenceInX**2+differenceInY**2)

    print("\n")
    print("The difference between the ideal and actual final positions is "+str(distanceDifference)+" cm")
    print("\n")
    time2.sleep(sleep)





    






























