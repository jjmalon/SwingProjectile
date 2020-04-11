import numpy as np
import matplotlib.pyplot as plt
import math
from statistics import mean

class swingProjectile:

    def __init__(self,x,l):
        self._g = 9.81
        self._x = x
        self._l = l
        self._d = 0
        self._phi =0
        self._theta =0
        self._firstderiv = 0

    def __repr__(self):
        return 'Phi : {1}, Theta: {2}, Distance: {0}, First Deriv" {3}'.format('{0:.02f}'.format(self._d), '{0:.02f}' .format(np.degrees(self._phi)), '{0:.02f}'.format( np.degrees(self._theta)), '{0:.02f}'.format(self._firstderiv) )

    def getDistance(self):
        return self._d

    def getPhi(self):
        return self._phi

    def getTheta(self):
        return self._theta

    def getFirstDeriv(self):
        return self._firstderiv

    def getLaunchVelocity(self,phi,theta):
        return (self._g*self._l*(math.cos(theta)- math.cos(phi)))**(0.5)
    
    

    def calcDistance(self,phi,theta,h):

        self._phi = phi
        self._theta = theta
        v = self.getLaunchVelocity(phi,theta)
        ky = (self._l - self._l*math.cos(theta) + self._x)
        kx = (self._l*math.sin(theta))

        fdlv = self.getLaunchVelocity(phi,theta + h)
        bdlv = self.getLaunchVelocity(phi,theta - h)
        fkx = (self._l*math.sin(theta + h))
        bkx = (self._l*math.sin(theta - h))
        fky = (self._l - self._l*math.cos(theta + h) + self._x)
        bky = (self._l - self._l*math.cos(theta - h) + self._x)

        fwddiff = fkx + fdlv*((fdlv*math.sin(theta + h) + ((fdlv**2)*(math.sin(theta+h))**2 + 2*self._g*fky)**0.5)/self._g)*math.cos(theta+h)
        backdiff = bkx + bdlv*((bdlv*math.sin(theta - h) + ((bdlv**2)*(math.sin(theta - h))**2 + 2*self._g*bky)**0.5)/self._g)*math.cos(theta - h)
        firstderiv =   (fwddiff - backdiff)/(2*h) 
        
        self._d = kx + v*((v*math.sin(theta) + ((v**2)*(math.sin(theta))**2 + 2*self._g*ky)**0.5)/self._g)*math.cos(theta)
        self._firstderiv = firstderiv

    def generateData():
        
        swinglengths = np.linspace(2.4,2.7,4)
        angles = np.linspace(math.pi/2,np.radians(10),5)
        x = 0.35
        h = 0.01
        swinglengthdict = {}
        for swinglength in swinglengths:

            swingprojectiles = []
            swinglengthdict[swinglength] = swingprojectiles
   
            for angle in angles:

                phi = angle
                thetas = np.linspace(0,phi,20)

                for theta in thetas:

                    swingProjectileobj = swingProjectile(x,swinglength)                   
                    swingProjectileobj.calcDistance(phi,theta,h)                    
                    swingprojectiles.append(swingProjectileobj)

        return swinglengthdict,angles

    def testExternalModel(swinglength):
        x = 0.35
        h = 0.01
        theta = np.radians(60)
        angles = np.linspace(math.pi/2,np.radians(60),30)
        
        a = swinglength * 0.867
        b = 1.75
        model_a = swinglength*math.sin(theta)
        projectileDistance = a + b
        
        print(a,model_a,projectileDistance)
        
        
        swingprojectiles = []
        
        
        for phi in angles:
            
            swingProjectileobj =  swingProjectile(x,swinglength) 
            swingProjectileobj.calcDistance(phi,theta,h)
            swingprojectiles.append(swingProjectileobj)
            
        print(swingprojectiles)  
            
        

    def getMaxItem(swingdatalst):

        maxitem = None

        for index in range(0,len(swingdatalst)-1):
                item = swingdatalst[index]
                itemfwd = swingdatalst[index +1]

                if (item.getFirstDeriv() > 0 and itemfwd.getFirstDeriv() <0):
                         maxitem = item
                         break

        if maxitem ==None:
            maxitem = swingdatalst[-1]

        return maxitem

    def linearRegress(xs,ys):
          m = (((mean(xs)*mean(ys)) - mean(xs*ys)) / ((mean(xs)*mean(xs)) - mean(xs*xs)))
    
          b = mean(ys) - m*mean(xs)
    
          return m, b

    def createMaxDistanceProjectilePlot(swingdatadict,angles):


        maxitems = list()
        lengths = sorted(swingdatadict.keys())
        xplot = lengths

        for length in lengths:
            
            swingprojectiles = swingdatadict[length]
            maxphi = angles[0]
                              
            swingdatafltd = [x for x  in swingprojectiles if maxphi == x.getPhi()]
            
            maxitem = None                    
            maxitem = swingProjectile.getMaxItem(swingdatafltd)
            maxitems.append(maxitem)

        
        yplot = [ x.getDistance() for x in maxitems ]

        m,b = swingProjectile.linearRegress(np.array(xplot),np.array(yplot))

        labelstr = r'$y$ = {0}$x$ + {1}'.format('{0:.02f}'.format(m),'{0:.02f}'.format(b))
        
        plt.plot(xplot,yplot,marker='o',linestyle='dashed',label=labelstr)
        plt.legend(loc="lower right")
        plt.xlabel(r'Swing Lengths ($m$)')
        plt.ylabel(r'Projectile Distance ($m$) ')
        plt.grid(axis='both')
        plt.title(r'Swing length ($m$) vs Projectile Distance ($m$) ')

        plt.show()
                       
        
    def createProjectilePlot(swingdatadict,angles,samplelength):

        for length,swingprojectiles in swingdatadict.items():

            if length == samplelength:

                for phi in angles:
                    swingdatafltd = [x for x  in swingprojectiles if phi == x.getPhi()]
                    xplot = [ np.degrees( x.getTheta()) for x in swingdatafltd ]
                    yplot = [ x.getDistance() for x in swingdatafltd ]

                    maxitem = None                    
                    maxitem = swingProjectile.getMaxItem(swingdatafltd)
                    
                    print(maxitem)

                    if maxitem is not None:
                          labelstr =  r'$\phi$ : {0}, Max : {1}'.format('{0:.02f}'.format(np.degrees(phi)), '{0:.02f}'.format(maxitem.getDistance()))
                    else:
                          labelstr=r'$\phi$ : {0}, '.format( '{0:.02f}'.format(  np.degrees(phi)))
                    
                    plt.plot(xplot,yplot,label=labelstr)
                    plt.legend(loc="lower right")
                    plt.xlabel(r'Launch Angle ($\theta$)')
                    plt.ylabel('Projectile Distance (M) ')
                    plt.grid(axis='both')
                    plt.title(r'Launch Angle ($\theta$) vs Projectile Distance ($m$) ')


        plt.show()
    

samplelength = 2.7
swingdatadict,angles = swingProjectile.generateData()

def plot1():
    swingProjectile.createProjectilePlot(swingdatadict,angles,samplelength)

def plot2():
    swingProjectile.createMaxDistanceProjectilePlot(swingdatadict,angles)

def evaluateModel():
    swingProjectile.testExternalModel(samplelength)


evaluateModel()




        
        


