# Illustration of the OOP  class construction.
from visual.graph import *
# Functional programming example.

#graph = gdisplay(x=0,y=0, width=500,height=500,
#                    title='Motion of a satellite around a planet',
#                    xtitle='x',ytitle='y',xmax=5.0,xmin=-5.0,ymax=5.0,ymin=-5.0,
#                    foreground=color.black, background=color.white)
#                    
#moonfun = gcurve(color=color.red) # 2D graph, curve color in red.
#Radius = 4.0 # Planet orbit radius
#wplanet = 2.0 # Planet angular velocity.
#radius = 1.0 # Moon orbit radius around planet. 
#wmoon = 14.0 # Moon angular velocity around planet.
#for time in arange(0.,3.2,0.02):
#    rate(20)
#    x = Radius * math.cos(wplanet*time) + radius * math.cos(wmoon * time)
#    y = Radius * math.sin(wplanet*time) + radius * math.sin(wmoon * time)
#    moonfun.plot(pos=(x,y))
    
    
# Class example.

class OOPPlanet:
    
    def __init__(self,Radius,wplanet): # Initialize planet parameters in init constructor. 
        self.Radius = Radius
        self.wplanet = wplanet
      
        
    def getX(self,time): # These take one parameter, time !
        return self.Radius * math.cos(self.wplanet*time) 
    def getY(self,time):
        return self.Radius * math.sin(self.wplanet*time) 
        
    def scenario(self,mytitle,myxtitle,myytitle,xma,xmi,yma,ymi):
                
        graph = gdisplay(x=0,y=0, width=500,height=500,
                            title='Motion of a satellite around a planet',
                            xtitle='x',ytitle='y',xmax=5.0,xmin=-5.0,ymax=5.0,ymin=-5.0,
                            foreground=color.black, background=color.white)
    def position(self):
        planetfun = gcurve(color=color.red)
        for time in arange(0.,3,0.02):
            xp = self.getX(time) # Call get X method.
            yp = self.getY(time)
            planetfun.plot(pos=(xp,yp))

class OOPMoon(OOPPlanet):

    def __init__(self,Radius,wplanet,radius,wmoon): # Class constructor of Moon and planet paramaters.
        
        OOPPlanet.__init__(self,Radius,wplanet) # Abstract initialization as mother class doesn't contain such paramaters.
                                                # This step sens Radius, wplanet to OOPPlanet but without any instance made.

        self.radius = radius # initialized instances 
        self.wmoon = wmoon # initialized instances 
        #self.xm = xm
        #self.ym = ym
        
    def position(self):
        moonfun = gcurve(color=color.blue)
        for time in arange(0.,3,0.02):
            xm = self.getX(time) +  self.radius * math.cos(self.wmoon * time)  # Due to class inheritance,and abstraction here position() method overrides OOPPlanet method
            ym = self.getY(time) +  self.radius * math.sin(self.wmoon * time)           
            moonfun.plot(pos=(xm,ym))
            rate(10)
                
moon = OOPMoon(4.0,2.0,1.0,14.0) # Simulataneously initialize Moon and planet instances.
moon.scenario('Satellite orbir around planet','x','y','5','-5','5','-5') # Call scenario method on moon, defined in OOPCLass.
moon.position()
planet = OOPPlanet(4.0,2.0)
planet.position()


#Personal note: If multiple constructors present and one mis typed, then the next constructors would run! 
