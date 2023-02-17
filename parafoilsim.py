import numpy as np
import matplotlib.pyplot as plt
import math as meth

t = []
vX = []
vY = []
headingAngle = []

class State:
    def __init__(self):
        self.pos = np.array((0,0,0),dtype='float64')
        self.vel = np.array((0,0,0),dtype='float64')
        self.dThetaDs = 1 + 0j
        self.heading = 1 + 0j
        
        self.time = 0
        self.dt = 0.01
        self.glideRatio = 6.5
        self.zVelLast = 0

    def update(self):
        self.time += self.dt
        dz = self.zVelLast * self.dt
        dx = dz * self.glideRatio
        self.heading *= self.dThetaDs
        self.vel[0] = np.real(self.heading) * dx
        self.vel[1] = np.imag(self.heading) * dx
        self.pos += self.vel * self.dt
        self.zVelLast = self.vel[2]

        t.append(self.time)
        vX.append(self.vel[0])
        vY.append(self.vel[1])
        headingAngle.append(np.angle(self.heading))

    def airDensityFromAltitude(self,z):  #https://en.wikipedia.org/wiki/Density_of_air REDO under troposphere
         p = 0 #air density [kg/m^3]
         StandAtmPres = 101325 #[Pa] at sea lvel
         StandTemp = 288.12 #[K] at sea level
         gravity = 9.80665 #[mn/s^2]
         tempLapseRate = 0.0065 #[K/m]
         gasConstant =  8.31432 #[N*m / mol*k]
         molarMass = 0.0289644 #molar mass of dry air [kg/mol]
         intermediate = (StandAtmPres*molarMass)/(gasConstant*StandTemp)
         intermediate2 = (1-tempLapseRate*z/StandTemp)**(((gravity*molarMass)/(gasConstant*tempLapseRate))-1)
         
         return (intermediate *intermediate2)
          
    def convertAirDensitytoImperial(self, airDensity): #converts airDensity to imperial units in slugs/ft^3 from kg/m^3
        convert = airDensity  #storing airDensity in new variable called convert so the imperial value is not over written when given a new altitude reading
        convert *= 0.00194032 #converting airDensity from metric to imperial unit by conversion factor 0.00194032
        return (convert) 
        
    def rateOfDescent(self,height): #Parachutes 101 in CRT drive
        # air density
        airDensity = self.convertAirDensitytoImperial(self.airDensityFromAltitude(height))

        WeightofLoad= 8.8 #load + parachute, (lbs) 8.8lbs (all according to Luke)
        SurfaceArea =1.9375 #canopy surface area (ft^2) 
        dragCoefficient = 1.6 #related to SA (unitless) ~1.6
        densityO = 0.237689 #Density at sea level
        calculation = WeightofLoad *2/(SurfaceArea*dragCoefficient*densityO)
        calculation = meth.sqrt(calculation)
        intermediate = 1/(meth.sqrt(airDensity/densityO))
        return (calculation*intermediate)/3.281
        # returns in metric again


rocket = State()
rocket.pos = np.array((0,0,15000),dtype='float64')
rocket.vel = np.array((0,0,0),dtype='float64')
rocket.acc =np.array((0,0,0),dtype='float64')

x = []
y = []
z = []

rocket.dThetaDs = np.exp(1j * 0.01)

while rocket.pos[2] > 0:
    x.append(rocket.pos[0])
    y.append(rocket.pos[1])
    z.append(rocket.pos[2])
        
    #fallRate = -rocket.rateOfDescent(rocket.pos[2])
    #rocket.vel[2] = fallRate
    rocket.vel[2] = -1000
    
    rocket.update()
    
pos = plt.figure().add_subplot(projection='3d')
pos.plot(x, y, z, label='parametric curve')

""""
vel = plt.figure().add_subplot(projection='3d')
vel.plot(vX, vY, t, label='parametric curve')

pos = plt.figure().add_subplot()
pos.plot(t, headingAngle, label='parametric curve')

plt.figure().add_subplot().plot(t, vX, label='parametric curve')
"""

plt.show()