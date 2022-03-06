import numpy as np

class Particle:
    
    def __init__(self, r0,v0,a0,t,m,radius,Id):
        
        # Variables

        self.dt = t[1]-t[0]
        self.r = r0
        self.v = v0
        self.a = a0
        self.m = m
        self.radius = radius
        self.Id = Id
        
        # Vectors

        self.rVector = np.zeros((len(t),len(r0)))
        self.vVector = np.zeros((len(t),len(v0)))
        self.aVector = np.zeros((len(t),len(a0)))
        self.L = np.zeros(len(r0))
        
        # Momentum

        self.MomentumVector = np.zeros((len(t),len(v0)))
        
        # Energias

        self.EpVector = np.zeros((len(t),1))
        self.EkVector = np.zeros((len(t),1))
        
        # Otther variables
        
        self.Ep = 0.
        self.Force = self.m * self.a
        self.G = 4*np.pi**2
        self.rp = r0
        self.vp = v0
        
    def Evolution(self,i):
        
        # Set vectors
        
        self.SetPosition(i,self.r)
        self.SetVelocity(i,self.v)
        self.SetMomentum(i,self.m*self.v)
        
        # Energy 
        
        self.SetEk(i,self.GetKineticEnergy())
        self.SetEp(i,self.GetPotentialEnergy())
        
        self.a = self.Force/self.m
        
        if i == 0:
            self.rp = self.r
            self.r = self.rp + self.dt * self.v
        else:
            self.rf = 2*self.r - self.rp + self.a * self.dt**2            
            self.v = ( self.rf - self.rp ) / (2*self.dt)
            self.rp = self.r
            self.r = self.rf
            
    def ResetForce(self):
        
        self.Force[:] = 0.
        self.a[:] = 0.
        self.Ep = 0.
        
    # Setters
    
    def SetPosition(self,i,r):
        self.rVector[i] = r
        
    def SetVelocity(self,i,v):
        self.vVector[i] = v   
        
    def SetMomentum(self,i,p):
        self.MomentumVector[i] = p
        
    def SetAngularMomentum(self,i,r,p):
        self.L[0] = r[1]*p[2] - r[2]*p[1]
        self.L[1] = -(r[0]*p[2] - r[2]*p[0])
        self.L[2] = r[0]*p[1] - r[1]*p[0]
            
    def SetEk(self,i,Ek):
        self.EkVector[i] = Ek
    
    def SetEp(self,i,Ep):
        self.EpVector[i] = Ep

    # Getters
    
    def GetForce(self,p):
        eps=0.1
        d = np.linalg.norm(self.r-p.GetPosition())
        Fn = -self.G*self.m*p.m/(d**2 + 0.1**2)**3
        self.Force = np.add(self.Force,Fn*(self.r-p.GetPosition()))
        self.Ep += - self.G * self.m * p.m / d
    
    def GetPosition(self):
        return self.r
    
    def GetPositionVector(self):
        return self.rVector
    
    def GetReducePosition(self):
        return self.RrVector
        
    def GetVelocityVector(self):
        return self.vVector    
    
    def GetMomentumVector(self):
        return self.MomentumVector
        
    def GetReduceVelocity(self):
        return self.RvVector
     
    def GetKineticEnergy(self):
        return 0.5*self.m*np.linalg.norm(self.v)**2
    
    def GetPotentialEnergy(self):
        return 0.5*self.Ep 
   
    def GetNetForce(self):
        return self.Force

    def GetR(self):
        return self.radius

    def ReduceSize(self,factor):
        self.RrVector = np.array([self.rVector[0]])
        for i in range(1,len(self.rVector)):
            if i%factor == 0:
                self.RrVector = np.vstack([self.RrVector,self.rVector[i]])      
        self.RvVector = np.array([self.vVector[0]])
        for i in range(1,len(self.vVector)):
            if i%factor == 0:
                self.RvVector = np.vstack([self.RvVector,self.vVector[i]])