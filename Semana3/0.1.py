import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from celluloid import Camera
from tqdm import tqdm


class Ball:
    
    # Init the object

    def __init__(self,x0,v0,a0,t,m,k,radius,Id):

        # Scalar variables

        self.dt = t[1] - t[0]
        self.m = m
        self.radius = radius
        self.Id = Id
        self.k = k

        # Vector variables

        self.x = x0
        self.v = v0
        self.a = a0
        self.f = self.m * self.a

        self.xVector = np.zeros((len(t),len(x0)))
        self.vVector = np.zeros((len(t),len(v0)))

        # Energy

        self.Ep = 0

        # Energy Vector

        self.EpV = np.zeros((len(t),1))   
        self.EkV = np.zeros((len(t),1))

    # Methods

    def Evolution(self,i):
        
        self.SetPositionVector(i,self.x)
        self.SetVelocityVector(i,self.v)
        self.SetKineticEnergyVector(i,self.GetKineticEnergy())
        self.SetPotentialEnergyVector(i,self.GetPotentialEnergy())

        self.a = self.f/self.m

        # Euler method

        self.x += self.dt * self.v
        self.v += self.dt * self.a
        
    def CheckWallLimits(self,limits,dim=2):
        
        for i in range(dim):
            if self.x[i] + self.radius > limits[i]:
                self.v[i] = - self.v[i]
            if self.x[i] - self.radius < - limits[i]:
                self.v[i] = - self.v[i]
        
    def Resetforce(self):
        self.f[:] = 0.
        self.a[:] = 0.
        self.Ep = 0.
        
    # Setters

    def SetPositionVector(self,i,x):
        self.xVector[i] = x
        
    def SetVelocityVector(self,i,v):
        self.vVector[i] = v

    def SetPotentialEnergyVector(self,i,u):
        self.EpV[i] = u

    def SetKineticEnergyVector(self,i,k):
        self.EkV[i] = k

    def Getforce(self, ball):

        d = np.linalg.norm(self.x - ball.GetPosition())
        compression = self.radius + ball.GetRadius() - d
        
        if compression > 0:
            F = self.k * compression**3 / d
            self.f = np.add(self.f,F*(self.x - ball.GetPosition()))
            self.Ep += self.k * compression**4 / 4

    # Getters

    def GetPosition(self):
        return self.x

    def GetPositionVector(self):
        return self.xVector

    def GetRxVector(self):
        return self.RxVector 
    
    def GetVelocityVector(self):
        print(self.vVector)
        return self.vVector

    def GetRvVector(self):
        return self.RvVector
    
    def GetAcelerationVector(self):
        return self.aVector
    
    def GetRadius(self):
        return self.radius

    def GetKineticEnergy(self):
        return 0.5*self.m*np.linalg.norm(self.v)**2

    def GetPotentialEnergy(self):
        return 0.5*self.Ep
    
    # Reduce Size

    def ReduceSize(self,factor):
        
        self.RxVector = np.array([self.xVector[0]])
        
        for i in range(1,len(self.xVector)):
            if i%factor == 0:
                self.RxVector = np.vstack([self.RxVector,self.xVector[i]])
        
        self.RvVector = np.array([self.vVector[0]])

        for i in range(1,len(self.vVector)):
            if i%factor == 0:
                self.RvVector = np.vstack([self.RvVector,self.vVector[i]])

# Discretización

t_max = 10
dt = 0.0001
t = np.arange(0,t_max+dt,dt)
Limits = [-20.,20.]

# Balls

ball_1 = Ball(np.array([-5.,-0.8]),np.array([20.,0.]),np.array([0.,0.]),t,1,100,2,2)
ball_2 = Ball(np.array([0.,-1.6]),np.array([0.,0.]),np.array([0.,0.]),t,1,100,2,1)
ball_3 = Ball(np.array([-15.,-15.]),np.array([0.,0.]),np.array([0.,0.]),t,1,100,2,3)

balls = [ball_1,ball_2,ball_3]

def RunSimulation(balls,t):
    
    for it in tqdm(range(len(t))):

        for i in range(len(balls)):
            for j in range(len(balls)):
                if i != j:
                    balls[i].Getforce(balls[j])
                    
            balls[i].Evolution(it)
            balls[i].Resetforce()
            balls[i].CheckWallLimits(Limits)
            
    return balls

runsimulation = RunSimulation(balls,t)

EnergyK = balls[0].EkV
EnergyP = balls[0].EpV
EnergyTotal = balls[0].EkV + balls[0].EpV

for i in range(1,len(balls)):
    EnergyK = np.add(EnergyK, balls[i].EkV)
    EnergyP = np.add(EnergyP, balls[i].EpV)
    EnergyTotal = np.add(EnergyTotal, balls[i].EkV + balls[i].EpV)

fig = plt.figure(figsize=(8,5))
ax1 = fig.add_subplot(1,1,1)
ax1.plot(t, EnergyTotal, label='Total Energy')
ax1.plot(t, EnergyK, label='Kinetic Energy')
ax1.plot(t, EnergyP, label='Potential Energy')
ax1.legend()
plt.show()
plt.savefig('Energy.png')

def ReduceTime(t,factor):

    for ball in balls:
        ball.ReduceSize(factor)

    Newt = []
    
    for i in range(len(t)):
        if i%factor == 0:
            Newt.append(t[i])

    return np.array(Newt)

redt = ReduceTime(t,100)

fig1 = plt.figure(figsize=(5,5))
ax = fig1.add_subplot(1,1,1)

def init():
    ax.set_xlim(-Limits[0],Limits[0])
    ax.set_ylim(-Limits[1],Limits[1])

def Update(i):

    plot = ax.clear()
    init()
    
    j = 0

    for ball in balls:
        x = ball.GetXxVector()[i,0]
        y = ball.GetXxVector()[i,1]

        vx = ball.GetXvVector()[i,0]
        vy = ball.GetXvVector()[i,1]

        circle = plt.Circle((x,y),ball.GetRadius)
        plot = ax.add_patch(circle)
        j += 1
    
    plot = ax.legend(loc=1)
    return plot

Animation = anim.FuncAnimation(fig1,Update,frames=len(redt),init_func=init)