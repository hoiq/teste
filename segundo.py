omega=(0,0,4)
m=1
k=2000
com2=omega(0.0333333, 0.0166667, 0.0333333)
trailpoints=50
massA=(pos=vector(0,0,0)-com2, radius=0.01, color=color.yellow, make_trail=True, retain=trailpoints)
massA.m=0.5
massA.v=(omega,massA.pos)
massA.p=massA.m*massA.v
massB=(pos=omega(0.1,0,0)-com2, radius=0.01, color=color.red, make_trail=True,retain=trailpoints)
massB.m=m
massB.v=(omega,massB.pos)
massB.p=massB.m*massB.v
massC=(pos=omega(0,.1,0)-com2, radius=0.01,color=color.cyan, make_trail=True,retain=trailpoints)
massC.m=m/2
massC.v=(omega,massC.pos)
massC.p=massC.m*massC.v
temprD=-massA.m*massA.pos-massB.m*massB.pos-massC.m*massC.pos
massD=(pos=temprD, radius=0.01, make_trail=True,retain=trailpoints)
massD.m=m
#massD.p=massD.m*vector(0,0,0)
massD.v=(omega,massD.pos)
massD.p=massD.m*massD.v
rAB=massB.pos-massA.pos
rBC=massC.pos-massB.pos
rCA=massA.pos-massC.pos

rDB=massB.pos-massD.pos
rCD=massD.pos-massC.pos
rDA=massA.pos-massD.pos

com=(massA.pos*massA.m+massB.m*massB.pos+massC.m*massC.pos+massD.m*massD.pos)/(massA.m+massB.m+massC.m+massD.m)


springAB=cylinder(pos=massA.pos, axis=rAB, radius=0.001)
springBC=cylinder(pos=massB.pos, axis=rBC, radius=0.001)
springCA=cylinder(pos=massC.pos, axis=rCA, radius=0.001)

springDB=cylinder(pos=massD.pos, axis=rDB, radius=0.001)
springCD=cylinder(pos=massC.pos, axis=rCD, radius=0.001)
springDA=cylinder(pos=massD.pos, axis=rDA, radius=0.001)

LAB=mag(rAB)
LBC=mag(rBC)
LCA=mag(rCA)
LDB=mag(rDB)
LCD=mag(rCD)
LDA=mag(rDA)


print("mass A = ",massA.pos," m")
print("mass B = ",massB.pos," m")
print("mass C = ",massC.pos," m")
print("mass D = ",massD.pos," m")

def I(ro,A,B,C,D):
  #this function takes 4 mass objects and 
  #returns the moment of inertia tensor with respect to point ro
  ms=[A,B,C,D]
  Ixx=0
  Ixy=0
  Ixz=0
  Iyy=0
  Iyz=0
  Izz=0
  for pp in ms:
    Ixx=Ixx+pp.m*((pp.y-ro.y)**2+(pp.z-ro.z)**2)
    Ixy=Ixy-pp.m*(pp.x-ro.x)*(pp.y-ro.y)
    Ixz=Ixz-pp.m*(pp.x-ro.x)*(pp.z-ro.z)
    Ix=vector(Ixx,Ixy,Ixz)
    Iyy=Iyy+pp.m*((pp.x-ro.x)**2+(pp.z-ro.z)**2)
    Iyz=Iyz-pp.m*(pp.y-ro.y)*(pp.z-ro.z)
    Iy=vector(Ixy,Iyy,Iyz)
    Izz=Izz+pp.m*((pp.x-ro.x)**2+(pp.y-ro.y)**2)
    Iz=vector(Ixz,Iyz,Izz)
  return([Ix,Iy,Iz])

print(I(vector(0,0,0),massA,massB,massC,massD))
  