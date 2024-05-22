import math
import numpy as np
def efficiency(om):
  m=200 #mass of car 
  u1=0.001 #static friction coefficient
  rho=1.225 #air density
  a=1 #frontal area of car
  Cd=0.092 #coefficient of drag 
  r_in=0.18 #inner radius of wheel
  r_out=0.27 #outer radius of wheel
  omega =float(input(om))
  Ta=300
  visc=1.48*10**-5 #kinematic viscosity of air
  v_rotor = omega*r_in #circumferential speed of rotor
  RN = v_rotor*r_out/visc #reynolds number
  g=1.5*10**-3 #air gap spacing between stator and rotor
  pi=math.pi
  if RN > 0.8*10**5:
      #Regime III
      Cf=0.08/(((g/r_out)**0.167)*(RN**0.25))
      #Regime I
  else:
      Cf=2*pi*r_out/(g*RN)
  t=r_out*((m*9.81*u1)+(0.5*Cd*a*rho*(omega**2)*(r_out**2)))+0.5*Cf*rho*pi*(omega**2)*((r_out**5)-(r_in**5))
  #torque in steady state condition
  print("Torque : ",t)
  def find_winding_temp(Tw):
    B = 1.32 - 1.2*10**-3 *(Ta/2 + Tw/2 - 293) #magnetic remanence
    i = 0.561*B*t #RMS phase current
    R = 0.0575 * (1+0.0039*(Tw-293)) #resistance of windings
    Pc = 3 * i**2 * R #copper (ohmic) losses
    Pe = 9.602*10**-6 * (B*omega)**2 / R #eddy current losses
    Tw2 = 0.455 * (Pc + Pe) + Ta

    if np.absolute(Tw2-Tw) < 0.001:
      return Tw2
    else:
      return find_winding_temp(Tw2) #iterating till difference becomes less than 0.001
  Tw = find_winding_temp(Ta)
  B = 1.32 - 1.2*10**-3 *(Ta/2 + Tw/2 - 293)
  i = 0.561*B*t
  R = 0.0575 * (1+0.0039*(Tw-293))
  Pc = 3 * i**2 * R
  Pe = 9.602*10**-6 * (B*omega)**2 / R
  Pw=(omega**2)*(170.4*10**-6) #windage losses 
  t_f=0.5*Cf*rho*pi*(omega**2)*((r_out**5)-(r_in**5)) #frictional torque due to dynamic rolling resistance
  Pf=t_f*omega #wheel losses (frictional losses)
  P_out=t*omega #output power
  P_in=P_out+Pf+Pc+Pe+Pw #input power
  P_loss=Pf+Pc+Pe+Pw #total power loss
  eta=P_out*100/P_in #efficiency of motor
  # print("Input Power : ",P_in)
  # print("Copper Loss : ",Pc)
  # print("Eddy Current Loss : ",Pe)
  # print("Windage Loss : ",Pw)
  # print("Wheel Loss (Due to air drag) : ",Pf)
  # print("Total Power Loss : ",P_loss)
  # print("Output Power : ",P_out)
  # print("Efficiency : ",eta)
  return eta
