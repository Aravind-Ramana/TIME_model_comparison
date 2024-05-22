import plotly.graph_objects as go
import numpy as np
import math

omega_range = np.linspace(20.0, 100.0, 1000)
tau_range = np.linspace(0.0, 80.0, 1000)
m=200
u1=0.001
rho=1.225
a=1
Cd=0.11
r_in=0.18
r_out=0.27
Ta=float(input("Enter Ambient Temperature in K : "))
visc=1.48*10**-5
g=1.5*10**-3
pi=math.pi

temps = []
Tm = []
B = []
R = []

for j in range(1000):
    omega = omega_range[j]  
    v_rotor = omega * r_in  
    
    RN = v_rotor * r_out / visc  
    
    if RN > 0.8 * 10**5:
        Cf = 0.08 / (((g / r_out)**0.167) * (RN**0.25))  
    else:
        Cf = (2 * pi * r_out) / (g * RN) 
    
    t = (r_out * ((m * 9.81 * u1) + (0.5 * Cd * a * rho * (omega**2) * (r_out**2))) +
         0.5 * Cf * rho * pi * (omega**2) * ((r_out**5) - (r_in**5)))
    def find_winding_temp(Tw):
      B = 1.32 - 1.2*10**-3 *(Ta/2 + Tw/2 - 293)
      i = 0.561*B*t
      R = 0.0575 * (1+0.0039*(Tw-293))
      Pc = 3 * i**2 * R
      Pe = 9.602*10**-6 * (B*omega)**2 / R
      Tw2 = 0.455 * (Pc + Pe) + Ta

      if np.absolute(Tw2-Tw) < 0.001:
        return Tw2
      else:
        return find_winding_temp(Tw2)
    Tw = find_winding_temp(Ta)
    temps.append(Tw)
    Tm.append(Ta/2 + temps[j]/2)
    B.append(1.32 - 1.2 * 0.001 * (Tm[j] - 293))
    R.append(0.0575*(1+0.0039*(temps[j]-293)))

[X, Y] = np.meshgrid(omega_range, tau_range)
motor_loss = (170.4 * 0.000001 * X**2 + ((9.602 * 0.000001 * (B * X)**2)/R) + (3 * (R * (0.561 * Y * B)**2))+
              (0.5 * Cf * rho * pi * (X**3) * ((r_out**5) - (r_in**5))))
n = 100*X*Y / (X*Y + motor_loss ) 
fig = go.Figure(data = go.Contour(x = omega_range, y = tau_range, z = n, 
                colorscale='Rainbow',
                contours=dict(
                start=75,
                end=100,
                size=0.5,
                ),
                ))
fig.update_layout(
    title="Efficiency vs torque vs angular speed of motor",
    xaxis_title="Angular speed (rad/sec)",
    yaxis_title="Torque (Nm)",
    legend_title="Efficiency (%)", )
fig.show()