# wave-function
This repository contains three different programs.  
* box_1st_energy  
This program solves Schroedinger equation under the condition of box potential whose magnitude is infinite.  
Particles cannot trespass into the potential walls, which means in the very edges of zero potential field, the wave function's boundary condition is set to zero (Neumann boundary condition).
* box_Nth_energy  
The program above solves Schroedinger equation, but for the 1st energy eigenvalue.  
In contrast to the previous program, this program can solve energy eigenvalue of any levels.
* harmonic_energy  
This program solves Schroedinger equation under the condition of potential 1/2 kx^2 which is often called as "Harmonic potential (well)".  
The boundary condition is the same as previous programs which is \phi(-a)=\phi(a)=0
![Wave function under harmonic potential (n=3)](https://raw.githubusercontent.com/StellaContrail/wave-function/master/harmonic_energy/harmonic_3rd.png)  
Wave function under harmonic potential (n=3)
