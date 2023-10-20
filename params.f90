MODULE PARAMS
! This module contains parameters of calculations
! and variables 

USE CONSTANT

  IMPLICIT NONE

    INTEGER :: i,j;  
    INTEGER :: time_flush;
    ! save data Time_flush number of time steps

    DOUBLE PRECISION :: wo,w1, w2, m1, m2, k1, k2, a1, a2;
    ! coeficients of kinetic equations

    DOUBLE PRECISION :: v, x, dx, t, dt, dv,dv2, v_Ti, v_T, vmin,v_s;
    ! intrinsic variables 

    DOUBLE PRECISION :: ne, tqv, nb, Eb;
    ! plasma density, quasilinear time, 
    ! beam density, and beam energy density        

    DOUBLE PRECISION:: NuPlasma, omega;
	!Plasma frequency - defines plasma density
	! Omega = NuPlasma *2*Pi
    DOUBLE PRECISION:: Nbeam;
	! Electron beam density
	DOUBLE PRECISION:: Tmax;
	! Time of the calculations
	DOUBLE PRECISION:: Vbeam
	! Velocity of electron beam
	DOUBLE PRECISION:: Xmin, Xmax;
	! define the geometry of calculation area
	! Starting and End of simulation box
    DOUBLE PRECISION:: d;
	! Spatial size of the initial cloud
	DOUBLE PRECISION:: Time_save;
	! Time in seconds of each save
    
   DOUBLE PRECISION, DIMENSION(-Nv:Nv),public:: kx=0.,vx=0.;
        ! spectral energy density on kx

	DOUBLE PRECISION:: d1,d2,d3,d0, DF1,df2,df3,aa;

        DOUBLE PRECISION::spont_f,spont_w;

END MODULE PARAMS
