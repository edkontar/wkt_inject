 Module constant

 ! This module contains constant of calculations
    
	DOUBLE PRECISION, PARAMETER::PI=3.141593

        DOUBLE PRECISION, PARAMETER::V_c=3e+10
        !speed of light    

	DOUBLE PRECISION, PARAMETER:: tmin=0;
	! Initial time

        DOUBLE PRECISION, PARAMETER:: vmax1=1.1e+9;
	! Maximal velocity

        DOUBLE PRECISION, PARAMETER:: vmax =20.e+9;
	! Maximal velocity

	DOUBLE PRECISION, PARAMETER:: e=4.8e-10;
        ! Electron charge

	DOUBLE PRECISION, PARAMETER::  m=9.10953447e-28;
	! electron mass

        DOUBLE PRECISION, PARAMETER::m_p=1.67264859e-24
	! proton mass

        DOUBLE PRECISION, PARAMETER:: T_e =1.5e+6

        DOUBLE PRECISION, PARAMETER:: T_i =1.5e+6
        ! Electron and ion temperatures

        DOUBLE PRECISION, PARAMETER:: k_b = 1.380662e-16
	!Boltzman constant

	INTEGER, PARAMETER:: nx = 160;
	! Number of x cells

	!INTEGER, PARAMETER:: nv1=400
	! Number of V cells

	INTEGER, PARAMETER:: nv =250*2
	! Number of V cells

        Character(*),parameter:: w_v_sequence ='vfw'

        Character(*),parameter:: e_sequence ='energy.dat'
       
        ! OUTPUT File names

        Character(*),parameter:: init_params = 'init.par'
        ! Input file name 

        CHARACTER (8):: xvfwn, st; 

End  Module constant