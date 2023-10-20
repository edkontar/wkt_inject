MODULE INIT_BEAM

USE CONSTANT, ONLY: Pi;

USE PARAMS, ONLY: Nbeam,vbeam,d;

 !MODULE contains initial distribution functions of
 !electron beam: 
 !1) FUNCTION F0 - monoenergetic electron beam
 !2) FUNCTION F1 - linealy growing with velocity, v< Vbeam
 !To change the using function from 1 to 2 =>
 ! replace F1 -> F0 and F0 -> F1
  
   IMPLICIT NONE

   CONTAINS
 DOUBLE PRECISION FUNCTION F0(vel)
 ! Maxwellian electron distribution function of the plasma ONLY
 
   USE PARAMS
   IMPLICIT NONE

   DOUBLE PRECISION, INTENT(IN) :: vel; 
   DOUBLE PRECISION :: a, delta_v, delta;
   	    
      delta_v = 0.3*Vbeam;
      delta= 4.
          	     
      a=Nbeam/(v_T*(delta-1));
		   
      if (vel > 0 ) then 
      F0=ne*exp(-vel*vel/(2.*v_t*v_t))/(sqrt(2.*pi)*v_t)
	   ! only maxwellian and NO beam
	   ! beam injected later

      else
      F0=ne*exp(-vel*vel/(2.*v_t*v_t))/(sqrt(2.*pi)*v_t)
      end if

      ! function returns background plasma only 
END FUNCTION F0


 DOUBLE PRECISION FUNCTION F0_beam(vel)

 ! electron distribution function of the beam and plasma 
   USE PARAMS
   IMPLICIT NONE

   DOUBLE PRECISION, INTENT(IN) :: vel; 
   DOUBLE PRECISION :: a, delta_v, delta;
   	    
      delta_v = 0.3*Vbeam;
      delta= 4.
          	     
      a=Nbeam/(v_T*(delta-1));
		   
      if (vel > 0 ) then 
      F0_beam=ne*exp(-vel*vel/(2.*v_t*v_t))/(sqrt(2.*pi)*v_t) &
      +a*((30.*V_T)**2/((30.*v_T)**2+Vel**2))**delta
	  ! maxwellian and beam

      else
      F0_beam=ne*exp(-vel*vel/(2.*v_t*v_t))/(sqrt(2.*pi)*v_t)
      end if

      
END FUNCTION F0_beam

DOUBLE PRECISION FUNCTION F2(vel)

 ! Monoenergetic electron distribution function of the beam 
   USE PARAMS
   IMPLICIT NONE

   DOUBLE PRECISION, INTENT(IN) :: vel; 
   DOUBLE PRECISION :: a, delta_v;
   	    
      delta_v = 0.3*Vbeam;
          	     
      a=Nbeam/(sqrt(pi)*delta_v);
		   
      !if (abs(vel)<=2.*vbeam) then
      F2=a*exp(-(vel-vbeam)**2/delta_v**2)&
       +ne*exp(-vel*vel/(v_t*v_t))/(sqrt(pi)*v_t)
      !else
      !F0=0.
      !end if
      		  
      
END FUNCTION F2

DOUBLE PRECISION FUNCTION F1(v1, x1)

!USE CONSTANT
!USE PARAMS
!, ONLY: Nbeam,vbeam,d;

 ! Linear electron distribution function of the beam 
 ! v < Vbeam
 
   IMPLICIT NONE

   DOUBLE PRECISION, INTENT(IN) :: v1, x1;
   DOUBLE PRECISION :: a
   
      a=2.*Nbeam/Vbeam
		   
      if ((v1<=Vbeam).and. (v1>0.)) F1=a*exp(-(x1/d)*(x1/d))*v1/Vbeam
      if ((v1>Vbeam) .or. (v1 .LE. 0.)) F1 = 0.
	  
		  
END FUNCTION F1

END MODULE INIT_BEAM
