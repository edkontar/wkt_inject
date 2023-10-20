MODULE NONLIN


CONTAINS

SUBROUTINE ION_SCATT(W_L,dt_ion)

 ! Langmuir wave scattering off ions

USE CONSTANT
USE PARAMS

   IMPLICIT NONE

   DOUBLE PRECISION, DIMENSION(-Nv:Nv),INTENT(INOUT) :: w_L;
   DOUBLE PRECISION, INTENT(INOUT) :: dt_ion;
   ! loval variables are below
   DOUBLE PRECISION, DIMENSION(-Nv:Nv):: a_ion;
   INTEGER:: ik;
   DOUBLE PRECISION :: k,k_1,dk,Omega_1,ksi,Omega_0, alpha,k_m;
!----------------------------------------------------

Alpha = sqrt(2*Pi)*Omega**2/(4.*(2.*Pi)**3*ne*v_Ti*(1+T_e/T_i)**2)

   a_ion =0. 


DO j=-Nv+1, Nv-1   

   k  = kx(j)*Omega/v_t
   !k vector of integration
   
   Omega_0 =Omega+3.*v_T**2*k**2/(2.*Omega)
  
   DO ik=-Nv+1,Nv-1
   
   k_1= kx(ik)*Omega/v_T
   dk=abs(kx(ik)-kx(ik-1))*Omega/v_T
   k_m =abs(k-k_1)

   Omega_1  =Omega+3.*v_T**2*k_1**2/(2.*Omega)
   !frequency of a scattered wave

   IF (ik==j) THEN
   a_ion(j)=a_ion(j)-dk*exp(-ksi**2)*(2.*Pi)**3*6.*v_t**2*k*W_L(j)*W_L(ik)*k2&
   /(k_b*T_i*Omega**2)
   ELSE
   ksi  = (3.*V_T*v_T*(k+k_1))/(2.*sqrt(2.)*v_Ti*Omega)    
   a_ion(j)=a_ion(j)+dk*exp(-ksi**2)*(Omega_0*W_L(ik)/Omega_1-W_L(j)-&
   (2.*Pi)**3*(Omega_0-Omega_1)*W_L(j)*W_L(ik)*k2/(k_b*T_i*Omega_1))/k_m
   END IF 
  
   END DO
END DO

W_L=W_L+alpha*a_ion*dt_ion

dt_ion =0.2/(Alpha*MaxVal(abs(a_ion)/W_L))
! scattering from ions

if(dt_ion>Tqv)  dt_ion=Tqv
if(dt_ion<dt)   dt_ion=dt

END SUBROUTINE ION_SCATT

SUBROUTINE DECAY(W_L,W_S,dt_decay)

 ! 3waves decay L->L+S
 
USE CONSTANT
USE PARAMS

   IMPLICIT NONE

   DOUBLE PRECISION, DIMENSION(-Nv:Nv),INTENT(INOUT) :: W_L,W_S;
   DOUBLE PRECISION, INTENT(INOUT) :: dt_decay;
   !below are local variables
   DOUBLE PRECISION, DIMENSION(-Nv:Nv)::N_L,N_S,a_S,a_L,gamma_s;
   INTEGER::kp
   DOUBLE PRECISION :: alpha_is,k_prime;
!----------------------------------------------------
alpha_is=Pi*Omega**3*v_s/(6.0*ne*k_b*T_e*v_T**2)

gamma_s=sqrt(Pi/2.)*abs(kx)*v_s*Omega/v_t*(V_s/V_t+&
(v_s/v_Ti)**3*exp(-v_s*v_s/(2.*v_Ti*v_Ti)))

k_prime =v_s*Omega/(3.*v_t**2)

kp=floor(1.+(2.*k_prime*v_t/Omega)/abs(kx(1)-kx(2)))


N_L=W_L/(Omega+3.*Omega*kx**2/2.)
WHERE (kx/=0.) N_s=W_s/(v_s*abs(kx)*Omega/v_t)
! the number of S and L waves

FORALL(j=1:Nv) a_s(j)=N_L((j+kp)/2)*N_L((-j+kp)/2)+N_s(j)*(N_L((j+kp)/2)-N_L((-j+kp)/2))

FORALL(j=-Nv:-1) a_s(j)=N_L((j-kp)/2)*N_L((-j-kp)/2)+N_s(j)*(N_L((j-kp)/2)-N_L((-j-kp)/2))

FORALL(j=kp/2:(Nv-kp)/2) a_L(j)=N_s(2*j-kp)*(N_L(-j+kp)-N_L(j))&
+N_L(j)*(N_L(-j-kp)-N_L(-j+kp))+N_s(-2*j-kp)*(N_L(-j-kp)-N_L(j))
    
FORALL(j=-kp/2:kp/2) a_L(j)=N_L(j)*(N_L(-j+kp)+N_L(-j-kp))&
+N_S(-2*j+kp)*(N_L(-j+kp)-N_L(j))+N_s(-2*j-kp)*(N_L(-j-kp)-N_L(j))

FORALL(j=-(Nv-kp)/2:-kp/2) a_L(j)=N_s(2*j+kp)*(N_L(-j-kp)-N_L(j))&
+N_L(j)*(N_L(-j+kp)-N_L(-j-kp))+N_s(-2*j+kp)*(N_L(-j+kp)-N_L(j))
    
FORALL(j=kp/2:kp/2) a_L(j)=N_L(j)*N_L(-j-kp)+N_S(-2*j-kp)*(N_L(-j-kp)-N_L(j))

FORALL(j=-kp/2:-kp/2) a_L(j)=N_L(j)*N_L(-j+kp)+N_S(-2*j+kp)*(N_L(-j+kp)-N_L(j))
!end calculation of coeficients for N_L


W_L=W_L+a_L*alpha_is*k2*Omega*(1+3.*kx**2/2.)*dt_decay
W_S=W_S-W_s*gamma_s*dt_decay+0.5*a_s*alpha_is*k2*v_s*Omega*abs(kx)*dt_decay/v_T

where(W_s<0)W_s=0.

dt_decay=0.2/(alpha_is*k2*MaxVal(a_L*Omega*(1+3.*kx**2/2.)/W_L))

if(dt_decay>Tqv)      dt_decay=Tqv
if(dt_decay<dt)       dt_decay=dt

END SUBROUTINE DECAY

END MODULE NONLIN