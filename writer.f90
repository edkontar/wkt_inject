MODULE writer

IMPLICIT NONE

CONTAINS

Subroutine  WriteEnergy (fv,wv,W_t1,W_t2,Tstep)
! The procedure is intended to write profiles
! of electron distribution function and  
! spectral energy density of Langmuir waves to the 
! respective disk files

  USE constant
! General definition of constants
  USE params 
! Load certain global variables

  IMPLICIT NONE

  DOUBLE PRECISION, dimension (-Nv:Nv), intent (in) :: fv, wv, W_t1,W_t2;
 

  integer :: WriteStat, filename  
! Service variables 
  integer, intent(in)::  Tstep
  DOUBLE PRECISION::E_beam,n_beam,E_L, E_F,E_H,dk,k;

  DOUBLE PRECISION ::E_beam_0, E_total,Sminus,Splus,Lminus,Lplus;

  character (LEN= 32):: E_flush;
 
! Findex - the output file numbers

IF(Tstep ==0) THEN
open (UNIT=19, FILE =e_sequence, STATUS = 'REPLACE', ACTION='WRITE')
ELSE
open (UNIT=19, FILE =e_sequence, STATUS = 'OLD',POSITION='APPEND', ACTION='WRITE')
END IF
! open file

E_beam   = 0.
n_beam   = 0.
E_total  = 0.
E_L      = 0.
E_F      = 0.
E_H      = 0.
Lplus=0.
Lminus=0.
Splus=0.
Sminus=0.


DO i=-Nv+1,Nv
          
 dk=abs(kx(i)-kx(i-1))*Omega/v_T
 dv=abs(vx(i)-vx(i-1))
 v=vx(i)
 
 k =kx(i)*Omega/V_t
 
 IF((abs(v)>4.*v_t).AND.(i>1)) THEN
 n_beam = n_beam + dv*(Fv(i)+Fv(i-1))*k1/(2.*nb)
 E_beam = E_beam + dv*(Fv(i)*vx(i)**2+Fv(i-1)*vx(i-1)**2)*k1/(4.*Eb)
 END IF
 
 E_L    = E_L +dk*(Wv(i)+Wv(i-1))*k2/(2.*Eb*m)
 E_F    = E_F +dk*(W_t1(i)+W_t1(i-1))*k2/(2.*Eb*m)
 E_H    = E_H +dk*(W_t2(i)+W_t2(i-1))*k2/(2.*Eb*m)

 IF (i<0) Lminus=Lminus+dk*(Wv(i)+Wv(i-1))*k2/(2.*Eb*m)
 IF (i>0) Lplus =Lplus+dk*(Wv(i)+Wv(i-1))*k2/(2.*Eb*m)
 IF (i<0) Sminus=Sminus+dk*(W_t2(i)+W_t2(i-1))*k2*1000./(2.*Eb*m)
 IF (i>0) Splus =Splus+dk*(W_t2(i)+W_t2(i-1))*k2*1000./(2.*Eb*m)
END DO
  


 E_total =E_beam+E_L+E_F+E_H
 
 if (abs(E_beam) <10d-99) E_beam  =0.
 if (abs(n_beam) <10d-99) n_beam  =0.
 if (abs(E_L)    <10d-99) E_L     =0.
 if (abs(E_F)    <10d-99) E_F     =0.
 if (abs(E_H)    <10d-99) E_H     =0.
 if (abs(E_total)<10d-99) E_total =0.

write (19,'(E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3)', IOSTAT= WriteStat )&
t,' ',n_beam,' ',E_beam,' ',E_L,' ',E_F,' ',E_H,' ',E_total,' ',Lminus,' ',Lplus,' ',Sminus,' ',Splus
  If(WriteStat /= 0) stop 'Error while writing F and W'

close (19)
write (*,'(A,E10.3,A,E10.3,A,E10.3,A,E10.3)')'Energy (E,L,S,Total) =>',E_beam,' ',E_L,' ',E_F,' ',E_total
write (*,'(A,ES10.3)')'Number (Electrons) =>',N_beam
! Close files

end subroutine  WriteEnergy

!-----------------------------------------------------

Subroutine  WriteProfiles (fv,wv,W_t1,W_t2,Tstep)
! The procedure is intended to write profiles
! of electron distribution function and  
! spectral energy density of Langmuir waves to the 
! respective disk files

  USE constant
! General definition of constants
  USE params 
! Load certain global variables

  IMPLICIT NONE

  DOUBLE PRECISION, dimension (-Nv:Nv), intent (in) :: fv, wv, W_t1,W_t2;
 

  integer :: WriteStat, filename  
! Service variables 
  integer, intent(in)::  Tstep
  DOUBLE PRECISION ::N_beam, Norm_x, Norm_v, Norm_fxyz, Norm_wxyz, Norm_dens,Norm_wt1, Norm_wt2;
  DOUBLE PRECISION ::Norm_eb, Norm_ew, w_norm;

  character (LEN= 32):: w_v_flush;
  character (LEN= 5) :: Findex
! Findex - the output file numbers


  filename =Tstep
  write (Findex,'(i5.5)') filename
  Findex=trim(Findex)
  w_v_flush= trim(w_v_sequence//Findex//'.dat')
! this file contains all data

  open (UNIT=12, FILE =w_v_flush, STATUS = 'REPLACE', ACTION='WRITE')

! open file with density profile

  do i=-Nv,Nv

	 Norm_v     = vx(i)/(v_T)
	 Norm_Fxyz  = Fv(i)*k1*vbeam/nb
	 Norm_Wxyz  = Wv(i)*k2*omega/(m*nb*vbeam**3)

         Norm_Wt1   = W_t1(i)*k2*omega/(m*nb*vbeam**3)
         Norm_Wt2   = W_t2(i)*k2*omega/(m*nb*vbeam**3)
	 	 
     if (abs(Norm_Fxyz)<10d-99) Norm_Fxyz=0.0
     if (abs(Norm_Wxyz)<10d-99) Norm_Wxyz=0.0
     if (abs(Norm_Wt1) <10d-99) Norm_Wt1 =0.0 
     if (abs(Norm_Wt2) <10d-99) Norm_Wt2 =0.0        	  
 
write (12,'(E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3)', IOSTAT = WriteStat )&
kx(i),' ',Norm_v,' ',Norm_Fxyz,' ',Norm_Wxyz,' ',Norm_Wt1,' ',Norm_Wt2

  If(WriteStat /= 0) stop 'Error while writing F and W'

  end do

  close (12)
! Close files 
     
end subroutine  WriteProfiles

END MODULE WRITER