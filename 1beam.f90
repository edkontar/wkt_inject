PROGRAM QBEAM_PLASMA

  ! This program solves the one-dimensional kinetic equations
  ! of quasilinear relaxation using implicit
  ! difference scheme (second order approximation over coordinate
  ! and velocity, and first order approximation over time)

  ! Created by :
  ! Eduard@ Glasgow 2009
  ! Eduard@ Glasgow, April 2023: minor changes to upload to github
  ! Eduard@ Glasgow, April 2023: changed alog to a log 
  !
 
  USE CONSTANT

  USE PARAMS  

  USE READER

  USE WRITER  
  
  USE INIT_BEAM

  USE NONLIN

  IMPLICIT NONE 

    integer:: step,  Time_control, Num_time_steps,NumPasses, WriteStat;
	! sevice variables
	DOUBLE PRECISION :: Wmax,Grad;
	! maximum of spectral energy density

    DOUBLE PRECISION:: time_end, time_begin,t_scat;
	! time diagnostic variables
    DOUBLE PRECISION::  remain_time, time_est, last_save,dt_deay,Dt_ion,LLL,dkx,dk;
	!service variables

    DOUBLE PRECISION::  t_deay,dt_decay,t_ion,t_decay,coll_ei,lnC;
	!service variables

    integer(4) :: itime, rate, max
	! CPU_time  variables
      
    DOUBLE PRECISION, DIMENSION(-Nv:Nv):: F=0.0, W=0.0, Dcol=0.0, W_F,W_H,gamma_L;
	! Electron distribution function and
	! Spectral energy density
    
    DOUBLE PRECISION, DIMENSION(nx):: E_Density=0.0;
	! Electron number density
    
!------------  Main plasma values definition --------------------
last_save=0.
T_scat=0.
v_T = sqrt(k_b*T_e/m);
v_Ti= sqrt(k_b*T_i/m_p);
vmin= V_T;
v_s=sqrt(k_b*T_e/m_p)*sqrt(1+3.*T_i/T_e);

  CALL ReadCONFIG;
  ! reading initial configuration 

  WRITE(*,*) "Initial configuration      -       OK";
  write(*, '(A,F5.1,A)')'Total time of calculations =', tmax - tmin, ' sec'

  ne=m*OMEGA*OMEGA/(4.*pi*e*e);

  !defines background plasma density

  Grad=V_t
  ;coll_ei=1d2

  lnC=20.
  coll_ei=4.*PI*e**4*ne*lnC/(m*m*V_T*V_T*V_T)


!local plasma gradient

! --- calculations of coordinate and velocity steps ---

    !dv=(vmax1 - vmin)/(nv1-1)
    !dv2=(vmax-vmax1)/(nv-nv1-1)
    
    !dk=2./float((2*Nv-1))
    
     !dx=(xmax - xmin)/(nx-1)

DO i=-Nv,Nv
kx(i)=0.5*float(i)/float(Nv)
if (i/=0) vx(i)=v_t/kx(i)
!IF (i.LE.Nv1) THEN
!kx(i)= vmin+dv*(i-1)
!kx(i)= v_t/vx(i)
!vx(2*nv-i+1)= -vx(i)
!kx(2*nv-i+1)= -kx(i)
!ELSE
!vx(i)= vx(i-1)+dv2*(i-nv1)
!kx(i)= v_t/vx(i)
!vx(2*nv-i+1)= -vx(i)
!kx(2*nv-i+1)= -kx(i)
!END IF

END DO 


dk=abs(kx(1)-kx(2))

DO i=-Nv,Nv
if (i/=0) Dcol(i)=coll_ei*kx(i)*kx(i)
END DO 
!Dcol=coll_ei*kx*kx
! collisional difusion in energy space


Nb=0.
Eb=0.

  DO  i=2,Nv
  
  v = vx(i)
  dv= abs(vx(i)-vx(i-1))
  
  IF(abs(v)>(4.*v_t)) THEN
  nb = nb + (F0(v)+F0(vx(i-1)))*dv/2.0
  Eb = Eb + (F0(v)*v*v+F0(vx(i-1))*vx(i-1)**2)*dv/4.
  END IF
  
  END DO

  ! calculation of electron beam density and 
  ! energy energy of electron beam

! --- calculation of dimensionless coeficients starts here ----
  
   ! k1=2*pi*vp*vp*vp/(vd*vd*vd*vd);
   ! k2=m*vd*vd/(2*pi*pi*pi)

   tqv=ne/(omega*nbeam);  
   ! the quasilinear time

   k1 = nbeam/vbeam;
   ! electron distribution dimensionl parameter
   k2 = m*nbeam*Vbeam**3/Omega
   ! spectral energy density dimensional parameter
    
   !a1=pi*OMEGA*k1/ne

    a1 = Vbeam **3 *PI/Tqv

    !spont_f   =  omega/(4.*PI*ne);
    !spont_w  =  spont_f/Vbeam**4 *V_T*V_T/omega*omega *k1/k2;

! 1D spontaneous terms by Olena
   spont_w =omega**3*m*V_T*k1/k2/(4.*PI*ne)

   ! First equation koeficients
  
    !a2=(4*pi*pi*e*e*k2/m)/m;     

     a2 = Pi/(Tqv*Vbeam)
    ! Second equation Koeficient
     
    !a_s =Pi*Omega**2*k2/(18.*m_p*ne*(v_t)**4)

!------------ Introduction of  values to calculate  ---------------------//
     
      t=tmin;
      DO i=-Nv,Nv
      v=vx(i)
      IF(i/=0) F(i) = F0(v)/k1
	  !
      !F(i) =0. ! starting without electrons, electrons are injected with time.
	  !W(i) =1e-12
      W(i) =0.
      !If ((v>0.).AND.(v<=vbeam)) W(i) =(v/vbeam)**4*(1-v/vbeam)
      W_H(i)=1e-16
      !W_H(i)=0.
      !if((kx(i)<0.3) .and. (kx(i)>0.2)) w_h(i)=-10.*(kx(i)-0.2)*(kx(i)-0.3)
      ! Spectral energy density of the thermal level
      If(i/=0) gamma_L(i) =2.*sqrt(Pi/8.)*Omega*exp(-0.5/kx(i)**2)/abs(kx(i))**3
      END DO
! --------BOUNDARY CONDITIONS IN VELOCITY SPACE ------------------

!F(Nv) =1.
!F(Nv-1) =1.
!W(Nv) =0.
!W(-NV)=0.


W_H(0) =0.
W_H(Nv)=0.
W_H(-Nv)=0.

! ---------------- Time step calculations ----------------

   ! according to stability criteria of the first equation  
   
         Wmax = 1
	 !dt = dv*tqv/(1000.0*PI*vbeam)
        dt= dv*tqv/(100.0*PI*vbeam)
		  
!dt  = Tqv*abs(kx(1)-kx(2))*v_T/(1.*pi*vbeam*MAXVAL(F))

dt=0.01*dk*dk*Tqv/PI*(vbeam/v_t)*(vbeam/v_t)

!write(*,*) dt
!stop;

dt = 1.01e-9

dt_decay=1.*dt
dt_ion=1.*dt

dk =abs(kx(1)-kx(2))

!write(*,*)dt2
! -----------------  data output to the screen ------------
   write(*,*) '------------------------------------------'
   write(*,*) 'Finite diference scheme parameters :' 
   write(*,'(A,F14.12,A)')'Time step        dt=', dt,  ' sec'
   write(*,'(A,F6.4,A)')'X-coodinate step dx=', dx/d,'*d cm'
   write(*,'(A,F6.4,A)')'Velocity step    dv=', dv/vbeam,'*Vbeam cm/sec'
   write(*,*) '------------------------------------------'
   Write(*,*)'ELECTRON BEAM & PLASMA PARAMETERS:'
   write(*,'(A,F5.1,A)')'Plasma frequency =', Omega/(2*PI*1e9), ' GHz'
   write(*,'(A,E10.3,A)')'Beam density     =', nb, ' cm^{-3}'
   write(*,'(A,E10.3,A)')'Plasma density   =', ne, 'cm^{-3}'
   write(*,'(A,E10.3,A)')'Beam_density/plasma_density=', nbeam/ne,'  '
   write(*,'(A,E9.2,A)')'Pl Freq. Omega_pe=', Omega,' s^{-1}'
   write(*,'(A,E9.2,A)')'Collisional Freq =', coll_ei,' s^{-1}'
   write(*,'(A,E10.3,A)')'Quasilinear time =', tqv,' sec'

! ----------------------------------------------
CALL WriteProfiles(F,W,W_F,W_H,0)
CALL WriteEnergy  (F,W,W_F,W_H,0)

! Initial configuration is saved
! ------------------------------------------------
	 
! ------- main loop parameters ----------------------------

	   step=0;	  	   
       Time_flush   = Time_save/dt;
       Time_control = 10;
       Num_time_steps = (tmax-tmin)/dt;

   write(*,*) '------------------------------------------'
   write(*,*) "Calculations started    -    OK";

   call system_clock(count=itime, count_rate=rate, count_max=max)
   time_begin = Float(itime)/rate


! MAIN LOOP starts here 

    DO WHILE (t<tmax)

!dt=MIN(dv*tqv/(1000.0*PI*vbeam), 0.01*dv/(a_s*vbeam*vbeam*MaxVal(w)))

!write(*,*) t, dt
	  t = t + dt
!step = step + 1

IF ((t-last_save).GE.time_save) THEN
                 
                     CALL WriteProfiles(F,W,W_F,W_H,FLOOR(t/time_save))
                     CALL WriteEnergy (F,W,W_F,W_H,FLOOR(t/time_save))
                     last_save = time_save*FLOOR(t/time_save)
                     write(*,'(A,F9.5,A)') 'All data saved at T=',t, ' sec'
                     call Check_time
                    
		     !saving of data to disk

     END IF

! ----------- BASIC CALCULATION ---------------------------- 

!FORALL(i=2:floor(0.3/dk-1)) f(i)=f(i)+dt*5e6*kx(i)**4/v_t*(0.05-t/kx(i))*exp(-(0.05-t/kx(i))**2/0.01**2)/k1/0.01

FORALL(i=2:floor(0.3/dk-1)) f(i)=f(i)+dt*Nbeam*exp(-(1./kx(i)-25.)**2/10.**2)*&
                                      exp(-(1.-t)**2/0.15**2)/k1/0.15/10./PI
!injection of electrons 
! change made after discussion with Adam K. and Meriem 
! Basically WKT code is modified similar to equations 1 and 2 in 
! https://ui.adsabs.harvard.edu/abs/2008AnGeo..26.2435S/abstract
!Note that factor is still arbitrary and needs fixing

FORALL(i=2:floor(0.4/dk-1)) w(i) = w(i)&
-a2*W(i)*(f(i+1)-f(i))*v_t*dt/(dk)+spont_w*f(i)/abs(kx(i))*log(1./abs(kx(i)))*dt
!- gamma_L(i)*W(i)*dt +spont_w*f(i)/abs(kx(i))*log(1./abs(kx(i)))*dt

FORALL(i=-floor(0.4/dk-1):-2) w(i) = w(i)&
+a2*W(i)*(f(i+1)-f(i))*v_t*dt/(dk) +spont_w*f(i)/abs(kx(i))*log(1./abs(kx(i)))*dt
! - gamma_L(i)*W(i)*dt+spont_w*f(i)/abs(kx(i))*log(1./abs(kx(i)))*dt

FORALL(i=-Nv+1:Nv-1)  w(i) = w(i)+v_t*(W(i+1)-W(i))*dt/(dk*1e+7)
! inhomogeneity influence 

FORALL(i=-Nv+1:Nv-1)  w(i) = w(i)-coll_ei*dt*w(i)/4.
! collisional absorption of Langmuir waves 
where( W < 1e-29 ) W=1e-29

FORALL(i=2:floor(0.4/dk-1)) f(i) =f(i) + a1*dt*((w(i)*kx(i)**3)*(f(i+1)-f(i))-&
(w(i-1)*kx(i-1)**3)*(f(i)-f(i-1)))*kx(i)*kx(i)/(dk*dk*v_t**3)

FORALL(i=-floor(0.4/dk-1):-2) f(i) =f(i) - a1*dt*((w(i)*kx(i)**3)*(f(i+1)-f(i))-&
(w(i-1)*kx(i-1)**3)*(f(i)-f(i-1)))*kx(i)*kx(i)/(dk*dk*v_t**3)

!FORALL(i=2:floor(0.4/dk-1)) f(i)=f(i)+dt*1e-1/k1/3.*(&
!exp(-(15.-1./kx(i))**2/3.**2)+1e6*(kx(i)/(1.+5.*kx(i)))**8)

!FORALL(i=2:floor(0.4/dk-1)) f(i)=f(i)+dt/k1*1e12*exp(-((.05-t/kx(i))/0.01)**2)*((.05-t/kx(i))/kx(i))*(kx(i)/(1.+10.*kx(i)))**8
! mimicking arrival and departure of particles into agiven point

!FORALL(i=2:floor(0.4/dk-1)) f(i)=f(i)+dt/k1*1e9*exp(-((.05-t/kx(i))/0.01)**2)*((.05)/kx(i))*(kx(i)/(1.+10.*kx(i)))**6
!FORALL(i=2:floor(0.4/dk-1))     f(i)=f(i)-f(i)*dt*(1e2/kx(i)) ! leaving particles
!FORALL(i=-floor(0.4/dk-1):-2)  f(i)=f(i)-f(i)*dt*(1e1/(kx(i))) ! leaving particles-opposite direction

!collisions for electrons
!FORALL(i=2:floor(0.4/dk-1)) f(i)=f(i)-dt*coll_ei*exp(-25.*kx(i)**2)*kx(i)*kx(i)*(F(i)*kx(i)*kx(i)-F(i-1)*kx(i-1)*kx(i-1))/dk
!FORALL(i=-floor(0.4/dk-1):-2) f(i)=f(i)+dt*coll_ei*exp(-25.*kx(i)**2)*kx(i)*kx(i)*(F(i+1)*kx(i+1)*kx(i+1)-F(i)*kx(i)*kx(i))/dk

FORALL(i=2:floor(0.4/dk-1)) f(i)=f(i) - dt*coll_ei*kx(i)*kx(i)*(F(i)*kx(i)*kx(i)-F(i-1)*kx(i-1)*kx(i-1))/dk

FORALL(i=-floor(0.4/dk-1):-2) f(i)=f(i)+dt*coll_ei*kx(i)*kx(i)*(F(i+1)*kx(i+1)*kx(i+1)-F(i)*kx(i)*kx(i))/dk

FORALL(i=2:floor(0.4/dk-1)) f(i) =f(i) + coll_ei*dt*((kx(i)**5)*(f(i+1)-f(i))-&
(kx(i-1)**5)*(f(i)-f(i-1)))*kx(i)*kx(i)/(dk*dk)

FORALL(i=-floor(0.4/dk-1):-2) f(i) =f(i) - coll_ei*dt*((kx(i)**5)*(f(i+1)-f(i))-&
(kx(i-1)**5)*(f(i)-f(i-1)))*kx(i)*kx(i)/(dk*dk)

! collisions stop here

where (F < 1e-20) F=1e-20

where( W < 1e-29 ) W=1e-29

where (abs(kx)>0.3) W=1e-29


!dt=MIN(0.01*(dk*dk*v_t**3)/(a1*MaxVal(W*kx*kx*kx*kx*abs(kx))),&
!0.01*dk/(a2*MaxVal(f)*v_t))
!write(*,*) dt
! ------------- ion scaterring ----------------------
 IF ((t-t_ion).GE.dt_ion) THEN
     dt_ion=t-t_ion
     !CALL ION_SCATT(W,dt_ion)
     t_ion = t
     !write(*,*)dt_ion
 END IF
! ----------- end ion scaterring ------------------- 

! ------------- 3 wave decay ----------------------
 IF ((t-t_decay).GE.dt_decay) THEN
     dt_decay=t-t_decay
     CALL DECAY(W,W_H,dt_decay)
     t_decay = t
     dt_decay =dt
!write(*,*) "heee"
 END IF
! ----------- end 3 wave decay  ------------------- 
	     
END DO
! End of main loop

WRITE(*,*) "CALCULATIONS COMPLETED     -    OK !"


CONTAINS
 
!-------------------------------------------------------------------------------------- 

 subroutine Check_time
 ! subroutine that calculates remaining time 

 IMPLICIT NONE

 DOUBLE PRECISION:: av_dens, timeav, timeacc;

    if ( Mod (Step, Time_control) == 0)   then
    !Check whether control event has happened and display current time 
       
  call system_clock(count=itime, count_rate=rate, count_max=max)
  
     time_end = float(itime)/rate
	 
	 time_est=time_end - time_begin

	 if (time_est<0) time_est = time_est + max/rate
     ! condition checks if counter started from zero again 
                                								
		remain_time= (Tmax-t)*time_est/t
!           Find remaining time   

Write (* , '(A,F6.3,A, F4.1,A,F6.2,A,I5)') 'Passed T=',t ,'sec (',t*100/(tmax-tmin),'% comleted) for ', time_est,' seconds/',&
 Time_control
    !Write (* , *) 'Passed T=',t ,'sec (',t*100/(tmax-tmin),'% comleted) for ', time_est,' seconds/',Time_control
	open (UNIT=33, FILE = 'progress.ttt', STATUS = 'REPLACE', ACTION='WRITE' )   
	
	
Write (33 , '(A,F6.3,A, F4.1, A, F6.2,A,I5)',IOSTAT = WriteStat) 'Passed T=',t ,'sec (',t*100/(tmax-tmin),'% comleted) for ',&
 time_est,' seconds/', Time_control

    If (WriteStat /= 0) stop 'Error while writing progress.ttt file'
    
	CALL WriteHRS(remain_time)
    
	close(33);

         						
!   Translate remaining time to hrs.mm.sec and write to screen
 	 time_begin = time_end

     end if
	   
end subroutine Check_time

!-----------------------------------------------------------

 subroutine WriteHRS(time0)
 !Translating remaining time to hrs.mm.sec and write to screen
 ! and progress file
    
      DOUBLE PRECISION, intent(in):: time0
      DOUBLE PRECISION :: timex

      integer:: hrs,mins,secs
      
	  timex=time0
      hrs=Floor(timex/3600)
      timex=timex-hrs*3600

      mins=Floor(timex/60)
      timex=timex-mins*60

      secs=timex
  
  Write(*,'(A,I5,A, I2, A, I2, A)')  'Remaining time:', hrs,' hrs ', mins,' mins ',secs,' sec.' 

  Write(33,'(A,I5,A, I2, A, I2, A)',IOSTAT = WriteStat) 'Remaining time:', hrs,' hrs ', mins,' mins ',secs,' sec.' 

  If (WriteStat /= 0) stop 'Error while writing REMAINING TIME to progress.ttt file'
  
  end subroutine WriteHRS

!-----------------------------------------------------------

END PROGRAM
