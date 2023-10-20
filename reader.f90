MODULE READER

 ! Module contains the subroutine reading
 ! initial configuration of the system 

 USE CONSTANT

 USE PARAMS

 IMPLICIT NONE

 CONTAINS

 SUBROUTINE ReadConfig
 
 IMPLICIT NONE

 INTEGER::ReadStatus, iacc
 ! service variable to control read status

 iacc=0

   OPEN (UNIT=10, FILE= init_params, STATUS='OLD', ACTION='READ' ) 
   ! open file for reading
 
	  READ (10,*, IOSTAT= ReadStatus ) Tmax
		If (ReadStatus /= 0) iacc=iacc+1
	   
          READ (10,*, IOSTAT=ReadStatus ) Nbeam
		If (ReadStatus /= 0) iacc=iacc+1		
       
          READ (10,*, IOSTAT=ReadStatus ) NuPlasma
		If (ReadStatus /= 0) iacc=iacc+1		
        
		NuPlasma = NuPlasma * 1.0e+6
		OMEGA= 2.0*Pi*NuPlasma
	    
	  READ (10,*, IOSTAT=ReadStatus ) d
		If (ReadStatus /= 0) iacc=iacc+1		
       
	        d = d * 1.0e+9

	  READ (10,*, IOSTAT=ReadStatus ) Vbeam
		If (ReadStatus /= 0) iacc=iacc+1		
       
	        Vbeam = Vbeam * 1.0e+9

	  READ (10,*, IOSTAT=ReadStatus ) Xmin
		If (ReadStatus /= 0) iacc=iacc+1		
       
                Xmin = Xmin * 1.0e+9

          READ (10,*, IOSTAT=ReadStatus ) Xmax
		If (ReadStatus /= 0) iacc=iacc+1		
      
                Xmax = Xmax * 1.0e+9

	  READ (10,*, IOSTAT=ReadStatus ) Time_save
		If (ReadStatus /= 0) iacc=iacc+1		
      
	 IF (iacc> 0) THEN
         ! check reading status
		
         Write (*,*) 'Error  reading file :', init_params
         STOP 
		
	 END IF		 
		            
	 CLOSE(10)
         ! close file 				        

 END SUBROUTINE ReadConfig

END MODULE READER