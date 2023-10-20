
pro  fw_v
; program to plot 1D images of
;density
;energy density of electrons and plasmons

nb=1e8
vb=5.e+9

sc_fv=[4.4e-8,5.1e-8,4.6e-8,1.5e-8,3.3e-9,9.3e-11,2e-9]
sc_v=[6.3,7.6,9.2,11.2,13.5,16.3]
col_ei=3e3 ; 1/s

Nk=1001
;the file to show

;i=FileNumber


;DEVICE, RETAIN=2,  DECOMPOSED =0
DEVICE, GET_DECOMPOSED=orig_decomposed, DECOMPOSED=0

LOADCT, 39  ;rainbow +white
; myMPEG = OBJ_NEW('IDLgrMPEG',  FILENAME='FLS.mpeg')

oVid = IDLffVideoWrite('eLS.mp4')

width =400
height=600
fps = 20
vidStream = oVid.AddVideoStream(width, height, fps)

FileNumber=n_elements(FILE_SEARCH('vfw*.dat'))
FOR i=0, FileNumber-1,1 DO begin
;  FOR i=0, 80 ,1 DO begin
;i=FileNumber

data_file    = 'vfw' +string(format='(I5.5)', i*1,/print)+'.dat'

openr,1, data_file

A = fltarr(6,Nk)
Vx=fltarr(Nk)
Felectron=fltarr(Nk)
Velocity=fltarr(Nk)
; describtion of the file content

READF, 1,  A

CLOSE, 1

A = transpose(A)
; conversion of data
vel  = A(*,1)
Kx   = A(*,0)
F_v  = A(*,2)
w_v  = A(*,3)
w_F  = A(*,4)
w_S  = A(*,5)


F_v=F_v*nb/vb


IF (i EQ 0) THEN BEGIN
WL_thermal=W_v
WS_thermal=W_s
F_V0=F_v
F_total=F_v*1e-10
ENDIF

;vel=1./kx

;FOR j=0,Nk/2-1 DO  BEGIN
;Felectron(j)=F_v(Nk/2-j-1)
;Velocity(j)=vel(Nk/2-j-1)
;Vx(j)=-2.-float(Nk/2-j)*40./float(Nk)
;END


;FOR j=0,Nk/2-1 DO  BEGIN
;Felectron(j+Nk/2)=F_v(Nk-j-1)
;Velocity(j+Nk/2)=vel(Nk-j-1)
;Vx(j+Nk/2)=2+float(j)*40/float(Nk)
;END

;Felectron=ALOG10(1+Felectron)
;printf, 33,  float(i)/1000.0 ,  u_t(i)

WL_thermal=1e-10

window, 0, xsize=400, ysize=600, TITLE='Distribution of Electrons & Langmuir waves'

!p.charsize=2
!P.thick=1.2
!P.MULTI=[0,1,3]

;Felectron =INTERPOL(Felectron, Velocity, Vx, /SPLINE)
;W_v =INTERPOL(W_v,kx,kx, /SPLINE)
;W_s =INTERPOL(W_s,kx,kx, /SPLINE)

;plot,Vel, (F_v)/sqrt(3.1416), /ylog,xrange=[-15.,15.], xticks=6,Thick=2,yrange=[1e-12,1e-6],BACKGROUND=255, COLOR=0, XTITLE ='V/V_beam', YTITLE ='Electrons F(V)'
;oplot,sc_v,sc_fv,line=1

;PLOT, Felectron
;F_v =INTERPOL(Felectron, Nk, /SPLINE)
;FVT(*,i)=ALOG10(1+felectron)

;plot,kx , W_v/WL_thermal,/ylog, ystyle=1,xstyle=1,xrange=[-0.5,0.5],xticks=10,thick=2,yrange=[1,1e6],BACKGROUND=255, COLOR=245, XTITLE ='k*d_De', YTITLE =' Langmuir W(k)'

;plot,kx , W_S/WS_thermal, /ylog,ystyle=1,xstyle=1,xrange=[-0.5,0.5],xticks=10,thick=2,yrange=[1,1e6],BACKGROUND=255, COLOR=50, XTITLE ='k*d_De', YTITLE ='Ion-sound Ws(k)x1E-7'

plot,Vel, F_v/sqrt(3.1416), xrange=[-40.,40.],xticks=6,Thick=2,/ylog,yrange=[1e-10,1], $
  XTITLE ='Velocity, v/v!BTe!N', YTITLE ='Electrons F(V)',BACKGROUND=255, COLOR=0,/xs,/ys
;oplot,sc_v,sc_fv,line=2,thick=2,color=245
oplot,vel,F_v0/sqrt(!PI),line=1,color=245
velp=vel(where(vel >0))
;oplot,velp,1e8*(1./(10.+velp))^8/sqrt(3.1416),line=2,color=45

plot,kx , W_v/WL_thermal,/ylog, ystyle=1,xstyle=1,xrange=[-0.5,0.5],xticks=10,thick=2,yrange=[1,1e14],XTITLE ='Wavenumber, k!4k!3!BDe!N', YTITLE ='Langmuir W!BL!N(k)/W!BT!N',BACKGROUND=255,COLOR=230

plot,kx , W_S/WS_thermal, /ylog,ystyle=1,xstyle=1,xrange=[-0.5,0.5],xticks=10,thick=2,yrange=[1.,1e6],XTITLE ='Wavenumber, k!4k!3!BDe!N', YTITLE ='Sound W!BS!N(k)/W!BT!N',BACKGROUND=255,COLOR=50

XYOUTS ,240,582,COLOR=50, string(format='(A, F5.3, A)', 'Time=', float(i)/(1e+3), 's', /print ), /device

!p.charsize=1
!P.thick=1
!P.MULTI=0

wait, 0.5

frame=TVRD(/true)
!NULL = oVid.Put(vidStream,frame )

image3d = TVRD(True=1)
image2d = Color_Quan(image3d, 1, r, g, b,CUBE=6)
;Write_GIF, 'Ftest.gif', image2d, r, g, b
IF N_ELEMENTS(ImageCube) LE 1 THEN ImageCube=image2D
imageCube=[[[imageCube]],[[image2D]]]
;write_gif,'FLS.gif', image2d, r, g, b, /multiple

;stop

;WDELETE,0
;window,0
;tv,image1,0,true=3
;stop

F_total=(F_total+F_v)
;window,2
;plot,Vel, F_total/sqrt(3.1416)/float(i), xrange=[-60.,60.],xticks=6,Thick=2,/ylog,yrange=[1e-10,1], XTITLE ='Velocity, v/v!BTe!N', YTITLE ='Electrons F(V)',BACKGROUND=255, COLOR=0

;WDELETE,0

;my_image =OBJ_NEW('IDLgrImage', image1)
;mypalette =OBJ_NEW('IDLgrPalette')
;mypalette-> LoadCT, 39
;my_image -> SetProperty, Palette=mypalette
;myMPEG -> Put, my_image, i

ps_file= 'FLS' +string(format='(I5.5)', i,/print)+'.ps'

oldDevice = !d.name
Set_plot, 'PS'
Device, /ENCAPSUL, Bits_per_Pixel =8, /COLOR,filename=ps_file, xsize=12,YSIZE=16

!p.charsize=2
!P.thick=1.2
!P.MULTI=[0,1,3]

plot,Vel, F_v/sqrt(3.1416), xrange=[-15.,15.], xticks=6,Thick=2,/ylog,yrange=[1e-12,1e-3], XTITLE ='Velocity, v/v!BTe!N', YTITLE ='Electrons F(V)'
oplot,sc_v,sc_fv,line=1


plot,kx , W_v/WL_thermal,/ylog,ystyle=1,xstyle=1,xrange=[-0.5,0.5],xticks=10,thick=2,yrange=[0.1,1e8],XTITLE ='Wavenumber, k!4k!3!BDe!N', YTITLE ='Langmuir W!BL!N(k)/W!BT!N'

plot,kx , W_S/WS_thermal,/ylog,ystyle=1,xstyle=1,xrange=[-0.5,0.5],xticks=10,thick=2,yrange=[0.1,1e3],XTITLE ='Wavenumber, k!4k!3!BDe!N', YTITLE ='Sound W!BS!N(k)/W!BT!N'

;XYOUTS ,240,582,COLOR=50, string(format='(A, F6.4, A)', 'Time=', float(i)/(1e+4), 's', /print ), /device

!p.charsize=1
!P.thick=1
!P.MULTI=0

DEVICE, /CLOSE
SET_PLOT, oldDevice


end
;myMPEG -> Save
;GIF file saved to disk
DEVICE, DECOMPOSED=orig_decomposed

For j=0, N_elements(ImageCube[0,0,*])-1 DO write_gif,'eLS.gif', imageCube[*,*,j],r, g, b,/multiple
write_gif,'eLS.gif', /close

oVid = 0
PRINT, 'image loaded & mpeg created ...................  ok !'

stop
end
