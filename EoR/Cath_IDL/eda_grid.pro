tvlct,[0,255,0,0],[0,0,255,0],[0,0,0,255]
device,decomposed=0

!P.CHARSIZE=2

common wblock,ww

file=read_ascii('ant_coords_EDA_east_north_metres.txt')

ant = 256
Nuu = 541

xpos = file.field1[0,*]
ypos = file.field1[1,*]
numpos = indgen(ant)

xpos = xpos-mean(xpos)
ypos = ypos-mean(ypos)


Nbaselines = long(float(ant)*(float(ant)-1.)/2.)

; assume the phase centre is at the zenith

; u and v are in wavelengths
xx = dblarr(Nbaselines)
yy = dblarr(Nbaselines)
numbase = intarr(Nbaselines,2)

lambda = 1.

counter = 0L
for i=0,ant-1 do begin
for j=i,ant-1 do begin

if ((i ne j)) then begin
xx[counter] = (xpos[i]-xpos[j])
yy[counter] = (ypos[i]-ypos[j])
numbase[counter,0] = numpos[i]
numbase[counter,1] = numpos[j]
counter = counter + 1L
;print,counter
endif

endfor
endfor

leng = xx^2+yy^2

xx2 = xx[sort(leng)]
yy2 = yy[sort(leng)]
numbase2 = numbase[sort(leng),*]
leng2 = leng[sort(leng)]

xx=xx2
yy=yy2
numbase = numbase2

sep = sqrt(leng)

;; define frequencies

Nchan = 150
Nchangrid = 150

nu = dindgen(Nchan)*1.e6+50.e6
nu_m = nu/1.e6
c = 2.997e8
lambda = c/nu
kboltz = 1380.

;;;; IMPORT COMPONENTS

S_21 = complexarr(Nuu,Nuu,Nchan)
S_GS = S_21*0.
S_PS = S_21*0.
weights = fltarr(Nuu,Nuu,Nchan)

;; S_21

restore,'beam.sav'
restore,'beam_150channels.sav'
norm = 8102.81     ; this is the FFT normalisation = N^2 Delta_l^2 from 'beam_extract_lm.pro'

; signal model
tb=150.*(tanh((nu/1.e6-80.)/20.))/1000.    ; K  (150mK is the level here)
tb = 1.e-4*nu/1.e6 - 0.15*exp(-(nu/1.e6-70.)^2/2./10.^2)

for ch=0,Nchangrid-1 do S_21[*,*,ch] = complex(beam_fft_out[*,*,ch],0)*tb[ch]*2.*kboltz/lambda[ch]^2*norm/1.5    ; should be Jy
for ch=0,Nchangrid-1 do weights[*,*,ch] = beam_fft_out[*,*,ch]/1.5*norm

;; S_GS

for ch=0,Nchangrid-1 do begin

fileout = 'sky_fft_Ksr_000MHz.dat'

str = nu[ch]/1.e6
if (nu[ch] lt 100.e6) then strput,fileout,string(str,format='(i2)'),13
if (nu[ch] ge 100.e6) then strput,fileout,string(str,format='(i3)'),12

S_GS[*,*,ch] = read_binary('../../../ska/sky_map/gsm2016-master/'+fileout,data_type=6,data_dims=[Nuu,Nuu])*2.*kboltz/lambda[ch]^2   ; Jy

endfor

;; S_PS

restore,'../../../ska/sky_map/gsm2016-master/s_ps_150_grid.sav'
S_PS = vis
; already in Jy

;;;; COMPUTE NOISE CHARACTERISTICS

time_tot = 3600.*2.   ; seconds
bw = nu[1]-nu[0]

noise = fltarr(Nchan)

tsys = 20.*(408.e6/nu)^2.75 + 2.73 + 288.*(0.005+0.1314*exp((alog(nu/1.e9)-alog(22.23))*8.))
FOV = beam_fft_out[0,0,*]*norm
noise = 2.*kboltz*tsys/sqrt(bw*time_tot)*FOV/sqrt(2.)   ; real part

; ________________________________________________________

;; GRID VIS AND WEIGHTS AT EACH BASELINE FOR EACH CHANNEL

S_TOT = S_21 + S_GS; + S_PS

; Nyquist sample data

Nnew = floor((Nuu-1.)/2.)
unew = fltarr(Nnew)

S_21_new = complexarr(Nnew,Nnew,Nchan)
S_PS_new = complexarr(Nnew,Nnew,Nchan)
S_GS_new = complexarr(Nnew,Nnew,Nchan)
weights_new = fltarr(Nnew,Nnew,Nchan)

count=0
count2=0

for i=0,Nnew-1 do begin
	unew[i] = uout[2.*i]
	for j=0,Nnew-1 do begin
	S_21_new[i,j,*] = S_21[2*i,2*j,*]
	S_PS_new[i,j,*] = S_PS[2*i,2*j,*]
	S_GS_new[i,j,*] = S_GS[2*i,2*j,*]
	weights_new[i,j,*] = weights[2*i,2*j,*]
	endfor
endfor

d_unew = unew[1]-unew[0]

visout = complexarr(Nchan)
vis21out = visout*0.
visgsout = visout*0.
vispsout = visout*0.
visnoise = visout*0.
numch = intarr(Nchan)
weightsout = complexarr(Nchan)
numbaseout = intarr(2,Nbaselines)

Nbaselinesclose = 32640
offset = 0

for i = offset , Nbaselinesclose+offset-1 do begin
	for ch = 0 , Nchan-1 do begin

	uu = abs(xx[i]/lambda[ch]/d_unew)
	vv = abs(yy[i]/lambda[ch]/d_unew)

	vnoise = randomn(seed)*noise[ch]

	test_weight = bilinear(weights_new[*,*,ch],uu,vv)
	if (abs(test_weight) gt 0.5) then begin

    vis21out[ch] += bilinear(s_21_new[*,*,ch],uu,vv)

	visgsout[ch] += bilinear(s_gs_new[*,*,ch],uu,vv)
	    vispsout[ch] += bilinear(s_ps_new[*,*,ch],uu,vv)
        visnoise[ch] += vnoise
        numch[ch] += 1
	if (ch eq 0) then begin
	numbaseout[*,count] = numbase[i,*]
	count = count + 1L
	endif
	
	weightsout[ch] += test_weight

	endif

	endfor

endfor

ww = weightsout

diffgs = fltarr(Nchan-1)
diffps = fltarr(Nchan-1)
for i=0,Nchan-2 do begin
diffgs[i] = visgsout[i+1]-visgsout[i]
diffps[i] = vispsout[i+1]-vispsout[i]
endfor

; smooth out the model errors in the real part of the GS signal

yf = poly_fit(nu[0:Nchan-1],real_part(visgsout[0:Nchan-1]),3,yfit=realvisgsout)
yf = poly_fit(nu[0:Nchan-1],imaginary(visgsout[0:Nchan-1]),3,yfit=imagvisgsout)

visgsoutnew = complex(realvisgsout,imagvisgsout)
stop
visgsout = visgsoutnew

S_TOT_new = S_21_new + S_GS_new + S_PS_new

visout = visgsout + vis21out + vispsout + visnoise

power = visout   ;*conj(visout)
power21 = vis21out   ;*conj(vis21out)
powergs = visgsout   ;*conj(visgsout)
powerps = vispsout   ;*conj(vispsout)

;yf = poly_fit(nu[0:Nchan-1],power[0:Nchan-1],2,yfit=powerfit)
;yf = poly_fit(nu[0:Nchan-1],power21[0:Nchan-1],2,yfit=power21fit)
;yf = poly_fit(nu[0:Nchan-1],powergs[0:Nchan-1],2,yfit=powergsfit)


; perform the fitting routine. We want to fit a polynomial-based power-law to the GS and PS and a tanh-like function, weighted by the beam weights, to the 21cm signal.


;; [A, B, C, D, nu_c, delt_nu]
;; [Jy, unitless, unitless, Jy, MHz, MHz]

;;nu_0 = 50.   ; MHz

;;sig1 = a[0]*(nu/nu_0)^(a[1]+a[2]*alog(nu))*nu^2*ww
;;sig2 = a[3]*tanh((nu-a[4])/a[5])*nu^2*ww

A0 = [1.e9,-2.72661,0.01,-1.,-52.1967,70.,10.,3.92]
A0=[3.19540e+09 ,    -1.38843 ,    0.178591  , 0.00143414 ,    -6964.04 ,    -193.757 ,     76.1507 ,     0.00000]

res=curvefit(nu_m,real_part(visout),fltarr(Nchan)+1.,A0,fita=[1,1,1,1,1,1,1,1],function_name='signal_fg',itmax=100,chisq=chisq)


;oplot,[xpos[numbaseout[0,0]],xpos[numbaseout[0,0]]],[ypos[numbaseout[0,0]],ypos[numbaseout[0,0]]],color=1,psym=6,symsize=2

; ________________________________________________________

;; histogram stuff

;mx = max(sep)
;mn = 0.

;nbins = 352

;xhist = findgen(nbins)/float(nbins)*mx

;hist = histogram(sep,min=mn,max=mx,nbins=nbins)

;bb = barplot(xhist[10:35],hist[10:35],xtitle='Separation (m)',ytitle='Number of baselines')

;cumul = hist*0.
;for i=1,nbins-1 do cumul[i] = cumul[i-1]+hist[i]

;bb = barplot(xhist[10:35],hist[10:35],nbars=2,index=0,fill_color='green',xtitle='Separation (m)',ytitle='Number of baselines',yrange=[0.,1150.],/ylog)
;bc = barplot(xhist[10:35],cumul[10:35],nbars=2,index=1,fill_color='gold',/overplot)
;text_pine = TEXT(115, 300, 'Histogram', /CURRENT,COLOR='green', /DATA)
;text_yel = TEXT(115, 275, 'Cumulative', /CURRENT, COLOR = 'gold', /DATA)

end
