
!P.CHARSIZE=2
tvlct,[0,255,0,0],[0,0,255,0],[0,0,0,255]
device,decomposed=0

Nchan = 150

nu = dindgen(Nchan)*1.e6+50.e6
c = 2.997e8
lambda = c/nu

x = 1.8
deltx=x

;sinc_func = fltarr(Nchan)

;for i=0,Nchan-1 do sinc_func[i] = 2.*sinc(2.*!pi*x/lambda[i])

tb=150.*(tanh((nu/1.e6-80.)/20.)-0.95)

; short dipole model for primary beam

Nsize = 127

dipole_height = 0.29
make_azza_arrays_fov,Nsize,az,za,mask,180.

beam_dipole = dblarr(Nsize,Nsize,Nchan)

for ch = 0 , Nchan-1 do begin
print,nu[ch]

  field_x = complex(za*0.0)
  field_y = complex(za*0.0)
  lambda = 300.0e6/nu[ch]
	print,lambda
  p = where(za lt !pi/2)
  gp_norm = end_fire_2element(2*(dipole_height/lambda),0.0)
  tile_norm = 1. ; how many dipoles in tile
  zaaz_to_xyz,za[p],az[p],x,y,z

    ground_plane_effect = end_fire_2element(2*dipole_height/lambda,za[p])/gp_norm

  proj_x = sqrt(1-(sin(az[p])*sin(za[p]))^2)
  proj_y = sqrt(1-(cos(az[p])*sin(za[p]))^2)

field_x[p] = ground_plane_effect*proj_x
field_y[p] = ground_plane_effect*proj_y

beam_dipole[*,*,ch] = float((field_x^2+field_y^2)/max(field_x^2 + field_y^2))

endfor

l = dindgen(Nsize+2)/((Nsize+1.)/2.)-1.
m = l

padsize=45.

; pad beams to obtain delta_u = 0.1

Nsizebig = floor((Nsize+1.)*padsize+1.)

beam_dipole_pad = fltarr(Nsizebig,Nsizebig)

inp = ((padsize-1.)/2.)
inp1 = ((padsize-1.)/2.+1.)

outp = ((padsize-1.)/2.+2.)

delta_u = 1./(2.)/padsize
u = (dindgen(Nsizebig))*delta_u

pl = where(abs(u-6.) eq min(abs(u-6.)))

beam_fft_out = fltarr(pl+1,pl+1,Nchan)

for ch = 0 , Nchan-1 do begin
;for ch = 0 , 19 do begin
print,ch

beam_dipole_pad[inp*Nsize+outp:inp1*Nsize+outp-2,inp*Nsize+outp:inp1*Nsize+outp-2] = beam_dipole[1:Nsize-1,1:Nsize-1,ch]

; Fourier Transform to uv-plane

beam_fft = fft(beam_dipole_pad[*,*],/double)

beam_fft_out[*,*,ch] = real_part(beam_fft[0:pl,0:pl])

endfor

; find where fourier beam is sampled for a given frequency and attenuate mock signal

atten = fltarr(Nchan)

for ch = 0 , Nchan-1 do begin
	loc = where(abs(uout-deltx*nu[ch]/c) eq min(abs(uout-deltx*nu[ch]/c)))
	atten[ch] = abs(beam_fft_out[0,loc,ch])
endfor


;plot,nu/1.e6,tb*abs(atten),xtitle='Frequency (MHz)',color=0,background=-1,ytitle='Observed Temperature (K)',title='Visibility response to T!dB!N at 1.5m',yrange=[-220.,50.],ystyle=1                                          
;oplot,nu/1.e6,tb*abs(atten),color=1,thick=2
;oplot,nu/1.e6,tb,color=3,thick=2 

end
