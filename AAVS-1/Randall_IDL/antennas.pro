################################
antennas.pro (IDL) from Randall
################################
; code for calculating the array antenna response for a square array
; of ideal (isotropic) antennas. Randall Wayth, May 2006 +


; project the power pattern from theta,phi coords into a plane
; cartesian coords using an azimuthal equal area projection.
; dat is the array of powers, assumed to be in theta,phi coords
; proj is the result in x,y cartesian.
pro project_azimuthal_equal_areas,dat,proj
  s = size(dat)
  proj = dat*0.
  scale = s[1]/2.0

  if (s[1] ne s[2]) then begin
    print,'ERROR: array must be square'
    return
  endif

  for x=0, s[1]-1 do begin
    for y=0, s[1]-1 do begin
      ; scale to be between -1 and 1
      xs = (x - scale)/scale
      ys = (y - scale)/scale
      r = sqrt(xs^2 + ys^2)

      if (r lt 1.0) then begin
        ; calculate the az and el of the coord
        az = atan(ys,xs)
        el = acos(r/sqrt(2.0))*2.0

        ; project x and y
        xp = cos(az)*sin(el)
        yp = sin(az)*sin(el)

        ; convert this to a theta, phi coords and scale to grid units
        th = acos(xp)
        ph = acos(yp)
        ths = (!pi-th)/!pi*s[1]
        phs = (!pi-ph)/!pi*s[1]

        ; sample nearest grid unit
        proj[x,y] = dat[round(ths),round(phs)]
      endif
    endfor
  endfor
end


; project the full sky from a plan-carre ra/dec coord system to a 
; Hammer-Aitoff all-sky projection
; for the x and y coords on a grid, the corresponding lon/lat can be calculated
; as lon = 2*atan2[2z^2-1,zx/2]
; lat = arcsin(zy), where
; z = sqrt((1+cos(lat)cos(lon/2))/2)
; or z = sqrt(1.0-x^2/8 - y^2/4)
; the extrema of lon at -2sqrt(2) and 2sqrt(2)
; for lat it is -sqrt(2) and sqrt(2)
function project_hammer_aitoff,sky_radec
  s = size(sky_radec)
  sx = s[1]
  sy = s[2]
  res = fltarr(sx,sy)
  ras = (findgen(sx)-sx/2)*2.0*sqrt(2)/(sx/2.0)
  decs= (findgen(sy)-sy/2)*sqrt(2)/(sy/2.0)
  ra_grid = res
  dec_grid= res
  for d=0,n_elements(decs)-1 do ra_grid[*,d]=ras
  for r=0,n_elements(ras)-1 do dec_grid[r,*]=decs
  ; make projections
  z = sqrt(1.0-(ra_grid*ra_grid/16.) - (dec_grid*dec_grid/4))
  lat = asin(z*dec_grid)
  lon = 2.0*atan(2.0*z*z-1.0,z*ra_grid/2.0)
  ;lon = 2.0*atan(z*ra_grid/(2.0*(2.0*z*z-1.0)))
  p = where(finite(z))
  res[p] = sky_radec[(lon[p]+!pi)*sx,(lat[p]+!pi/2)*sy]
  return,res
end

; make az,za arrays for a field of view (degrees)
; the ZEA option make the image in Zenithal Equal Area projection
pro make_azza_arrays_fov,gsize,az,za,mask,fov,ZEA=zea
  if fov gt 180.0 then message,'FOV is too large. max: 180 deg'

  c = (indgen(gsize,gsize,/long) mod (gsize) - gsize/2.0)/(gsize/2)
  r = (indgen(gsize,gsize,/long) / (gsize) - gsize/2.0)/(gsize/2)
  myfov = sin(fov/2*!pi/180.)
  mask = fltarr(gsize,gsize)
  dsqu = (c^2+r^2)*(myfov)^2
  p = where(dsqu lt 1.0)
  za = fltarr(gsize,gsize)
  if Keyword_Set(ZEA) eq 1 then begin
    ; for ZEA, R = sin((FOV/2)/2)
    r_fac = sin(!pi/180*fov/2/2)
    za[p] = 2*asin(sqrt(dsqu[p])/2.0/r_fac)
  endif else begin
    za[p] = asin(sqrt(dsqu[p]))
  endelse

  az = atan(c,r)
  mask[p] =1.0
  p = where(dsqu ge 1.0)
  if (p[0] ne -1) then za[p] = !pi/2.0;

end


; calculate the E field response of a 2-element end-fire antenna with
; isotropic receivers separated by d, and opposite phase
; which is the same as a single element with 
; a groundplane and distance d/2 above ground.
; separation of antennas: d (wavelengths)
; angle of incoming radiation = th (radian, from line joining
; antennas)
; see Krauss equation 4-10, section 4-2.
function end_fire_2element,d,th
  return,sin(!pi*d*cos(th))*2.0
end


; calculate the E field of a finite length dipole
; l = length of dipole (wavelengths), t = angle from dipole axis (== elevation for horizontal dipole)
; this function can handle only 1 length param but arrays of angle are supported
function dipole_response,l,t
    res = t*0.0
    p = where(sin(t) ne 0.0) ; handle divide by zero
    if p[0] ne -1 then res[p] = (cos(!pi*l*cos(t[p])) - cos(!pi*l))/sin(t[p])
    return, res
end


; write an (assumed) all-sky beam image to a FITS file
; with a zenith orthographic (sine) coord system.
pro write_beam_fits,img,filename,point_az,point_za,freq_MHz,lst_hours
    mkhdr,hdr,float(img)
    dec=-26.7033
    siz = size(img)
    pixscale=180.0/(siz[1]/2) ; for all-sky
    sxaddpar,hdr,"CRPIX1",siz[1]/2+1
    sxaddpar,hdr,"CDELT1",pixscale/!pi
    sxaddpar,hdr,"CRVAL1",lst_hours*15.0
    sxaddpar,hdr,"CTYPE1","RA---SIN"
    sxaddpar,hdr,"CRPIX2",siz[2]/2+1
    sxaddpar,hdr,"CDELT2",pixscale/!pi
    sxaddpar,hdr,"CRVAL2",dec
    sxaddpar,hdr,"CTYPE2","DEC--SIN"    
    sxaddpar,hdr,"BEAM_AZ",point_az,' [deg]'
    sxaddpar,hdr,"BEAM_ZA",point_za,' [deg]'
    sxaddpar,hdr,"FREQ",freq_MHz*1e6,' [Hz]'
    writefits,filename,float(img),hdr
end


;============================================================
; more generic code for antennas with imperfect dipoles.
; Dipoles can be twisted or have gains that differ to 1.

; EXAMPLE:
;make_azza_arrays_fov,256,az,za,mask,180.0
;tiledata = setup_dipoles(0.1,0,0.00)
;tiledata = setup_dipoles(0.0,1,0.0) ; ideal dipoles.
;point_za=0.0  ; degs
;point_az=0.0  ; degs
;freq=150.0    ; MHz
;setup_tile_pointing_direction,point_za*!pi/180,point_az*!pi/180,tiledata
;calc_directional_response,za,az,freq,tiledata,fieldx,fieldy
;tvscl,sqrt((abs(fieldx)/max(abs(fieldx)))^2+1e-3)
;


; setup for MWA dipoles. gainerror is rms error relative to
; 1.0 angle_error is rms angle error in deg, height error is rms
; height offset in m
function setup_dipoles,gainerror,angle_error,height_error,NDIPOLES=ndipoles
  n_dipoles=4
  if keyword_set(NDIPOLES) then n_dipoles=ndipoles
  relative_gain_error = gainerror
  dipole_separation = 1.1 ; meters.
  dipole_height = 0.29 ; meters
  ; params for dipoles: x_pos,y_pos: x,y position of dipole relative to centre
  ; theta_x,_theta_y: angle offset of x,y feeds relative to E, N.
  ; gain_x,gain_y: gain of x,y feed (1=normal)
  ; gain_xy: crosstalk between x and y feeds
  ; phase: delay for tile relative to fixed point on tile
  ; height error
  dat = fltarr(n_dipoles,n_dipoles,9)
  offset_x = dat[*,*,0]
  offset_y = dat[*,*,0]
  theta_x =  dat[*,*,0]
  theta_y =  dat[*,*,0]
  gain_x =   fltarr(n_dipoles,n_dipoles)+1.0
  gain_y =   fltarr(n_dipoles,n_dipoles)+1.0
  gain_xy =  dat[*,*,0]

  for i=0,n_dipoles-1 do begin
    dat[i,*,0] = (i-n_dipoles/2)*dipole_separation + dipole_separation/2    ; X position
    dat[*,i,1] = (i-n_dipoles/2)*dipole_separation + dipole_separation/2    ; Y position
    ; make phase centre at corner of tile not middle
    ;dat[i,*,0] = (i)*dipole_separation
    ;dat[*,i,1] = (i)*dipole_separation
  endfor

  dat[*,*,0] += offset_x    ; dipole offsets
  dat[*,*,1] += offset_y
  dat[*,*,2] = theta_x      ; dipole rotations (not fully supported)
  dat[*,*,3] = theta_y
  dat[*,*,4] = gain_x + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error   ; gains
  dat[*,*,5] = gain_y + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error
  dat[*,*,6] = gain_xy      ; crosstalk
  dat[*,*,8] = dipole_height*(1.0 + randomn(seed,n_dipoles,n_dipoles)*height_error) ; height in meters

  return, dat
end


; setup model for AAVS0.5 prototype station
function setup_dipoles_aavs05,freq_MHz
  n_dipoles=4
  relative_gain_error = 0.0
  height_error = 0.0
  dipole_height =  300./freq_MHz/4. ; meters
  ; params for dipoles: x_pos,y_pos: x,y position of dipole relative to centre
  ; theta_x,_theta_y: angle offset of x,y feeds relative to E, N.
  ; gain_x,gain_y: gain of x,y feed (1=normal)
  ; gain_xy: crosstalk between x and y feeds
  ; phase: delay for tile relative to fixed point on tile
  ; height error
  dat = fltarr(n_dipoles,n_dipoles,9)

  dat[*,*,0] = [[-0.3070,-2.3760,1.6641,-0.0841],[-3.8039,3.0916,-1.9308,0.6160],[3.7586,-3.5161,-0.6200,1.0324],[2.6793,-2.9657, 0.4716,-1.6010]]
  dat[*,*,1] = [[-3.8874,-2.5737,-2.4956,-2.0228],[-1.2085,-1.1204,-0.6729,-0.2155],[0.5343,0.5463,1.2723,1.5772],[2.0953,2.2106,3.2399,3.5317]]

  theta_x =  [ [0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0] ]
  theta_y =  [ [0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0] ]
  gain_x =   [ [1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0] ]
  gain_y =   [ [1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0] ]
  gain_xy =  [ [0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0] ]

  dat[*,*,2] = theta_x      ; dipole rotations (not fully supported)
  dat[*,*,3] = theta_y
  dat[*,*,4] = gain_x + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error   ; gains
  dat[*,*,5] = gain_y + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error
  dat[*,*,6] = gain_xy      ; crosstalk
  dat[*,*,8] = dipole_height*(1.0 + randomn(seed,n_dipoles,n_dipoles)*height_error) ; height in meters

  return, dat
end



function setup_lfaa,station_diam,freq_mhz,OFFSET_ERROR=offseterror
  offset_error = 0.0
  if keyword_set(OFFSET_ERROR) then offset_error=offset_error
  if offset_error eq 0.0 then offset_error = 0.01
  dipole_separation = 1.9 ; meters.
  n_ant = 256L 	; target number of antennas
  n_dipoles=2*fix(station_diam/2.0/dipole_separation)
  print,'Estimated number of elements: '+string(n_dipoles*n_dipoles*0.785398)
  relative_gain_error = 0.0
  height_error = 0.0
  dipole_height = 300./freq_mhz/4. ; meters
  ; params for dipoles: x_pos,y_pos: x,y position of dipole relative to centre
  ; theta_x,_theta_y: angle offset of x,y feeds relative to E, N.
  ; gain_x,gain_y: gain of x,y feed (1=normal)
  ; gain_xy: crosstalk between x and y feeds
  ; phase: delay for tile relative to fixed point on tile
  ; height error
  dat = fltarr(n_dipoles,n_dipoles,9)
  offset_x = randomu(seed,n_dipoles,n_dipoles)*offset_error-offset_error/2.0
  offset_y = randomu(seed,n_dipoles,n_dipoles)*offset_error-offset_error/2.0
  theta_x =  reform(dat(*,*,0))*0.0
  theta_y =  reform(dat(*,*,0))*0.0
  gain_x =   reform(dat(*,*,0))*0.0 + 1.0
  gain_y =   gain_x
  gain_xy =  reform(dat(*,*,0))*0.0

  for i=0,n_dipoles-1 do begin
    dat[i,*,0] = (i-n_dipoles/2)*dipole_separation + dipole_separation/2
    dat[*,i,1] = (i-n_dipoles/2)*dipole_separation + dipole_separation/2
    ; make phase centre at edge of tile not middle
    ;dat[i,*,0] = (i)*dipole_separation
    ;dat[*,i,1] = (i)*dipole_separation
  endfor

  dat[*,*,0] += offset_x    ; dipole offsets
  dat[*,*,1] += offset_y
  dat[*,*,2] = theta_x      ; dipole rotations (not fully supported)
  dat[*,*,3] = theta_y
  dat[*,*,4] = gain_x + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error   ; gains
  dat[*,*,5] = gain_y + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error
  dat[*,*,6] = gain_xy      ; crosstalk
  dat[*,*,8] = dipole_height*(1.0 + randomn(seed,n_dipoles,n_dipoles)*height_error) ; height in meters

  dist = sqrt(dat[*,*,0]^2 + dat[*,*,1]^2)
  p = where(dist gt (station_diam/2.0))
  gainmask = reform(dat[*,*,4])*0.0 + 1.0
  gainmask[p] = 0.0
  dat[*,*,4] *= gainmask
  dat[*,*,5] *= gainmask

  n_good =  fix(total(gainmask))

  print,'There are '+string(n_good)+' non-zero antennas'
  if n_good gt n_ant then begin
  endif

  return, dat
end

function setup_Ph3_Hex96,OFFSET_ERROR=offset_error,STWIDTH=stwidth
  dip_offset_error = 0.0
  if keyword_set(OFFSET_ERROR) then dip_offset_error =offset_error
  dips = read_ascii(GETENV("HOME")+"/Dropbox/MWA/Phase3/STATIONS/hex96_regular.txt",COMMENT_SYMBOL="#")
  n_dipoles = n_elements(dips.field1)/2
  relative_gain_error = 0.0
  height_error = 0.0
  dipole_height = 0.29 ; meters
  width = 14.0 ; meters
  if keyword_set(STWIDTH) then width=stwidth
  ; params for dipoles: x_pos,y_pos: x,y position of dipole relative to centre
  ; theta_x,_theta_y: angle offset of x,y feeds relative to E, N.
  ; gain_x,gain_y: gain of x,y feed (1=normal)
  ; gain_xy: crosstalk between x and y feeds
  ; phase: delay for tile relative to fixed point on tile
  ; height error
  dat = fltarr(n_dipoles/16,16,9)
  offset_x = randomu(seed,n_dipoles/16,16)*dip_offset_error -dip_offset_error/2.0
  offset_y = randomu(seed,n_dipoles/16,16)*dip_offset_error -dip_offset_error/2.0
  theta_x =  reform(dat(*,*,0))*0.0
  theta_y =  reform(dat(*,*,0))*0.0
  gain_x =   reform(dat(*,*,0))*0.0 + 1.0
  gain_y =   gain_x
  gain_xy =  reform(dat(*,*,0))*0.0

  dat[*,*,0] = reform(dips.field1[0,*],n_dipoles/16,16) * (width/2)
  dat[*,*,1] = reform(dips.field1[1,*],n_dipoles/16,16) * (width/2)

  dat[*,*,0] += offset_x    ; dipole offsets
  dat[*,*,1] += offset_y
  dat[*,*,2] = theta_x      ; dipole rotations (not fully supported)
  dat[*,*,3] = theta_y
  dat[*,*,4] = gain_x + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error   ; gains
  dat[*,*,5] = gain_y + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error
  dat[*,*,6] = gain_xy      ; crosstalk
  dat[*,*,8] = dipole_height*(1.0 + randomn(seed,n_dipoles,n_dipoles)*height_error) ; height in meters

  gainmask = reform(dat[*,*,4])*0.0 + 1.0
;  gainmask[p] = 0.0
;  dat[*,*,4] *= gainmask
;  dat[*,*,5] *= gainmask

  return, dat
end

function setup_Ph3_Hex60,OFFSET_ERROR=offset_error,STWIDTH=stwidth
  dip_offset_error = 0.0
  if keyword_set(OFFSET_ERROR) then dip_offset_error =offset_error
  dips = read_ascii(GETENV("HOME")+"/Dropbox/MWA/Phase3/STATIONS/hex60_regular.txt",COMMENT_SYMBOL="#")
  n_dipoles = n_elements(dips.field1)/2
  relative_gain_error = 0.0
  height_error = 0.0
  dipole_height = 0.29 ; meters
  width = 10.0 ; meters
  if keyword_set(STWIDTH) then width=stwidth
  ; params for dipoles: x_pos,y_pos: x,y position of dipole relative to centre
  ; theta_x,_theta_y: angle offset of x,y feeds relative to E, N.
  ; gain_x,gain_y: gain of x,y feed (1=normal)
  ; gain_xy: crosstalk between x and y feeds
  ; phase: delay for tile relative to fixed point on tile
  ; height error
  dat = fltarr(n_dipoles,1,9)
  offset_x = randomu(seed,n_dipoles,1)*dip_offset_error -dip_offset_error/2.0
  offset_y = randomu(seed,n_dipoles,1)*dip_offset_error -dip_offset_error/2.0
  theta_x =  reform(dat(*,*,0))*0.0
  theta_y =  reform(dat(*,*,0))*0.0
  gain_x =   reform(dat(*,*,0))*0.0 + 1.0
  gain_y =   gain_x
  gain_xy =  reform(dat(*,*,0))*0.0

  dat[*,*,0] = reform(dips.field1[0,*],n_dipoles,1) * (width/2)
  dat[*,*,1] = reform(dips.field1[1,*],n_dipoles,1) * (width/2)

  dat[*,*,0] += offset_x    ; dipole offsets
  dat[*,*,1] += offset_y
  dat[*,*,2] = theta_x      ; dipole rotations (not fully supported)
  dat[*,*,3] = theta_y
  dat[*,*,4] = gain_x + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error   ; gains
  dat[*,*,5] = gain_y + randomn(seed,n_dipoles,n_dipoles)*relative_gain_error
  dat[*,*,6] = gain_xy      ; crosstalk
  dat[*,*,8] = dipole_height*(1.0 + randomn(seed,n_dipoles,n_dipoles)*height_error) ; height in meters

  gainmask = reform(dat[*,*,4])*0.0 + 1.0
;  gainmask[p] = 0.0
;  dat[*,*,4] *= gainmask
;  dat[*,*,5] *= gainmask

  return, dat
end


; setup properties of the dipoles. gainerror is rms error relative to
; 1.0 angle_error is rms angle error in deg, height error is rms
; height offset in m
function setup_dipoles_readlayout,filename
  
  t = read_ascii(filename,comment_symbol='#')
  input = t.field1
  n_dipoles=4

  ; params for dipoles: x_pos,y_pos: x,y position of dipole relative to centre
  ; theta_x,_theta_y: angle offset of x,y feeds relative to E, N.
  ; gain_x,gain_y: gain of x,y feed (1=normal)
  ; gain_xy: crosstalk between x and y feeds
  ; delay: delay (in meters) for tile relative to fixed point on tile
  ; height error
  dat = fltarr(n_dipoles,n_dipoles,9)
  offset_x = [ [0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0] ]
  offset_y = [ [0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0] ]
  theta_x =  [ [0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0] ]
  theta_y =  [ [0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0] ]
  gain_x =   [ [1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0] ]
  gain_y =   [ [1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0] ]
  gain_xy =  [ [0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0] ]

  for n_ind=0,n_dipoles-1 do begin
    for e_ind=0,n_dipoles-1 do begin
        inp_index = e_ind+n_ind*n_dipoles
        dat[e_ind,n_ind,0] = input[1,inp_index]
        dat[e_ind,n_ind,1] = input[0,inp_index]
        dat[e_ind,n_ind,4] = input[4,inp_index]
        dat[e_ind,n_ind,8] = input[7,inp_index]
    endfor
  endfor

  dat[*,*,0] += offset_x
  dat[*,*,1] += offset_y
  dat[*,*,2] = theta_x
  dat[*,*,3] = theta_y
  dat[*,*,6] = gain_xy

  return, dat
end


; adjust the phases of the tile (setup already) to point to a direction in the sky.
; inputs: za,az: zenith angle, azimuth (radians)
;         tiledata is an array from setup_dipoles()
;         z optional keyword to supply impedance scattering matrix
;         freq optional keyword (required for z) to turn delays into phases
pro setup_tile_pointing_direction,za,az,tiledata,IMPEDANCE=z,FREQ=freq
    delay_x = sin(az)*sin(za)
    delay_y = cos(az)*sin(za)
    tiledata[*,*,7] = delay_x*tiledata[*,*,0] + delay_y*tiledata[*,*,1]

    ; if an impedance matrix and freq have been supplied, then solve the port
    ; currents matrix to adjust the gains of each dipole depending on the look direction
    if Keyword_Set(IMPEDANCE) then begin
      s = size(tiledata)
      n_dipoles=s[1]
      lambda = 299.8/freq
      phases = tiledata[*,*,7]*(-2.0)*!pi/lambda
      complex_voltages = complex(cos(phases),sin(phases))
      ; unpack the voltages to the same ordering as the Z matrix
      v = complexarr(n_elements(phases))
      for i=0,n_dipoles-1 do v[i*n_dipoles:(i+1)*n_dipoles-1]=complex_voltages[*,n_dipoles-1-i]
      iz = invert(z)
      gains=iz##[v,v]   ; duplicate voltage matrix for Y and X dipoles
      ; unpack gains into tiledata
      gainx=gains[16:31]/max(abs(gains[16:31]))
      gainy=gains[0:15]/max(abs(gains[0:15]))
;      for i=0,n_dipoles-1 do tiledata[*,n_dipoles-1-i,4]
    end
end


; function to convert zenith angle (za), azimuth(az) coords in radians
; to xyz coords on a unit sphere. 
pro zaaz_to_xyz,za,az,x,y,z
  x = sin(az)*sin(za)
  y = cos(az)*sin(za)
  z = cos(za)
  return
end


; calculate the array factor field pattern of an array of non equal isotropic antennas ("dipoles").
; requires "tiledata" to be set up using the "setup_dipoles"
; procedure.
; za and az are arrays, of the same size, of points of zenith angle
; and azimuth in radian. freq is the frequency of observation (MHz).
pro calc_directional_response,za,az,freq,tiledata,field_x,field_y
  field_x = complex(za*0.0)
  field_y = complex(za*0.0)
  lambda = 300.0/freq
  p = where(za lt !pi/2)
  gp_norm = end_fire_2element(2*(tiledata[*,*,8]/lambda),0.0)
  s = size(tiledata)
  tile_norm = s[1]*s[2] ; how many dipoles in tile
  zaaz_to_xyz,za[p],az[p],x,y,z

  for i=0L,n_elements(p)-1 do begin

    ; due to ground plane, each dipole acts as a 2-element out of phase end-fire.
    ; this reduces sensitivity towards the horizon. the effect will be
    ; different for each dipole if they are not all the same height.
    ground_plane_effect = end_fire_2element(2*(tiledata[*,*,8]/lambda),za[p[i]])/gp_norm
    ;print,"GP effect: ",ground_plane_effect
    ;ground_plane_effect = 1.0

    ; phases are dipole position offsets dotted with sky position
    ; normalised by wavelength. Pointing direction phases must
    ; also be added.
    phases = -((x[i]*tiledata[*,*,0]) + (y[i]*tiledata[*,*,1])-tiledata[*,*,7])*(2.0*!pi)/(lambda)
    ;print,"phases: ",phases

    ; multiply dipole phase by gain and sum over tile
    field_x[p[i]] = total(complex(cos(phases),sin(phases))*tiledata[*,*,4]*ground_plane_effect)/tile_norm
    field_y[p[i]] = total(complex(cos(phases),sin(phases))*tiledata[*,*,5]*ground_plane_effect)/tile_norm
  endfor

  ; reduction due to geometrical projection of dipole onto sky
  ; for a NS dipole, the angle from the axis of the dipole, alpha
  ; is cos(alpha) = sin(za)cos(az)
  ; for EW cos(alpha) = sin(za)sin(za)
  ; fraction of dipole seen from sky pos is then sin(alpha) = sqrt(1-cos(alpha)^2)
  proj_x = sqrt(1-(sin(az[p])*sin(za[p]))^2)
  proj_y = sqrt(1-(cos(az[p])*sin(za[p]))^2)
  field_y[p] *= proj_y
  field_x[p] *= proj_x

  return
end

; make a bunch of GIFs that can be turned into a movie using the
; command: convert /tmp/img*.gif ant_lobes.mpeg
; (requires the mpeg2encode program for convert)
; gsize is the number of pixels along the edge of the image
pro make_freq_dep_movie,gsize
  window,xsize=gsize,ysize=gsize+50
  make_azza_arrays_fov,gsize,az,za,mask,180.0
  tiledata = setup_dipoles(0.0,1,0.0) ; ideal dipoles.
  loadct,39
  n_frames=250;
  freqs = findgen(n_frames)/n_frames*(330.-80.)+80.0
  dummy = fltarr(gsize,gsize)
  for i=0,n_frames-1 do begin
    contour,dummy,title=strtrim(freqs[i],2)+' MHz',xmargin=[0,0],ymargin=[0,2.5],/nodata
    
    ;  setup_tile_pointing_direction,45*!pi/180,20*!pi/180,200,tiledata ; point off zenith
    setup_tile_pointing_direction,0.0*!pi/180,0.0*!pi/180,tiledata ; point at zenith
    CALC_DIRECTIONAL_RESPONSE,za,az,freqs[i],tiledata,fieldx,fieldy
    tvscl,alog10((abs(field)/max(abs(field)))^2+1e-3)

    filename = '/tmp/img'+strtrim(i+100,2)+'.gif'
    write_gif,filename,tvrd()
  endfor
end



; take a rectangular image of the sky in RA/DEC coords, rotate it and project it
; into the sky location defined by the az and za arrays for a location
; defined by the latitude (degrees) and LST (hours)
; e.g:
; make_azza_arrays_fov,200,az,za,mask,180.0
; sky = mrdfits('/data2/BIGHORNS/radio408MHz_fixheader.fits')
; tvscl,map_sky(sky,0,-26.7033,az,za)
; h = findgen(24)
; for i=0,23 do tvscl,map_sky(sky,h[i],-26,az,za)
function map_sky,sky,lst_hr,lat,az,za

  s = size(sky)
  size_ra = s[1]
  size_dec = s[2]

  res = az*0.0

  ; use ALTAZ2HADEC in idl astro library
  ; takes AZ and ALT is degrees. resulting HA
  ; and DEC are also in degrees
  altaz2hadec,(!pi/2 - za)*180/!pi,az*180/!pi,lat,ha,dec

  ; the following assumes RA=0 in centre
  ; of the sky image and increases to the left.
  ra = (lst_hr - ha/15.0) mod 24.0
  ra_index = (((36-ra) mod 24)/24)*size_ra
  dec_index = (dec/180.0+0.5)*size_dec

  p = where(za lt !pi/2.0)

  res[p]=sky[ra_index[p],dec_index[p]]
;  p1 = where(za ge !pi/2)
;  res[p1] = min(res[p])

  return,res
  
end


pro make_sky_movie
  siz=400
  make_azza_arrays_fov,siz,az,za,mask,180.0
  window,xsize=siz,ysize=siz
  sky = readfits(getenv('HOME')+'/Dropbox/MWA/Images/radio408MHz_fixheader.fits')
  ;sky = sky*(-1.46579605785) + 24125
  tvscl,sqrt(map_sky(sky,0,-26.5,az,za))
  n_frames=48
  h = (findgen(n_frames)/n_frames*24.0 + 0) mod 24
  smax = max(sky)
  smin = min(sky)
;  loadct,5
  for i=0, n_frames-1 do begin
     s=map_sky(sky,h[i],-26.5,az,za)
     tvscl,sqrt((s-smin)/smax)
     filename = '/tmp/sky_ha'+strtrim(i+100,2)+'.png'
     write_png,filename,tvrd()
  endfor
end


; make a movie what a survey would see.
; convert resulting images using:
; mencoder "mf:///tmp/survey*.png" -mf fps=10 -ovc lavc -lavcopts vhq:vbitrate=1500 -o /tmp/survey_movie.avi
; mencoder "mf:///tmp/survey*.png" -mf fps=10 -ovc lavc -lavcopts vcodec=libx264 -o /tmp/survey_movie.mp4
pro make_survey_movie
  siz=720/2
;  siz=1280/2
  ;make_azza_arrays_fov,siz,az,za,mask,180.0,/ZEA   ; use Zenithal equal area projection
  make_azza_arrays_fov,siz,az,za,mask,180.0         ; use standard orthographic (sine) projection.
  DEVICE, DECOMPOSED=0  ; force color table to work for true color display
  loadct,4
  window,xsize=siz*2,ysize=480
  piccy = fltarr(siz*2,siz+100)
;  sky = readfits('/home/rwayth/Ubuntu One/MWA/Images/radio408MHz_fixheader.fits')
  sky = readfits('/tmp/radio408MHz_fixheader.fits')
  ;sky = readfits('/data2/BIGHORNS/radio408MHz_fixheader.fits')
  ; proposed scans 1 minute, 3 freqs per scan
  n_frames=1440/3
  h = (findgen(n_frames)/n_frames*24.0 + 0.0) mod 24
  freqs = [120.0,150.0,180.0]
  freqs = [100.0]
  decs = [ 18.6, 1.6, -13.0, -26.7, -40.2, -55.0, -72.0]
  decs = [ -17.0 ]
  tiledata = setup_dipoles(0.0,1,0.0) ; ideal dipoles.
  smax = max(sky)
  smin = min(sky)
  max_pix_val = 1000.0

  ; pre-compute the response of the tiles for all freqs and pointings
  tilebeam = fltarr(n_elements(freqs),n_elements(decs),siz,siz)
  for dec=0, n_elements(decs)-1 do begin
    for freq=0, n_elements(freqs)-1 do begin
        point_az = 0.0
        if decs[dec] lt -26.7 then point_az = !pi
        point_za = abs(decs[dec] + 26.7)*!pi/180.0
        setup_tile_pointing_direction,point_za,point_az,tiledata
        CALC_DIRECTIONAL_RESPONSE,za,az,freqs[freq],tiledata,fieldx,fieldy
        tilebeam[freq,dec,*,*] = fieldx^2
    endfor
  endfor

  for i=0L, n_frames-1 do begin
     s=reverse(map_sky(sky,h[i],-26.7,az,za))
     s_clip = s
     p = where(s gt max_pix_val)
     if p[0] ne -1 then s_clip[p] = max_pix_val
     smax_img = max_pix_val
     piccy[0:siz-1,0:siz-1] = sqrt((s_clip-smin)/smax_img)*255.
     dec_ind = i mod n_elements(decs)
     point_az = 0.0
     if decs[dec_ind] lt -26.7 then point_az = !pi
     point_za = abs(decs[dec_ind] + 26.7)*!pi/180.0
     for freq=0, n_elements(freqs)-1 do begin
       beamsky = tilebeam[freq,i mod n_elements(decs),*,*]*s
       p = where (beamsky gt max_pix_val)
       if p[0] ne -1 then beamsky[p] = maxbeamsky
       maxbeamsky = max(beamsky)
       piccy[siz:*,0:siz-1] = sqrt(beamsky/maxbeamsky)*255.0
       tv,piccy
       lst_text = string(FORMAT='(%"%05.2f")', h[i])
       text = 'Freq: '+strtrim(fix(freqs[freq]),2)+' MHz. LST: '+lst_text
       xyouts,siz,siz+25,text,alignment=0.5,/device
       ;return
       filename = '/tmp/survey_lst'+lst_text+'_freq_'+strtrim(fix(freqs[freq]),2)+'.png'
       write_png,filename,tvrd(/true)
     endfor
  endfor
end



pro make_sky_movie_steve
  wsiz=1024
  window,xsize=wsiz/2-50,ysize=wsiz/2
  files = ['movie_087.770_00.fits','movie_096.541_00.fits','movie_106.188_00.fits','movie_116.799_00.fits','movie_128.470_00.fits' ,'movie_141.307_00.fits' ,'movie_155.427_00.fits' ,'movie_170.959_00.fits' ,'movie_188.042_00.fits' ,'movie_206.832_00.fits' ,'movie_227.500_00.fits' ,'movie_250.233_00.fits' ,'movie_275.238_00.fits' ,'movie_302.741_00.fits' ,'movie_332.993_00.fits' ,'movie_366.267_00.fits' ,'movie_402.867_00.fits' ,'movie_443.124_00.fits' ,'movie_487.403_00.fits' ,'movie_536.107_00.fits' ,'movie_589.678_00.fits' ,'movie_648.602_00.fits' ,'movie_713.414_00.fits' ,'movie_784.702_00.fits' ,'movie_863.114_00.fits' ,'movie_949.362_00.fits' ,'movie_1044.227_00.fits' ,'movie_1148.573_00.fits' ,'movie_1263.345_00.fits' ,'movie_1389.585_00.fits' ,'movie_1528.441_00.fits' ,'movie_1681.171_00.fits' ,'movie_1849.163_00.fits' ,'movie_2033.942_00.fits' ,'movie_2237.185_00.fits' ,'movie_2460.738_00.fits' ,'movie_2706.629_00.fits' ,'movie_2977.091_00.fits' ,'movie_3274.579_00.fits' ,'movie_3601.793_00.fits' ,'movie_3961.705_00.fits' ,'movie_4357.582_00.fits' ,'movie_4793.016_00.fits' ,'movie_5271.962_00.fits' ,'movie_5798.766_00.fits' ,'movie_6378.212_00.fits' ,'movie_7015.560_00.fits' ,'movie_7716.595_00.fits' ,'movie_8487.681_00.fits' ,'movie_9335.819_00.fits' ,'movie_10268.707_00.fits' ,'movie_11294.816_00.fits' ,'movie_12423.458_00.fits' ,'movie_13664.881_00.fits' ,'movie_15030.355_00.fits' ,'movie_16532.274_00.fits' ,'movie_18184.274_00.fits' ,'movie_20001.350_00.fits' ,'movie_22000.000_00.fits' ]
  freqs = [087.770 ,096.541 ,106.188 ,116.799 ,128.470 ,141.307 ,155.427 ,170.959 ,188.042 ,206.832 ,227.500 ,250.233 ,275.238 ,302.741 ,332.993 ,366.267 ,402.867 ,443.124 ,487.403 ,536.107 ,589.678 ,648.602 ,713.414 ,784.702 ,863.114 ,949.362 ,1044.22 ,1148.57 ,1263.34 ,1389.58 ,1528.44 ,1681.17 ,1849.16 ,2033.94 ,2237.18 ,2460.73 ,2706.62 ,2977.09 ,3274.57 ,3601.79 ,3961.70 ,4357.58 ,4793.01 ,5271.96 ,5798.76 ,6378.21 ,7015.56 ,7716.59 ,8487.68 ,9335.81 ,10268.7 ,11294.8 ,12423.4 ,13664.8 ,15030.3 ,16532.2 ,18184.2 ,20001.3 ,22000.0]
  siz = size(files)
  nfiles=siz[1]
  mins = fltarr(nfiles)
  maxs = fltarr(nfiles)
  mask = circmask(wsiz,wsiz,wsiz/2,wsiz/2,wsiz/2/(!pi/2),/aperture)
  sm_mask = rebin(mask,wsiz/2,wsiz/2)
  p = where(mask gt 0.)
  sm_p = where(sm_mask gt 0.)
  ; first pass through- get min/max
  for i=11, nfiles-1 do begin
    sky = mrdfits(files[i])+2.7
    mins[i] = min(sky[p])
    maxs[i] = max(sky)
  endfor
  smax = max(maxs)
  smin = min(mins)
  print,maxs,mins
  print,smax,smin
  dummy = fltarr(wsiz/2,wsiz/2)
  loadct,5

  ; second pass through- make pics
  for i=11, nfiles-1 do begin
    sky = mrdfits(files[i])+2.7
    s = rebin(sky,wsiz/2,wsiz/2)
    rescale_s = s*.0
    ;rescale_s[sm_p] = alog10((s[sm_p])/smin)/alog10(smax/smin)
    rescale_s[sm_p] = alog10((s[sm_p])/mins[i])/alog10(maxs[i]/mins[i])
    print,min(rescale_s),max(rescale_s)
	;contour,dummy,title=string(fix(freqs[i]),format='(%"Frequency %d MHz")'),xmargin=[0,0],ymargin=[0,2.5],/nodata
    tv,rescale_s[25:wsiz/2-26,25:wsiz/2-26]*255
    ;colorbar,divisions=5,min=alog10(smin),max=alog10(smax),position=[0.10, 0.85, 0.90, 0.90],title="Log_10 sky temp (K)",charsize=1.5,color=250
    colorbar,divisions=5,min=alog10(mins[i]),max=alog10(maxs[i]),position=[0.10, 0.85, 0.90, 0.90],title="Log_10 sky temp (K)",charsize=1.5,color=250
	filename = string(fix(freqs[i]),format='(%"sky_%05d.gif")')
    print,'making sky image for i:',i,',file ',filename
    write_gif,filename,tvrd()
  endfor
  
end

; make some plots of the angle between "feed" and sky l,m:
; make_azza_arrays_fov,200,az,za,mask,180.0
; lat = -26.
; altaz2hadec,(!pi/2 - za)*180/!pi,az*180/!pi,lat,ha,dec
; ha = (ha+180.)mod 360. - 180.  ; convert to ha=0 at zenith
; ha = ha*!pi/180.0
; dec = dec*!pi/180
; phi_q = atan(cos(ha),sin(ha)*sin(dec))
; lat_r = lat*!pi/180.0
; phi_p = atan(-sin(lat_r)*sin(ha),(cos(ha)*sin(dec)*sin(lat_r)+cos(lat_r)*cos(dec)))
; 


; test how small we can make a convolution kernel for a primary beam before
; the difference between the true beam and one represented by a kernel goes
; above some limit, epsilon.
pro test_pb_kernel_size,freq
    bigsize = 99

    make_azza_arrays_fov,bigsize,az,za,mask,180.0
    tiledata = setup_dipoles(0.0,1,0.0) ; ideal dipoles.
    setup_tile_pointing_direction,0.0*!pi/180,0.0*!pi/180,tiledata
    CALC_DIRECTIONAL_RESPONSE,za,az,freq,tiledata,fieldx,fieldy
    
    ; make shifted FFT of field
    sf = shift(field,-bigsize/2,-bigsize/2) ; shift so that center is at 0,0
    ff = fft(sf)
    sff = shift(ff,bigsize/2,bigsize/2)   ; shift back so peak of fft is in centre of image
    
    bigind = (findgen(bigsize)-bigsize/2)/(bigsize-1.0)*2
    plot,bigind,field[*,bigsize/2]/max(field),xtitle = 'l (direction cosine)',ytitle='voltage beam response OR diff X50',charsize=1.5
    ;tvscl,field
    legend,['Volage beam','True minus 9x9 model X50','True minus 15x15 model X50','True minus 21x21 model X50'],linestyle=findgen(4),/right,charsize=1.0
    
    ;loadct,12
    
    symindex=1
    for siz=9,21,6 do begin
        xind  = findgen(siz)/(siz-1.0)
        mask = fltarr(bigsize,bigsize)
        mask[bigsize/2-siz/2:bigsize/2+siz/2,bigsize/2-siz/2:bigsize/2+siz/2] = 1.0
        ftest = sff*mask
        ;tvscl,ftest
        
        ; make a PB from the truncated kernel
        sftest = shift(ftest,-bigsize/2,-bigsize/2)
        btest = fft(sftest)
        test = shift(btest,bigsize/2,bigsize/2)
        print,'plotting kernel size: ',siz
        oplot,bigind,50*(test[*,bigsize/2]/max(test)-field[*,bigsize/2]/max(field)),linestyle=symindex
        ; stop
        symindex+=1
    endfor
    
end

pro make_pb_kernel_alias_plot,freq
    bigsize = 99
    lilsize = 15

    make_azza_arrays_fov,bigsize,az,za,mask,180.0
    tiledata = setup_dipoles(0.0,1,0.0) ; ideal dipoles.
    setup_tile_pointing_direction,0.0*!pi/180,0.0*!pi/180,tiledata
    CALC_DIRECTIONAL_RESPONSE,za,az,freq,tiledata,fieldx,fieldy
    
    bigind = (findgen(bigsize)-bigsize/2)/(bigsize-1.0)*2
    plot,bigind,field[*,bigsize/2]/max(field),xtitle = 'l (direction cosine)',ytitle='response (Volts)',charsize=1.5,linestyle=1,title=string(format='(%"MWA primary beam at %5.1f MHz.")',freq)
    sm_ind = (findgen(2*(lilsize/2))-lilsize/2)/(lilsize/2)
    sfield = float(field[indgen(2*(lilsize/2))*7,bigsize/2])
    oplot,sm_ind,sfield,linestyle=0

    fsm = fft(shift(sfield,-lilsize/2),/double)
    sfsm= shift(fsm,lilsize/2)
    fb = fltarr(bigsize)
    fb[bigsize/2-lilsize/2:bigsize/2+lilsize/2-1] = sfsm
    sfb = shift(fb,-bigsize/2)
    newb = shift(fft(sfb,/inverse,/double),bigsize/2)
    oplot,bigind,newb,linestyle=2
    oplot,bigind,10.0*(field[*,bigsize/2]-newb),linestyle=3
    legend,['Sampled beam','True beam','Beam after sample aliasing','Error x10'],linestyle=findgen(4),/right,charsize=1.25

end

pro make_pb_spectrum_alias_plot,freq
    bigsize = 99
    lilsize = 15

    make_azza_arrays_fov,bigsize,az,za,mask,180.0
    tiledata = setup_dipoles(0.0,1,0.0) ; ideal dipoles.
    setup_tile_pointing_direction,0.0*!pi/180,0.0*!pi/180,tiledata
    CALC_DIRECTIONAL_RESPONSE,za,az,freq,tiledata,fieldx,fieldy

    vb = field[*,bigsize/2]
    fvb = fft(vb,/double)
    sfield = float(field[indgen(2*(lilsize/2))*7,bigsize/2])
    
    fsm = fft(shift(sfield,-lilsize/2),/double)
    plot,10.*alog10(abs(fvb[0:20])),yrange=[-40,0],charsize=1.5,xtitle="Wavenumber",ytitle="Magnitude of FT of voltage (dB)",title=string(format='(%"MWA primary beam at %5.1f MHz.")',freq)
    oplot,10.*alog10(abs(fsm[0:7])),linestyle=1
    oplot,10.*alog10(abs(abs(fsm[0:7])-abs(fvb[0:7]))),linestyle=2
    legend,['True beam','Sampled beam','Difference'],linestyle=findgen(3),/right,charsize=1.5
end


; stolen from notes:  ~/datadisk/LFD/MWA_RTS/test/TEST_PB_GCF/notes_test.txt
pro make_dipole_volt_plot,freq
    bigsize=100
    l=300./freq
    m = findgen(bigsize)/bigsize*2.0-1.0
    z = findgen(bigsize)/bigsize*!pi-!pi/2
    za=asin(m)
    
    ; response- equally spaced in zenith angle
    e2 = end_fire_2element(2*0.3/l,z)
    e1 = cos(z)
    e3 = e1*e2
    
    ; response- equally spaced in l,m
    el2 = end_fire_2element(2*0.3/l,za)
    el1 = cos(za)
    el3 = el1*el2
    pl1 = el1^2
    pl2 = el2^2
    pl3 = el3^2
    plot,m,el1/max(el1),linestyle=2,xrange=[-1.0,1.0],xtitle='l (direction cosine w.r.t zenith)',ytitle='Normalised voltage response',charsize=2.0,title=string(format='(%"dipole over groundplane voltage response at %5.1f MHz.")',freq)
    oplot,m,el2/max(el2),linestyle=1
    oplot,m,el3/max(el3),linestyle=0
    ;oplot,m,e3/max(e3),linestyle=3
    ;legend,['Total','Groundplane','Dipole','Total (wrong coords)'],linestyle=findgen(4),/right,charsize=1.0
    legend,['Total','Groundplane','Dipole'],linestyle=findgen(3),/center,/bottom,charsize=2.0
end
; postscript,'volt_vs_za.ps'
; make_dipole_volt_plot,140.0
; postscript,/close


function jinc1,x,fwhm_radian
    result=x*0.0
    p = where(abs(x) lt 1e-6)
    if p[0] ne -1 then result[p] = 1.0
    p = where(abs(x) ge 1e-6)
    if p[0] ne -1 then result[p] = 2.0*beselj(!pi/(fwhm_radian)*x[p],1)/x[p]/(!pi/(fwhm_radian))
    return, result
end

pro make_pb_kernel_alias_plot_circularaperture,fov,fwhm
    bigsize=99
    lilsize=16

    fwhm_l = asin(fwhm*!pi/180.0)
    max_l = sin(fov/2.0*!pi/180)
    l = (findgen(2*(bigsize/2))-bigsize/2)/(bigsize/2)*max_l
    ;field = jinc(l,fwhm_l)
    sq_field = fltarr(2*(bigsize/2),2*(bigsize/2))
    xs = fltarr(2*(bigsize/2),2*(bigsize/2))
    ys = fltarr(2*(bigsize/2),2*(bigsize/2))
    for i=0,2*(bigsize/2)-1 do xs[*,i]=l
    for i=0,2*(bigsize/2)-1 do ys[i,*]=l
    sq_field = jinc1(sqrt(xs^2+ys^2),fwhm_l)
    field = sq_field[*,bigsize/2]
    
    ; make shifted FFT of field
    sf = shift(sq_field,-bigsize/2,-bigsize/2) ; shift so that center is at 0,0
    ff = fft(sf,/double)
    sff = shift(ff,bigsize/2,bigsize/2)   ; shift back so peak of fft is in centre of image
    sm_ind = (findgen(2*(lilsize/2))-lilsize/2)/(lilsize/2)*max_l
    sq_sfield = fltarr(2*(lilsize/2),2*(lilsize/2))
    for i=0, 2*(lilsize/2)-1 do sq_sfield[*,i] = sq_field[indgen(2*(lilsize/2))*7,i*7]
    sfield = sq_sfield[*,lilsize/2]
    ls=0

    plot,l,field,yrange=[-0.1,1.0],xtitle = 'l (direction cosine)',ytitle='Normalised response (Volts)',thick=1.0,charsize=2.0,linestyle=ls++,title=string(format='(%"Primary beam FOV: %0.1f degs. FWHM: %0.1f degs.")',fov,fwhm)
    ;plot,sm_ind,sfield/max(field),yrange=[-0.1,1.0],xtitle = 'l (direction cosine)',ytitle='response (Volts)',charsize=1.5,linestyle=ls++,title=string(format='(%"Primary beam at %0.1f GHz. FOV: %5.1f degs. FWHM: %0.2f degs.")',freq/1000.,fov,fwhm)
    oplot,sm_ind,sfield/max(field),psym=1

    fsm = fft(shift(sq_sfield,-lilsize/2,-lilsize/2),/double)
    sfsm= shift(fsm,lilsize/2,lilsize/2)
    fb = fltarr(2*(bigsize/2),2*(bigsize/2))
    fb[bigsize/2-lilsize/2:bigsize/2+lilsize/2-1,bigsize/2-lilsize/2:bigsize/2+lilsize/2-1] = sfsm
    sfb = shift(fb,-bigsize/2,-bigsize/2)
    sq_newb = shift(fft(sfb,/inverse,/double),bigsize/2,bigsize/2)
    newb = sq_newb[*,bigsize/2]
    oplot,l,newb,linestyle=ls++
    oplot,l,10.0*(field-newb),linestyle=ls++
    ;legend,['Sampled beam','True beam','Beam after sample aliasing','Error x10'],linestyle=findgen(4),/right,charsize=2
    ;legend,['True beam','Beam after sample aliasing','Error x10'],linestyle=findgen(ls),/right,/center,charsize=2.0

end
; postscript,'circ_pb_fov1.0_fwhm_1.0.ps'
; make_pb_kernel_alias_plot_circularaperture,1.0,1.0
; postscript,/close
; postscript,'circ_pb_fov2.5_fwhm_1.0.ps'
; make_pb_kernel_alias_plot_circularaperture,2.5,1.0

; postscript,'circ_pb_fov3.25_fwhm_1.0.ps'
; make_pb_kernel_alias_plot_circularaperture,3.25,1.0

; postscript,'circ_pb_fov4.5_fwhm_1.0.ps'
; make_pb_kernel_alias_plot_circularaperture,4.5,1.0


pro make_pb_kernel_alias_plot_2D
    bigsize=256
    lilsize=16
    lat=-26.0
    tile=0

    MAKE_AZZA_ARRAYS_fov,bigsize,az,za,mask,180.0
    MAKE_AZZA_ARRAYS_fov,lilsize,l_az,l_za,l_mask,180.0
    MAKE_JONES_MATRIX_PLOTS,lat,az,za,ra_e,ra_n,dec_e,dec_n
    MAKE_JONES_MATRIX_PLOTS,lat,l_az,l_za,l_ra_e,l_ra_n,l_dec_e,l_dec_n
    gp_effect = end_fire_2element(2*0.3/(300./140.),za)
    gp_effect = gp_effect/max(gp_effect)
    l_gp_effect = end_fire_2element(2*0.3/(300./140.),l_za)
    l_gp_effect = l_gp_effect/max(l_gp_effect)

    tiledata = setup_dipoles(0.0,1,0.0) ; ideal dipoles.
    setup_tile_pointing_direction,0.0*!pi/180,0.0*!pi/180,tiledata
    CALC_DIRECTIONAL_RESPONSE,za,  az,  180.,tiledata,tile_fieldx,fieldy
    CALC_DIRECTIONAL_RESPONSE,l_za,l_az,180.,tiledata,l_tile_fieldx,l_tile_fieldy


    for i=0,3 do begin
        if i eq 0 then begin
            sq_field = ra_e*gp_effect
            l_sq_field = l_ra_e*l_gp_effect
            if tile eq 0 then begin
            fname_true='dipole_ra_e.ps' &  filename = 'dipole_ra_e_pix.ps' & fname_diff='dipole_ra_e_diff.ps'
            endif else begin
            fname_true='tile_ra_e.ps' &  filename = 'tile_ra_e_pix.ps' & fname_diff='tile_ra_e_diff.ps'
            endelse
        endif
        if i eq 1 then begin
            sq_field = ra_n*gp_effect
            l_sq_field = l_ra_n*l_gp_effect
            if tile eq 0 then begin
            fname_true='dipole_ra_n.ps' &  filename = 'dipole_ra_n_pix.ps' & fname_diff='dipole_ra_n_diff.ps'
            endif else begin
            fname_true='tile_ra_n.ps' &  filename = 'tile_ra_n_pix.ps' & fname_diff='tile_ra_n_diff.ps'
            endelse
        endif
        if i eq 2 then begin
            sq_field = dec_e*gp_effect
            l_sq_field = l_dec_e*l_gp_effect
            if tile eq 0 then begin
            fname_true='dipole_dec_e.ps' &  filename = 'dipole_dec_e_pix.ps' & fname_diff='dipole_dec_e_diff.ps'
            endif else begin
            fname_true='tile_dec_e.ps' &  filename = 'tile_dec_e_pix.ps' & fname_diff='tile_dec_e_diff.ps'
            endelse
        endif
        if i eq 3 then begin
            sq_field = dec_n*gp_effect
            l_sq_field = l_dec_n*l_gp_effect
            if tile eq 0 then begin
            fname_true='dipole_dec_n.ps' &  filename = 'dipole_dec_n_pix.ps' & fname_diff='dipole_dec_n_diff.ps'
            endif else begin
            fname_true='tile_dec_n.ps' &  filename = 'tile_dec_n_pix.ps' & fname_diff='tile_dec_n_diff.ps'
            endelse
        endif

        if tile ne 0 then begin
            sq_field *= tile_field
            l_sq_field*=l_tile_field
        endif

        ; make shifted FFT of field
        sf  = shift(l_sq_field,-lilsize/2,-lilsize/2) ; shift so that center is at 0,0
        ff  = fft(sf,/double)
        sff = shift(ff,lilsize/2,lilsize/2)         ; shift back so peak of fft is in centre of image

        big = complexarr(bigsize,bigsize)
        big[(bigsize/2-lilsize/2):(bigsize/2+lilsize/2)-1,(bigsize/2-lilsize/2):(bigsize/2+lilsize/2)-1]=sff
        sbig = shift(big,-bigsize/2,-bigsize/2)
        fbig = fft(sbig,/inverse)
        sfbig = shift(fbig,bigsize/2,bigsize/2)
    
        window,xsize=bigsize,ysize=bigsize+75
        tvscl,sfbig
        colorbar,min=min(sfbig),max=max(sfbig),charsize=1.5,format='(F6.1)',divisions=1
        postscript,filename
        tvscl,sfbig
        colorbar,min=min(sfbig),max=max(sfbig),charsize=1.5,format='(F6.1)',divisions=1
        postscript,/close

        window,xsize=bigsize,ysize=bigsize+75
        ;contour,sq_field,xmargin=[0,0],ymargin=[0,4],xrange=[0,256],xstyle=4+1,yrange=[0,256],ystyle=4+1,/nodata
        tvscl,sq_field
        ;contour,sq_field+(1.0-mask)*10.,xmargin=[0,0],ymargin=[0,4],xrange=[0,256],xstyle=4+1,yrange=[0,256],ystyle=4+1,levels=[0.0],/overplot
        colorbar,min=min(sq_field),max=max(sq_field),charsize=1.5,format='(F6.1)',divisions=1
        postscript,fname_true
        tvscl,sq_field
        colorbar,min=min(sq_field),max=max(sq_field),charsize=1.5,format='(F6.1)',divisions=1
        postscript,/close
        
        window,xsize=bigsize,ysize=bigsize+75
        tvscl,sfbig-sq_field
        colorbar,min=min(sfbig-sq_field),max=max(sfbig-sq_field),charsize=1.5,format='(F6.2)',divisions=1
        postscript,fname_diff
        tvscl,sfbig-sq_field
        colorbar,min=min(sfbig-sq_field),max=max(sfbig-sq_field),charsize=1.5,format='(F6.2)',divisions=1
        postscript,/close
        
    endfor
end


; make plots of the direction dependent Jones matrix elements
; example:
; MAKE_AZZA_ARRAYS_fov,200,az,za,mask,180.0
; MAKE_JONES_MATRIX_PLOTS,-26.7,az,za,ra_e,ra_n,dec_e,dec_n
; freq = 140.0
; calc_directional_response,za,az,freq,tiledata,fieldx,fieldy
; gp_effect = end_fire_2element(2*0.3/(300./freq),za)
; gp_effect = gp_effect/max(gp_effect)
; tvscl,ra_e*gp_effect*mask
; t=ra_e*gp_effect*mask
; colorbar,
pro make_jones_matrix_plots,lat,az,za,ra_e,ra_n,dec_e,dec_n
    s = size(az)
    ha_img  = fltarr(s[1],s[2])
    dec_img = ha_img
    
    for j=0,s[2]-1 do begin
        for i=0, s[1]-1 do begin
            ALTAZ2HADEC,90.0-za[i,j]*180/!pi,az[i,j]*180/!pi,lat,ha,dec
            ha_img[i,j] = ha
            dec_img[i,j]=dec
        endfor
    endfor

    ra_e = cos(ha_img*!pi/180.)
    ra_n = -sin(lat*!pi/180)*sin(ha_img*!pi/180.)
    dec_e= sin(dec_img*!pi/180.0)*sin(ha_img*!pi/180.)
    dec_n= cos(lat*!pi/180)*cos(dec_img*!pi/180.0) + sin(lat*!pi/180)*sin(dec_img*!pi/180.0)*cos(ha_img*!pi/180.)
end


; load impedance simulation data. The dipoles are in the following order:
; 1-16 Y (north-south) dipoles, 17-32 X (east-west) dipoles.
; looking down on tile, it runs from top left across to top right, then down through
; the rows left to right.
pro load_MWA_z,z
  freq='216'
  filename='/tmp/TotalZ_real_'+freq+'.dat'
  r = read_ascii(filename)
  filename='/tmp/TotalZ_imag_'+freq+'.dat'
  i = read_ascii(filename)
  z = complex(r.field1,i.field1)
  ; invert the ordering on the Y axis to match labelling convention in this code ?
end

; calculate zenith-pointed directivity
function calc_directivity,za,az,field
  p = where(za lt !pi/2)
  sky_area_factor = za*0.0
  sky_area_factor[p] = 1./cos(za[p])
  pow = float(field*conj(field))
  m = max(pow)
  av = total(pow*sky_area_factor)/(2.0*total(sky_area_factor))
  return,m/av
end

; calculate zenith-pointed effective area
function calc_Aeff,za,az,field,freq
  lambda=300./freq
  d = calc_directivity(za,az,field)
  return,lambda*lambda*d/(4.0*!pi)
end

