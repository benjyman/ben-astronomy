
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

;extra idl stuff
pro mkhdr, header, im, naxisx, IMAGE = image, EXTEND = extend
;+
; NAME:
;       MKHDR
; PURPOSE:
;       Make a minimal primary (or IMAGE extension) FITS header
; EXPLANATION:
;       If an array is supplied,  then the created FITS header will be 
;       appropriate to the supplied array.  Otherwise, the user can specify 
;       the dimensions and datatype.
;
;       To update an *existing* FITS header with a new image array, instead 
;       use check_FITS, /Update 
;
; CALLING SEQUENCE:
;       MKHDR, header                   ;Prompt for image size and type
;               or
;       MKHDR, header, im, [ /IMAGE, /EXTEND ]
;               or
;       MKHDR, header, type, naxisx, [/IMAGE, /EXTEND ]         
;
; OPTIONAL INPUTS:
;       IM - If IM is a vector or array then the header will be made
;               appropriate to the size and type of IM.  IM does not have
;               to be the actual data; it can be a dummy array of the same
;               type and size as the data.    Set IM = '' to create a dummy
;               header with NAXIS = 0. 
;       TYPE - If 2 parameters are supplied, then the second parameter
;               is interpreted as an integer giving the IDL datatype e.g. 
;               1 - Byte, 2 - 16 bit integer, 4 - float, 3 - Long
;       NAXISX - Vector giving the size of each dimension (NAXIS1, NAXIS2, 
;               etc.).  
;
; OUTPUT:
;       HEADER - image header, (string array) with required keywords
;               BITPIX, NAXIS, NAXIS1, ... Further keywords can be added
;               to the header with SXADDPAR. 
;
; OPTIONAL INPUT KEYWORDS:
;       /IMAGE   = If set, then a minimal header for a FITS IMAGE extension
;               is created.    An IMAGE extension header is identical to
;               a primary FITS header except the first keyword is 
;               'XTENSION' = 'IMAGE' instead of 'SIMPLE  ' = 'T'
;       /EXTEND  = If set, then the keyword EXTEND is inserted into the file,
;               with the value of "T" (true).    The EXTEND keyword can 
;               optionally be included in a primary header, if the FITS file 
;               contains extensions.
;
; RESTRICTIONS:
;       (1)  MKHDR should not be used to make an STSDAS header or a FITS
;               ASCII or Binary Table extension header.   Instead use
;
;               SXHMAKE - to create a minimal STSDAS header
;               FXBHMAKE - to create a minimal FITS binary table header
;               FTCREATE - to create a minimal FITS ASCII table header
;
;       (2)  Any data already in the header before calling MKHDR
;               will be destroyed.
; EXAMPLE:
;       Create a minimal FITS header, Hdr, for a 30 x 40 x 50 INTEGER*2 array
;
;             IDL> mkhdr, Hdr, 2, [30,40,50]
;
;       Alternatively, if the array already exists as an IDL variable, Array,
;
;              IDL> mkhdr, Hdr, Array
;
; PROCEDURES CALLED:
;       SXADDPAR, GET_DATE
;
; REVISION HISTORY:
;       Written November, 1988               W. Landsman
;       May, 1990, Adapted for IDL Version 2.0, J. Isensee
;       Aug, 1997, Use SYSTIME(), new DATE format  W. Landsman
;       Allow unsigned data types    W. Landsman   December 1999
;       Set BZERO = 0 for unsigned integer data  W. Landsman January 2000
;       EXTEND keyword must immediately follow last NAXISi W. Landsman Sep 2000
;       Add FITS definition COMMENT to primary headers W. Landsman Oct. 2001
;       Allow (nonstandard) 64 bit integers   W. Landsman  Feb. 2003
;       Add V6.0 notation W. Landsman July 2012
;       Support unsigned 64 bit integers W. Landsman January 2018 
;-                          
 compile_opt idl2

 npar = N_params()
 if npar LT 1 then begin
   print,'Syntax:  MKHDR, header, [ im, /IMAGE, /EXTEND ]'
   print,'    or   MKHDR, header, [ type, naxisx, /IMAGE, /EXTEND ]'
   print,'   header - output FITS header to be created'
   return
 endif
 
 Catch, theError
IF theError NE 0 then begin
	Catch,/Cancel
	void = cgErrorMsg(/quiet)
	RETURN
ENDIF

 if (npar eq 1) then begin               ;Prompt for keyword values
    read,'Enter number of dimensions (NAXIS): ',naxis
    s = lonarr(naxis+2)
    s[0] = naxis
    if ( naxis GT 0 ) then begin       ;Make sure not a dummy header
    for i = 1,naxis do begin       ;Get dimension of each axis
      keyword = 'NAXIS' + strtrim(i,2)
      read,'Enter size of dimension '+ strtrim(i,2) + ' ('+keyword+'): ',nx
      s[i] = nx                            
    endfor
  endif

  print,'Allowed datatypes are (1) Byte, (2) 16 bit integer, (3) 32 bit integer'
  print,'                      (4) 32bit floating, (5) 64 bit double precision' 
  print,'                  or (14) 64bit integer'
  read,'Enter datatype: ',stype
  s[s[0] + 1] = stype

 endif else $
     if ( npar EQ 2 ) then s = size(im) $  ;Image array supplied
          else  s = [ N_elements(naxisx),naxisx, im ] ;Keyword values supplied

 stype = s[s[0]+1]              ;Type of data    
        case stype of
        0:      message,'ERROR: Input data array is undefined'
        1:      bitpix = 8  
        2:      bitpix = 16  
        3:      bitpix = 32  
        4:      bitpix = -32 
        5:      bitpix = -64 
        6:      message,'Complex types not allowed as FITS primary arrays'  
        7:      bitpix = 8
       12:      bitpix = 16
       13:      bitpix = 32
       14:      bitpix = 64
       15:      bitpix = 64
        else:   message,'ERROR: Illegal Image Datatype'
        endcase

 header = strarr(s[0] + 7) + string(' ',format='(a80)')      ;Create empty array
 header[0] = 'END' + string(replicate(32b,77))

 if keyword_set( IMAGE) then $
    sxaddpar, header, 'XTENSION', 'IMAGE   ',' IMAGE extension' $
 else $
    sxaddpar, header, 'SIMPLE', 'T',' Written by IDL:  '+ systime()

 sxaddpar, header, 'BITPIX', bitpix, ' Number of bits per data pixel'
 sxaddpar, header, 'NAXIS', S[0],' Number of data axes'       ;# of dimensions

 if ( s[0] GT 0 ) then begin
   for i = 1, s[0] do sxaddpar,header,'NAXIS' + strtrim(i,2),s[i]
 endif

 if keyword_set( IMAGE) then begin
     sxaddpar, header, 'PCOUNT', 0, ' No Group Parameters'
     sxaddpar, header, 'GCOUNT', 1, ' One Data Group'
 endif else begin
     if keyword_set( EXTEND) or (s[0] EQ 0) then $
          sxaddpar, header, 'EXTEND', 'T', ' FITS data may contain extensions'
     Get_date, dte                       ;Get current date as CCYY-MM-DD
     sxaddpar, header, 'DATE', dte, $
           ' Creation UTC (CCCC-MM-DD) date of FITS header'
  endelse

 if stype EQ 12 then sxaddpar, header,'O_BZERO',32768, $
            ' Original Data is Unsigned Integer'
 if stype EQ 13 then sxaddpar, header,'O_BZERO',2147483648, $
            ' Original Data is Unsigned Long'
 if stype EQ 15 then sxaddpar, header,'O_BZERO',ulong64(2)^63, $
            ' Original Data is Unsigned 64 bit Long'           
 header = header[0:s[0]+7]

 if ~keyword_set(IMAGE) then begin   ;Add FITS definition for primary header
     sxaddpar,header,'COMMENT ', $
      "FITS (Flexible Image Transport System) format is defined in 'Astronomy"
     sxaddpar,header,'COMMENT ', $
      "and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H"
 endif 
 end
 
 
 ; docformat = 'rst'
;
; NAME:
;   cgErrorMsg
;
; PURPOSE:
;   The purpose of this function is to have a device-independent error messaging function. 
;   The error message is reported to the user by using DIALOG_MESSAGE if widgets are
;    supported. Otherwise, it is just printed to standard out.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2013, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; The purpose of this function is to have a device-independent error messaging function.  
; The error message is reported to the user by using DIALOG_MESSAGE if widgets are
; supported. Otherwise, the error message is just printed to standard out.
;
; :Categories:
;    Utilities
;
; :Params:
;    themessage: in, optional, type=string
;       This is a string argument containing the error message you want reported. If undefined, 
;       this variable is set to the string in the !Error_State.Msg system variable.
;       
; :Keywords:
;    error: in, optional, type=boolean, default=1
;       Set this keyword to cause Dialog_Message to use the ERROR reporting dialog. 
;       Note that a longstanding bug in IDL causes the ERROR dialog to be used whether 
;       this keyword is set to 0 or 1!
;    informational: in, optional, type=boolean, default=0
;       Set this keyword to cause Dialog_Message to use the INFORMATION dialog instead of the 
;       WARNING dialog. Note that a bug in IDL causes the ERROR dialog to be used if this keyword 
;       is set to 0!
;    noname: in, optional, type=boolean, default=0
;       Normally, the name of the routine in which the error occurs is prepended to the error
;       message. Setting this keyword will suppress this behavior and no routine name will be used.
;    quiet: in, optional, type=boolean, default=0
;       Set this keyword to suppress the DIALOG_MESSAGE pop-up dialog.
;    title: in, optional, type=string
;       Set this keyword to the title of the DIALOG_MESSAGE window. By default the keyword is set to 
;       'System Error' unless !ERROR_STATE.NAME equals "IDL_M_USER_ERR", in which case it is set to 
;       "Trapped Error'.
;     traceback: in, optional, type=boolean, default=1
;        Setting this keyword results in an error traceback being printed to standard output with the 
;        PRINT command. Use TRACEBACK=0 to turn this functionality off.
;        
; :Examples:
;   In general, the cgErrorMsg function is not called directly. Rather, it is used in a 
;   CATCH error handler. Errors are thrown to cgErrorMsg with the MESSAGE command. A typical 
;   CATCH error handler is shown below::
;
;       Catch, theError
;       IF theError NE 0 THEN BEGIN
;          Catch, /Cancel
;          void = cgErrorMsg()
;          RETURN
;       ENDIF
;
;   Error messages will get into the cgErrorMsg function by throwing an error with the 
;   MESSAGE command, like this::
;
;       IF test NE 1 THEN Message, 'The test failed.'
;
; :Author:
;    FANNING SOFTWARE CONSULTING::
;       David W. Fanning
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: david@idlcoyote.com
;       Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;    Change History::
;      Written by: David W. Fanning, 27 April 1999.
;      Added the calling routine's name in the message and NoName keyword. 31 Jan 2000. DWF.
;      Added _Extra keyword. 10 February 2000. DWF.
;      Forgot to add _Extra everywhere. Fixed for MAIN errors. 8 AUG 2000. DWF.
;      Adding call routine's name to Traceback Report. 8 AUG 2000. DWF.
;      Added ERROR, INFORMATIONAL, and TITLE keywords. 19 SEP 2002. DWF.
;      Removed the requirement that you use the NONAME keyword with the MESSAGE
;        command when generating user-trapped errors. 19 SEP 2002. DWF.
;      Added distinctions between trapped errors (errors generated with the
;        MESSAGE command) and IDL system errors. Note that if you call ERROR_MESSAGE
;        directly, then the state of the !ERROR_STATE.NAME variable is set
;        to the *last* error generated. It is better to access ERROR_MESSAGE
;        indirectly in a Catch error handler from the MESSAGE command. 19 SEP 2002. DWF.
;      Change on 19 SEP 2002 to eliminate NONAME requirement did not apply to object methods.
;        Fixed program to also handle messages from object methods. 30 JULY 2003. DWF.
;      Removed obsolete STR_SEP and replaced with STRSPLIT. 27 Oct 2004. DWF.
;      Made a traceback the default case without setting TRACEBACK keyword. 19 Nov 2004. DWF.
;      Added check for window connection specifically for CRON jobs. 6 May 2008. DWF.
;      Added QUIET keyword. 18 October 2008. DWF.
;      The traceback information was bypassed when in the PostScript device. Not what I
;        had in mind. Fixed. 6 July 2009. DWF.
;      The QUIET keyword was clearing traceback information. Fixed with help from Phillip Bitzer. 2 Oct 2012. DWF.
;      Name changed to cgErrorMsg from Error_Message, 4 November 2013, by David W. Fanning.
;
; :Copyright:
;     Copyright (c) 2013, Fanning Software Consulting, Inc.
;-
FUNCTION cgErrorMsg, theMessage, Error=error, Informational=information, $
   Traceback=traceback, NoName=noname, Title=title, Quiet=quiet, _Extra=extra

   On_Error, 2
   
   ; Check for presence and type of message.
   
   IF N_Elements(theMessage) EQ 0 THEN theMessage = !Error_State.Msg
   s = Size(theMessage)
   messageType = s[s[0]+1]
   IF messageType NE 7 THEN BEGIN
      Message, "The message parameter must be a string.", _Extra=extra
   ENDIF
   IF N_Elements(traceback) EQ 0 THEN traceback = 1
   
   ; Get the call stack and the calling routine's name.
   Help, Calls=callStack
   callingRoutine = (StrSplit(StrCompress(callStack[1])," ", /Extract))[0]
   
   ; Are widgets supported?
   IF !D.Name EQ 'PS' OR !D.Name EQ 'Z' THEN BEGIN
      widgetsSupported = 1
   ENDIF ELSE BEGIN
      widgetsSupported = ((!D.Flags AND 65536L) NE 0)
   ENDELSE

   ; Is the QUIET keyword set? Then no dialogs.
   IF Keyword_Set(quiet) THEN widgetsSupported = 0
   
   ; It is not enough to know if widgets are supported. In CRON jobs, widgets are
   ; supported, but there is no X connection and pop-up dialogs are not allowed.
   ; Here is a quick test to see if we can connect to a windowing system. If not,
   ; then we are going to assume widgets are not supported.
   Catch, theError
   IF theError NE 0 THEN BEGIN
      Catch, /CANCEL
      widgetsSupported = 0
      GOTO, testWidgetSupport
   ENDIF
   theWindow = !D.Window
   IF (!D.Flags AND 256) NE 0 THEN Window, /FREE, XSIZE=5, YSIZE=5, /PIXMAP
   Catch, /CANCEL
   
   testWidgetSupport: ; Come here if you choke on creating a window.
   IF !D.Window NE theWindow THEN BEGIN
      WDelete, !D.Window
      IF theWindow GE 0 THEN WSet, theWindow
   ENDIF
   
   IF widgetsSupported THEN BEGIN
   
      ; If this is an error produced with the MESSAGE command, it is a trapped
      ; error and will have the name "IDL_M_USER_ERR".
      IF !ERROR_STATE.NAME EQ "IDL_M_USER_ERR" THEN BEGIN
   
         IF N_Elements(title) EQ 0 THEN title = 'Trapped Error'
   
         ; If the message has the name of the calling routine in it,
         ; it should be stripped out. Can you find a colon in the string?
   
         ; Is the calling routine an object method? If so, special processing
         ; is required. Object methods will have two colons together.
         doublecolon = StrPos(theMessage, "::")
         IF doublecolon NE -1 THEN BEGIN
   
            prefix = StrMid(theMessage, 0, doublecolon+2)
            submessage = StrMid(theMessage, doublecolon+2)
            colon = StrPos(submessage, ":")
            IF colon NE -1 THEN BEGIN
   
               ; Extract the text up to the colon. Is this the same as
               ; the callingRoutine? If so, strip it.
               IF StrMid(theMessage, 0, colon+StrLen(prefix)) EQ callingRoutine THEN $
                  theMessage = StrMid(theMessage, colon+1+StrLen(prefix))
            ENDIF
         ENDIF ELSE BEGIN
   
            colon = StrPos(theMessage, ":")
            IF colon NE -1 THEN BEGIN
   
               ; Extract the text up to the colon. Is this the same as
               ; the callingRoutine? If so, strip it.
               IF StrMid(theMessage, 0, colon) EQ callingRoutine THEN $
                  theMessage = StrMid(theMessage, colon+1)
            ENDIF
   
         ENDELSE
   
   
         ; Add the calling routine's name, unless NONAME is set.
         IF Keyword_Set(noname) THEN BEGIN
            answer = Dialog_Message(theMessage, Title=title, _Extra=extra, $
               Error=error, Information=information)
         ENDIF ELSE BEGIN
            answer = Dialog_Message(StrUpCase(callingRoutine) + ": " + $
               theMessage, Title=title, _Extra=extra, $
               Error=error, Information=information)
         ENDELSE
   
      ENDIF ELSE BEGIN
   
         ; Otherwise, this is an IDL system error.
         IF N_Elements(title) EQ 0 THEN title = 'System Error'
   
         IF StrUpCase(callingRoutine) EQ "$MAIN$" THEN $
            answer = Dialog_Message(theMessage, _Extra=extra, Title=title, $
               Error=error, Information=information) ELSE $
         IF Keyword_Set(noname) THEN BEGIN
            answer = Dialog_Message(theMessage, _Extra=extra, Title=title, $
               Error=error, Information=information)
         ENDIF ELSE BEGIN
            answer = Dialog_Message(StrUpCase(callingRoutine) + "--> " + $
               theMessage, _Extra=extra, Title=title, $
               Error=error, Information=information)
         ENDELSE
      ENDELSE
   ENDIF ELSE BEGIN
         Help, /Last_Message, Output=traceback_msg ; Required because following MESSAGE call clears traceback info.
         Message, theMessage, /Continue, /NoPrint, /NoName, /NoPrefix, _Extra=extra
         IF Keyword_Set(noname) THEN $
            Print, theMessage ELSE $
            Print, '%' + callingRoutine + ': ' + theMessage 
         answer = 'OK'
   ENDELSE
   
   ; Provide traceback information if requested and this is NOT an informational message.
   IF Keyword_Set(traceback) AND ~Keyword_Set(informational)THEN BEGIN
      IF N_Elements(traceback_msg) NE 0 THEN traceback = traceback_msg ELSE Help, /Last_Message, Output=traceback
      Print,''
      Print, 'Traceback Report from ' + StrUpCase(callingRoutine) + ':'
      Print, ''
      FOR j=0,N_Elements(traceback)-1 DO Print, "     " + traceback[j]
   ENDIF
   
   RETURN, answer
END ; ----------------------------------------------------------------------------


Pro sxaddpar, Header, Name, Value, Comment, Location, before=before, $
                 savecomment = savecom, after=after , format=format, pdu = pdu, $
                 missing = missing, null = null
;+
; NAME:
;       SXADDPAR
; PURPOSE:
;       Add or modify a parameter in a FITS header array.
;
; CALLING SEQUENCE:
;       SXADDPAR, Header, Name, Value, [ Comment,  Location, /SaveComment, 
;                               BEFORE =, AFTER = , FORMAT= , /PDU
;                               /SAVECOMMENT, Missing=, /Null
; INPUTS:
;       Header = String array containing FITS header.    The
;               length of each element must be 80 characters.    If not 
;               defined, then SXADDPAR will create an empty FITS header array.
;
;       Name = Name of parameter. If Name is already in the header the value 
;               and possibly comment fields are modified.  Otherwise a new 
;               record is added to the header.  If name is equal to 'COMMENT'
;               or 'HISTORY' or a blank string then the value will be added to 
;               the record without replacement.  For these cases, the comment 
;               parameter is ignored.
;
;       Value = Value for parameter.  The value expression must be of the 
;               correct type, e.g. integer, floating or string.  String values
;                of 'T' or 'F' are considered logical values.
;
; OPTIONAL INPUT PARAMETERS:
;       Comment = String field.  The '/' is added by this routine.  Added 
;               starting in position 31.    If not supplied, or set equal to 
;               '', or /SAVECOMMENT is set, then the previous comment field is 
;               retained (when found) 
;
;       Location = Keyword string name.  The parameter will be placed before the
;               location of this keyword.    This parameter is identical to
;               the BEFORE keyword and is kept only for consistency with
;               earlier versions of SXADDPAR.
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       BEFORE  = Keyword string name.  The parameter will be placed before the
;               location of this keyword.  For example, if BEFORE='HISTORY'
;               then the parameter will be placed before the first history
;               location.  This applies only when adding a new keyword;
;               keywords already in the header are kept in the same position.
;
;       AFTER   = Same as BEFORE, but the parameter will be placed after the
;               location of this keyword.  This keyword takes precedence over
;               BEFORE.
;
;       FORMAT  = Specifies FORTRAN-like format for parameter, e.g. "F7.3".  A
;               scalar string should be used.  For complex numbers the format
;               should be defined so that it can be applied separately to the
;               real and imaginary parts.  If not supplied then the default is
;               'G19.12' for double precision, and 'G14.7' for floating point.
;       /NULL   = If set, then keywords with values which are undefined, or
;                 which have non-finite values (such as NaN, Not-a-Number) are
;                 stored in the header without a value, such as
;
;                       MYKEYWD =                      /My comment
;
;       MISSING = A value which signals that data with this value should be
;                 considered missing.  For example, the statement
;
;                       FXADDPAR, HEADER, 'MYKEYWD', -999, MISSING=-999
;
;                 would result in the valueless line described above for the
;                 /NULL keyword.  Setting MISSING to a value implies /NULL.
;                 Cannot be used with string or complex values.
;       /PDU    = specifies keyword is to be added to the primary data unit
;               header. If it already exists, it's current value is updated in
;               the current position and it is not moved.
;       /SAVECOMMENT = if set, then any existing comment is retained, i.e. the
;               COMMENT parameter only has effect if the keyword did not 
;               previously exist in the header.
; OUTPUTS:
;       Header = updated FITS header array.
;
; EXAMPLE:
;       Add a keyword 'TELESCOP' with the value 'KPNO-4m' and comment 'Name
;       of Telescope' to an existing FITS header h.
;
;       IDL> sxaddpar, h, 'TELESCOPE','KPNO-4m','Name of Telescope'
; NOTES:
;       The functions SXADDPAR() and FXADDPAR() are nearly identical, with the
;       major difference being that FXADDPAR forces required FITS keywords
;       BITPIX, NAXISi, EXTEND, PCOUNT, GCOUNT to appear in the required order
;       in the header, and FXADDPAR supports the OGIP LongString convention.   
;       There is no particular reason for having two nearly identical 
;       procedures, but both are too widely used to drop either one.
;
;       All HISTORY records are inserted in order at the end of the header.
;
;       All COMMENT records are also inserted in order at the end of the header
;       header, but before the HISTORY records.  The BEFORE and AFTER keywords
;       can override this.
;
;       All records with no keyword (blank) are inserted in order at the end of
;       the header, but before the COMMENT and HISTORY records.  The BEFORE and
;       AFTER keywords can override this.

; RESTRICTIONS:
;       Warning -- Parameters and names are not checked
;               against valid FITS parameter names, values and types.
;
; MODIFICATION HISTORY:
;       DMS, RSI, July, 1983.
;       D. Lindler Oct. 86  Added longer string value capability
;       Added Format keyword, J. Isensee, July, 1990
;       Added keywords BEFORE and AFTER. K. Venkatakrishna, May '92
;       Pad string values to at least 8 characters   W. Landsman  April 94
;       Aug 95: added /PDU option and changed routine to update last occurrence
;               of an existing keyword (the one SXPAR reads) instead of the
;               first occurrence.
;       Comment for string data can start after column 32 W. Landsman June 97
;       Make sure closing quote supplied with string value  W. Landsman  June 98
;       Increase precision of default formatting of double precision floating
;               point values.   C. Gehman, JPL  September 1998
;       Mar 2000, D. Lindler, Modified to use capital E instead of lower case
;               e for exponential formats.
;       Apr 2000, Make user-supplied format upper-case  W. Landsman 
;       Oct 2001, Treat COMMENT or blank string like HISTORY keyword W. Landsman
;       Jan 2002, Allow BEFORE, AFTER to apply to COMMENT keywords W. Landsman
;       June 2003, Added SAVECOMMENT keyword    W. Landsman
;       Jan 2004, If END is missing, then add it at the end W. Landsman
;       May 2005 Fix SAVECOMMENT error with non-string values W. Landsman
;       Oct 2005 Jan 2004 change made SXADDPAR fail for empty strings W.L.
;       May 2011 Fix problem with slashes in string values W.L. 
;       Aug 2013 Only use keyword_set for binary keywords W. L. 
;       Sep 2015 Added NULL and MISSING keywords W.L.
;       Sep 2016 Allow writing of byte or Boolean variables  W.L.
;       Nov 2016 Allow value to be a 1 element scalar  W.L.
;       Jun 2018 W. Thompson, for backward compatibility, save non-finite
;                values (e.g. NaN) as strings if /NULL not set
;       
;-
 compile_opt idl2
 if N_params() LT 3 then begin             ;Need at least 3 parameters
      print,'Syntax - Sxaddpar, Header, Name,  Value, [Comment, Postion'
      print,'                      BEFORE = ,AFTER = , FORMAT =, /SAVECOMMENT'
      print,'                      MISSING =, /NULL'
      return
 endif

; Define a blank line and the END line

 ENDLINE = 'END' +string(replicate(32b,77))     ;END line
 BLANK = string(replicate(32b,80))             ;BLANK line
;
;  If Location parameter not defined, set it equal to 'END     '
;
 if ( N_params() GT 4 ) then loc = strupcase(location) else $
 if N_elements( BEFORE) GT 0 then loc = strupcase(before) else $
 if N_elements( AFTER) GT 0  then loc = strupcase(after) else $
 if N_elements( PDU) GT 0  then loc = 'BEGIN EX' else $
                             loc = 'END'

 while strlen(loc) lt 8 do loc += ' '

 if N_params() lt 4 then comment = ''      ;Is comment field specified?

 n = N_elements(header)                  ;# of lines in FITS header
 if (n EQ 0) then begin                  ;header defined?
          header=strarr(10)              ;no, make it.
          header[0]=ENDLINE
          n=10
 endif else begin
          s = size(header)               ;check for string type
              if (s[0] ne 1) || (s[2] ne 7) then $
                  message,'FITS Header (first parameter) must be a string array'
 endelse

;  Make sure Name is 8 characters long

        nn = string(replicate(32b,8))   ;8 char name
        strput,nn,strupcase(name) ;insert name
;
;  Check to see if the parameter should be saved as a null value.
;
        stype = size(value,/type)
        save_as_null = 0
        if stype EQ 0 then begin        
            if (n_elements(missing) eq 1) || keyword_set(null) then $
              save_as_null = 1 else $
                message,'Keyword value (third parameter) is not defined'
        endif else if (stype NE 6) && (stype NE 7) && (stype NE 9) then begin
            if N_elements(missing) eq 1 then $
              if value eq missing then save_as_null = 1
              if ~save_as_null then if ~finite(value) then begin
                if ((n_elements(missing) eq 1) || keyword_set(null)) then $
                  save_as_null = 1 else begin
                    message,/continue,'Keyword value (third parameter) ' + $
                            'is not finite, saving as string.'
                    stype = 7
                    save_as_string = 1
                endelse
            endif
        endif
;
;  Extract first 8 characters of each line of header, and locate END line

 keywrd = strmid(header,0,8)                 ;Header keywords
 iend = where(keywrd eq 'END     ',nfound)
;
;  If no END, then add it.  Either put it after the last non-null string, or
;  append it to the end.
;
        if nfound EQ 0 then begin
                ii = where(strtrim(header) ne '',nfound)
                ii = max(ii) + 1
                if ii eq n_elements(header) then begin
                        header = [header,endline]
                        n++ 
                endif else header[ii] = endline
                keywrd = strmid(header,0,8)
                iend = where(keywrd eq 'END     ',nfound)
        endif
;
        iend = iend[0] > 0                      ;make scalar

;  History, comment and "blank" records are treated differently from the
;  others.  They are simply added to the header array whether there are any
;  already there or not.

 if (nn EQ 'HISTORY ') || (nn EQ 'COMMENT ') || $
    (nn EQ '        ')  then begin             ;add history record?
;
;  If the header array needs to grow, then expand it in increments of 5 lines.
;

     if iend GE (n-1) then begin
                 header = [header,replicate(blank,5)] ;yes, add 5.
                 n = N_elements(header)
      endif

; Format the record

      newline = blank
      strput,newline,nn+string(value),0

;
;  If a history record, then append to the record just before the end.
;
      if nn EQ 'HISTORY ' then begin
             header[iend] = newline             ;add history rec.
             header[iend+1] = endline
;
;  The comment record is placed immediately after the last previous comment
;  record, or immediately before the first history record, unless overridden by
;  either the BEFORE or AFTER keywords.
;
      endif else if nn EQ 'COMMENT ' then begin
            if loc EQ 'END     ' then loc = 'COMMENT '
            iloc = where(keywrd EQ loc, nloc)
            if nloc EQ 0 then iloc = where(keywrd EQ 'HISTORY ', nloc)
            if nloc gt 0 then begin
               i = iloc[nloc-1]
               if keyword_set(after) or (loc EQ 'COMMENT ') then i = i+1 < iend 
               if i gt 0 then header=[header[0:i-1],newline,header[i:n-1]] $
                        else header=[newline,header[0:n-1]]
            endif else begin
                header[iend] = newline
                header[iend+1] = endline
            endelse

;
;  The "blank" record is placed immediately after the last previous "blank"
;  record, or immediately before the first comment or history record, unless
;  overridden by either the BEFORE or AFTER keywords.
;
          ENDIF ELSE BEGIN
            if loc EQ 'END     ' then loc = '       '
            iloc = where(keywrd[0:iend] EQ loc, nloc)
            if nloc gt 0 then begin
               i = iloc[0]
               if keyword_set(after) and loc ne 'HISTORY ' then i = i+1 < iend 
               if i gt 0 then header=[header[0:i-1],newline,header[i:n-1]] $
                        else header=[newline,header[0:n-1]]
            endif else begin
                iloc = where(keywrd EQ 'COMMENT ', nloc)
                if nloc Eq 0 then iloc = where(keywrd EQ 'HISTORY ', nloc)
                if nloc GT 0 then begin
                   i = iloc[0]
                   if i gt 0 then header=[header[0:i-1],newline,header[i:n-1]] $
                        else header=[newline,header[0:n-1]]
                endif else begin
                  header[iend] = newline
                  header[iend+1] = endline
            endelse
            endelse
           endelse
            RETURN
 endif

; Find location to insert keyword.   Save the existing comment if user did
; not supply a new one.   Comment starts after column 32 for numeric data,
; after the slash (but at least after final quote) for string data. 

 ncomment = comment
 ipos  = where(keywrd eq nn,nfound)
 if nfound gt 0 then begin
         i = ipos[nfound-1]
         if comment eq '' or keyword_set(savecom) then begin  ;save comment?
         if strmid(header[i],10,1) NE "'" then $
                 ncomment=strmid(header[i],32,48) else begin
		 quote = strpos(header[i],"'",11)
		
                 if quote EQ -1 then slash = -1 else $
		       slash = strpos(header[i],'/',quote)  		
                 if slash NE -1 then $
                        ncomment =  strmid(header[i], slash+1, 80) else $
                        ncomment = string(replicate(32B,80))
                endelse
        endif 
         goto, REPLACE    
 endif

 if loc ne '' then begin
          iloc =  where(keywrd eq loc,nloc)
          if nloc gt 0 then begin
             i = iloc[0]
             if keyword_set(after) && (loc ne 'HISTORY ') then i = i+1 < iend 
             if i gt 0 then header=[header[0:i-1],blank,header[i:n-1]] $
                        else header=[blank,header[0:n-1]]
             goto, REPLACE  
          endif
 endif

; At this point keyword and location parameters were not found, so a new
; line is added at the end of the FITS header

        if iend lt (n-1) then begin     ;Not found, add more?
                header[iend+1] = ENDLINE        ;no, already long enough.
                i = iend                ;position to add.
           endif else begin             ;must lengthen.
                header = [header,replicate(blank,5)] ;add an element on the end
                header[n]=ENDLINE               ;save "END"
                i =n-1                  ;add to end
        end

; Now put value into keyword at line i

REPLACE:    
        h=blank                 ;80 blanks
        strput,h,nn+'= '        ;insert name and =.
        apost = "'"             ;quote a quote
        if N_elements(value) NE 1 then $
                message,'Keyword Value (third parameter) must be a scalar'

        case stype of         ;which type?

7:      begin
          upval = strupcase(value)      ;force upper case.
          if (upval eq 'T') || (upval eq 'F') then begin
                strput,h,upval,29  ;insert logical value.
            end else begin              ;other string?
                if keyword_set(save_as_string) then $
                  svalue = strtrim(value,2) else svalue = value
                if strlen(svalue) gt 18 then begin       ;long string
                    strput, h, apost + strmid(svalue,0,68) + apost + $
                        ' /' + ncomment,10
                    header[i] = h
                    return
                endif
                strput, h, apost + svalue,10       ;insert string val
                strput, h, apost, 11 + (strlen(svalue)>8)  ;pad string vals
          endelse                                          ;to at least 8 chars
          endcase

5:      BEGIN
        IF (N_ELEMENTS(format) EQ 1) THEN $             ; use format keyword
            v = string(value, FORMAT='('+strupcase(format)+')') $
        ELSE v = STRING(value, FORMAT='(G19.12)')
        s = strlen(v)                                   ; right justify
        strput, h, v, (30-s)>10
        END               

 else:  begin
        if ~save_as_null then begin
        if stype EQ 1 then begin
             if !VERSION.RELEASE GE '8.4' && ISA(value,/boolean) then begin
                upval = ['F','T']
                strput,h,upval[value],29 
                break
             endif else v = strtrim(fix(value),2) 
        endif else begin
        if (N_elements(format) eq 1) then $            ;use format keyword
            v = string(value, FORMAT='('+strupcase(format)+')' ) else $
            v = strtrim(strupcase(value),2)      
                                      ;convert to string, default format
        endelse                              
        s = strlen(v)                 ;right justify
        strput,h,v,(30-s)>10          ;insert
        endif
        end
 endcase

 if (~save_as_null) || (strlen(strtrim(comment)) GT 0) then begin
   strput,h,' /',30       ;add ' /'
   strput, h, ncomment, 32 ;add comment
 endif  
   header[i] = h          ;save line
 
 return
 end
 
 pro writefits, filename, data, header, heap, Append = Append, Silent = Silent, $
       compress = compress, CheckSum = checksum, NaNValue = NaNvalue
       
;+
; NAME:
;       WRITEFITS
; PURPOSE:
;       Write IDL array and header variables to a disk FITS file.    
;
; EXPLANATION:
;       A minimal FITS header is created if not supplied.
;       WRITEFITS works for all types of FITS files except random groups
;
; CALLING SEQUENCE:
;       WRITEFITS, filename, data [, header, /APPEND, /COMPRESS, /CHECKSUM] 
;
; INPUTS:
;       FILENAME = String containing the name of the file to be written.
;
;       DATA = Image array to be written to FITS file.    If DATA is 
;              undefined or a scalar, then only the FITS header (which
;              must have NAXIS = 0) will be written to disk
;
; OPTIONAL INPUT:
;       HEADER = String array containing the header for the FITS file.
;                If variable HEADER is not given, the program will generate
;                a minimal FITS header.
;       HEAP -   A byte array giving the heap area following, e.g. a variable
;                length binary table
;
; OPTIONAL INPUT KEYWORD:
;       /APPEND - If this keyword is set then the supplied header and data
;                array are assumed to be an extension and are appended onto
;                the end of an existing FITS file.    If the file does not 
;                exist, then WRITEFITS will create one with a minimal primary
;                header (and /EXTEND keyword) and then append the supplied
;                extension header and array.     Note that the primary
;                header in an existing file must already have an EXTEND
;                keyword to indicate the presence of an FITS extension.
;       /COMPRESS - If this keyword is set, then the FITS file is written as
;                a gzip compressed file.   An extension '.gz' is appended to
;                to the file name if it does not already exist.   The /COMPRESS
;                option is incompatible with the /APPEND option.
;      /Checksum - If set, then the CHECKSUM keywords to monitor data integrity
;                 will be included in the FITS header.    For more info, see
;                 http://fits.gsfc.nasa.gov/registry/checksum.html
;                 By default, checksum keywords will updated if they are already
;                 in the FITS header.
;       NaNvalue - Value in the data array which represents missing pixels.
;		 This keyword should only used when missing pixels are not
;		 represented by NaN values in the input array.
; OUTPUTS:
;       None
;
; RESTRICTIONS:
;       (1) It recommended that BSCALE and BZERO not be used (or set equal
;           to 1. and 0) except with integer data
;       (2) WRITEFITS will remove any group parameters from the FITS header
;       (3) As of Feb 2008, WRITEFITS no longer requires the primary header of a
;           FITS file with extensions to contain the EXTEND keyword, consistent 
;           with Section 4.4.2.1 of the FITS 3.0 standard.    A warning is still 
;           given.  See http://fits.gsfc.nasa.gov/fits_standard.html
;
; EXAMPLE:
;       Write a randomn 50 x 50 array as a FITS file creating a minimal header.
;
;       IDL> im = randomn(seed, 50, 50)        ;Create array
;       IDL> writefits, 'test', im             ;Write to a FITS file "test"
;
; PROCEDURES USED:
;       CHECK_FITS, FITS_ADD_CHECKSUM, MKHDR, MRD_HREAD, SXDELPAR, SXADDPAR, 
;       SXPAR()
;
; MODIFICATION HISTORY:
;       WRITTEN, Jim Wofford, January, 29 1989
;       Added call to IS_IEEE_BIG()  W. Landsman  Apr 96
;       Make sure SIMPLE is written in first line of header  W. Landsman Jun 97
;       Use SYSTIME() instead of !STIME    W. Landsman  July 97
;       Create a default image extension header if needed W. Landsman June 98
;       Write unsigned data types W. Landsman       December 1999
;       Update for IDL V5.3, add /COMPRESS keyword W. Landsman  February 2000
;       Correct BZERO value for unsigned data  W. Landsman   July 2000
;       Eliminate duplication of input array if possible W. Landsman April 2001
;       Use FILE_SEARCH for V5.5 or later     W. Landsman    April 2002
;       Create the file if not already present and /APPEND is set
;                                             W. Landsman    September 2002
;       Proper call to MRD_HREAD if /APPEND is set  W. Landsman December 2002 
;       Added /CHECKSUM keyword              W. Landsman     December 2002
;	    Restored NANvalue keyword, William Thompson,	     October 2003
;       Write BZERO in beginning of header for unsigned integers WL April 2004
;       Added ability to write heap array       WL             October 2004
;       Correct checksum if writing heap array   WL           November 2004
;       Assume since V5.5, no VMS support, use file_search() WL   September 2006
;       Set nbytes variable to LONG64 for very large files WL  May 2007
;       Update CHECKSUM keywords if already present  WL   Oct 2007
;       EXTEND keyword no longer required in FITS files with extensions WL Feb 2008
;       Bug fix when filename ends with '.gz' and COMPRESS is used,
;            the output file must be compressed          S. Koposov June 2008
;       Introduce V6.0 notation                W.L. Nov. 2010 
;       Set /APPEND if XTENSION specifies a table   W.L.  July 2012
;       Bug fix when /CHECKSUM used with unsigned data  W.L. June 2013
;       June 2013 bug fix introduced problem when NAXIS=0  W.L. July 2013
;       Added /Silent keyword W.L. April 2016
;		Support unsigned 64 bit data type  W.L.  January 2018
;       Fix case of header with no data and checksum  W.L.   August 2018
;-

  compile_opt idl2  

  if N_params() LT 2 then begin 
       print,'Syntax - WRITEFITS, filename, data,[ header, /APPEND, /CHECKSUM]'
       return
  endif
  
  Catch, theError
  IF theError NE 0 then begin
	Catch,/Cancel
	void = cgErrorMsg(/quiet)
	RETURN
  ENDIF

; Get information about data

  siz = size( data )      
  naxis = siz[0]                    ;Number of dimensions
  if naxis GT 0 then nax = siz[ 1:naxis ]              ;Vector of dimensions
  lim = siz[ naxis+2 ]              ;Total number of data points
  type = siz[naxis + 1]             ;Data type

;Create a primary or image extension header if not supplied by the user

        if N_elements(header) LT 2 then begin 
                if keyword_set(append) then mkhdr, header, data, /IMAGE  $
                                       else mkhdr, header, data, /EXTEND
        endif else if naxis GT 0 then $         
              check_FITS, data, header, /UPDATE, Silent= silent

  hdr = header     ;Don't modify supplied header
  
;If header indicates a table extension then set the append keyword  
  if ~keyword_set( APPEND) && ( strmid(hdr[0],0,8) EQ 'XTENSION' ) then begin
	 xten = strtrim(sxpar(hdr,'XTENSION'),2)
	 if (xten EQ 'TABLE') || (xten Eq 'BINTABLE') || (xten Eq 'A3DTABLE') $
	     then begin 
	     append = 1
	     message,'Writing FITS table extension',/INF,NoPrint = silent
	 endif    
   endif	     

  if ~keyword_set( APPEND) then begin 
         simple = 'SIMPLE  =                    T / Written by IDL:  ' $
                        + systime()  
         hdr[0] =  simple + string( replicate(32b,80-strlen(simple) ) )
         sxdelpar, hdr, [ 'GCOUNT', 'GROUPS', 'PCOUNT', 'PSIZE' ]     ;Remove random groups keywords
  endif
  
; If necessary,convert unsigned to signed.    Do not destroy the original data

  unsigned = 0
  if naxis NE 0 then begin
              
        unsigned = (type EQ 12) || (type EQ 13) || (type EQ 15)
        if  unsigned then begin
             if type EQ 12 then begin
                     sxaddpar,hdr,'BZERO',32768,' Data is Unsigned Integer', $
                              before = 'DATE'
                     newdata = fix(data - 32768)
             endif else if type EQ 13 then begin 
                    sxaddpar,hdr,'BZERO',2147483648,' Data is Unsigned Long', $
                              before = 'DATE'
                    newdata = long(data - 2147483648)
             endif else if type EQ 15 then begin
                offset = ulong64(2)^63
             	sxaddpar,hdr,'BZERO',offset,' Data is 64 bit Unsigned Long', $
             			before = 'DATE'
             	newdata = long64(data - offset )
             endif	
        endif


; For floating or double precision test for NaN values to write

  NaNtest = keyword_set(NaNvalue) && ( (type EQ 4) || (type EQ 5) )
  if NaNtest then begin
     NaNpts = where( data EQ NaNvalue, N_NaN)
     if (N_NaN GT 0) then begin
         if type EQ 4 then data[NaNpts]  = !Values.F_NaN	$
     else if type EQ 8 then data[NaNpts] = !Values.D_NaN
     endif
  endif 
  endif

; Open file and write header information

        if keyword_set( APPEND) then begin
            if (strmid( hdr[0],0,8 ) NE 'XTENSION') then begin
                   message, $
            'ERROR - "XTENSION" must be first keyword in header extension',/CON
                  return
            endif
            if ~file_test(filename)  then  begin       ;Create default primary header
                 mkhdr,h0,0b,/exten
                 writefits,filename,0b,h0, checksum = checksum
                 openu, unit, filename, /GET_LUN, /swap_if_little_endian
             endif else begin
            openu, unit, filename, /GET_LUN, /swap_if_little_endian
            mrd_hread, unit, hprimary
            extend = where( strcmp(hprimary,'EXTEND  ',8), Nextend)
            if Nextend EQ 0 then $
               message,'WARNING - EXTEND keyword not found in primary FITS header',/CON,NoPrint=silent
            endelse
                   
            file = fstat(unit)
            nbytes  = file.size
            point_lun, unit, nbytes
            npad = nbytes mod 2880
            if npad NE 0 then writeu, unit, replicate(32b, 2880 - npad)

    endif else begin

        ext = ''
        if keyword_set(COMPRESS) then begin 
            if strlowcase(strmid(filename,2,3,/reverse)) NE '.gz' $
               then ext = '.gz' 
        endif else compress = 0


       openw, unit, filename + ext, /GET_LUN, /swap_if_little_endian, $
                             compress = compress

    endelse

; Determine if an END line occurs, and add one if necessary

       endline = where( strcmp(hdr, 'END     ', 8), Nend)
     if Nend EQ 0 then begin

     message,'WARNING - An END statement has been appended to the FITS header',/INF,NoPrint=silent
     hdr = [ hdr, 'END' + string( replicate(32b,77) ) ]
     endline = N_elements(hdr) - 1 

   endif

; Add any CHECKSUM keywords if desired or already present
   
    do_Checksum = keyword_set(checksum)
    if ~do_checksum then test = sxpar(hdr,'CHECKSUM',count=do_checksum)
     if do_checksum then begin 
               if unsigned then begin 
	       if N_elements(heap) GT 0 then $
	         FITS_ADD_CheckSum, hdr, [newdata,heap] else $
		 FITS_Add_CheckSum, hdr, newdata
	       endif else begin 	 
               if N_elements(heap) GT 0 then $
	         FITS_ADD_CHECKSUM, hdr, [data,heap] else $
                 FITS_ADD_CHECKSUM, hdr, data
	       endelse	 
               endline = where( strcmp(hdr,'END     ',8), Nend)
       endif
       nmax = endline[0] + 1

; Convert to byte and force into 80 character lines

       bhdr = replicate(32b, 80l*nmax)
       for n = 0l, endline[0] do bhdr[80*n] = byte( hdr[n] )
       npad = 80l*nmax mod 2880
       writeu, unit, bhdr
       if npad GT 0 then writeu, unit,  replicate(32b, 2880 - npad)

; Write data
       if naxis EQ 0 then goto, DONE
        bitpix = sxpar( hdr, 'BITPIX' )
        nbytes = long64( N_elements( data)) * (abs(bitpix) / 8 )
        npad = nbytes mod 2880

        if unsigned then writeu, unit, newdata $
                    else writeu, unit, data 

; Write optional heap area
        if N_elements(heap) GT 0 then begin
              theap = sxpar(hdr,'THEAP', Count=N_Theap)
              if N_Theap GT 0 then begin
                  offset = theap - nbytes
                  if offset GT 0 then begin
                      writeu, unit, bytarr(offset)
                      npad = (npad + offset) mod 2880
                  endif
                  writeu, unit, heap
                  npad = (npad + N_elements(heap)) mod 2880
              endif
         endif

; ASCII Tables padded with blanks (32b) otherwise pad with zeros
        if keyword_set( APPEND) then begin
             exten = sxpar( header, 'XTENSION')
             padnum =  exten EQ 'TABLE   ' ? 32b : 0b
        endif else padnum = 0b
         
        if npad GT 0 then writeu, unit, replicate( padnum, 2880 - npad)
DONE:
        free_lun, unit  

  return
  end
  
  pro get_date, dte, in_date, OLD = old, TIMETAG = timetag
;+
; NAME:
;       GET_DATE
; PURPOSE:
;       Return the (current) UTC date in CCYY-MM-DD format for FITS headers
; EXPLANATION:
;       This is the format required by the DATE and DATE-OBS keywords in a 
;       FITS header.  
;
; CALLING SEQUENCE:
;       GET_DATE, FITS_date, [ in_date, /OLD, /TIMETAG ]
; OPTIONAL INPUTS:
;       in_date - string (scalar or vector) containing dates in IDL
;            systime() format (e.g. 'Tue Sep 25 14:56:14 2001')
; OUTPUTS:
;       FITS_date = A scalar character string giving the current date.    Actual
;               appearance of dte depends on which keywords are supplied.
;       
;       No Keywords supplied - dte is a 10 character string with the format
;               CCYY-MM-DD where <CCYY> represents a calendar year, <MM> the
;               ordinal number of a calendar month within the calendar year, 
;               and <DD> the ordinal number of a day within the calendar month.
;       /TIMETAG set - dte is a 19 character string with the format
;               CCYY-MM-DDThh:mm:ss where <hh> represents the hour in the day,
;                <mm> the minutes, <ss> the seconds, and the literal 'T' the 
;               ISO 8601 time designator
;       /OLD set - dte is an 8 character string in DD/MM/YY format
;
; INPUT KEYWORDS:
;       /TIMETAG - Specify the time to the nearest second in the DATE format
;       /OLD - Return the DATE format formerly (pre-1997) recommended for FITS
;               Note that this format is now deprecated because it uses only
;               a 2 digit representation of the year. 
; EXAMPLE:
;       Add the current date to the DATE keyword in a FITS header,h
;     
;       IDL> GET_DATE,dte
;       IDL> sxaddpar, h, 'DATE', dte, 'Date header was created'
;
; NOTES:
;       (1) A discussion of the DATExxx syntax in FITS headers can be found in
;       http://www.cv.nrao.edu/fits/documents/standards/year2000.txt
;
;       (2) Those who wish to use need further flexibility in their date 
;       formats (e.g. to use TAI time) should look at Bill Thompson's time
;       routines in http://sohowww.nascom.nasa.gov/solarsoft/gen/idl/time
;
; PROCEDURES USED:
;       DAYCNV - Convert Julian date to Gregorian calendar date
; REVISION HISTORY:
;       Written      W. Landsman          March 1991
;       Major rewrite to write new DATExxx syntax  W. Landsman  August 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Work after year 2000 even with /OLD keyword W. Landsman January 2000
;       Don't need to worry about TIME_DIFF since V5.4 W. Landsman July 2001
;       Assume since V5.4, remove LOCAL_DIFF keyword  W. Landsman April 2006
;-
 On_error,2
 compile_opt idl2
 
 if N_params() LT 1 then begin
     print,'Syntax - Get_date, FITS_date, [ in_date, /TIMETAG, /OLD ]'
     print,'  FITS_date - output string giving date(s) in FITS format'
     print,'  in-date - Optional input string giving date in systime() format'
     return
 endif

 if N_elements(in_date) GT 0 then begin
     mn = strmid(in_date,4,3)
     month = month_cnv(mn)
     day = fix(strmid(in_date,8,2))
     ihr = fix(strmid(in_date,11,2))
     imn = fix(strmid(in_date,14,2))
     sec = fix(strmid(in_date,17,2))
     yr = fix(strmid(in_date,20,4))
 endif else begin
     seconds = systime(1)          ;Number of seconds since Jan 1, 1970
     dayseconds = 86400.D0               ;Number of seconds in a day
     mjd = seconds/dayseconds + 40587.0D
     jd =  2400000.5D + mjd
     DAYCNV, jd, yr, month, day, hr
 endelse 

 if keyword_set(old) then begin

 if yr GE 2000 then yr = yr - 100
 dte =  string(day,f='(I2.2)') + '/' + string(month,f='(i2.2)') +  $
        '/' + string( yr-1900,f='(I2.2)')                          

 endif else $ 

 dte =  string(yr,f='(I4.4)') + '-' + string(month,f='(i2.2)') + '-' + $
        string(day,f='(I2.2)')

 if keyword_set(TIMETAG) then begin
 if N_elements(in_date) EQ 0 then begin
   ihr = fix(hr)
   mn = (hr - ihr)*60.
   imn = fix(mn)
   sec = round((mn - imn)*60.)
 endif 

 dte =  dte + 'T' + string(ihr,f='(I2.2)') + ':' + string(imn,f='(I2.2)') +  $
               ':' + string(round(sec),f='(I2.2)')
 endif
         
 return
 end
 
 
 PRO DAYCNV, XJD, YR, MN, DAY, HR
;+
; NAME:
;       DAYCNV
; PURPOSE:
;       Converts Julian dates to Gregorian calendar dates
;
; EXPLANATION:
;       Duplicates the functionality of the intrinsic JUL2GREG procedure
;       which was introduced in V8.2.1
; CALLING SEQUENCE:
;       DAYCNV, XJD, YR, MN, DAY, HR
;
; INPUTS:
;       XJD = Julian date, positive double precision scalar or vector
;
; OUTPUTS:
;       YR = Year (Integer)
;       MN = Month (Integer)
;       DAY = Day (Integer)
;       HR = Hours and fractional hours (Real).   If XJD is a vector,
;               then YR,MN,DAY and HR will be vectors of the same length.
;
; EXAMPLE:
;       IDL> DAYCNV, 2440000.D, yr, mn, day, hr    
;
;       yields yr = 1968, mn =5, day = 23, hr =12.   
;
; WARNING:
;       Be sure that the Julian date is specified as double precision to
;       maintain accuracy at the fractional hour level.
;
; METHOD:
;       Uses the algorithm of Fliegel and Van Flandern (1968) as reported in
;       the "Explanatory Supplement to the Astronomical Almanac" (1992), p. 604
;       Works for all Gregorian calendar dates with XJD > 0, i.e., dates after
;       -4713 November 23.
; REVISION HISTORY:
;       Converted to IDL from Yeoman's Comet Ephemeris Generator, 
;       B. Pfarr, STX, 6/16/88
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2
 compile_opt idl2

 if N_params() lt 2 then begin
    print,"Syntax - DAYCNV, xjd, yr, mn, day, hr'
    print,'  Julian date, xjd, should be specified in double precision'
    return
 endif

; Adjustment needed because Julian day starts at noon, calendar day at midnight

 jd = long(xjd)                         ;Truncate to integral day
 frac = double(xjd) - jd + 0.5          ;Fractional part of calendar day
 after_noon = where(frac ge 1.0, Next)
 if Next GT 0 then begin                ;Is it really the next calendar day?
      frac[after_noon] = frac[after_noon] - 1.0
      jd[after_noon] = jd[after_noon] + 1
 endif
 hr = frac*24.0
 l = jd + 68569
 n = 4*l / 146097l
 l = l - (146097*n + 3l) / 4
 yr = 4000*(l+1) / 1461001
 l = l - 1461*yr / 4 + 31        ;1461 = 365.25 * 4
 mn = 80*l / 2447
 day = l - 2447*mn / 80
 l = mn/11
 mn = mn + 2 - 12*l
 yr = 100*(n-49) + yr + l
 return
 end
 
 pro check_FITS, im, hdr, dimen, idltype, UPDATE = update, NOTYPE = notype, $
                SDAS = sdas, FITS = fits, SILENT = silent, ERRMSG = errmsg, $
                ALLOW_DEGEN=ALLOW_DEGEN
;+
; NAME:
;       CHECK_FITS
; PURPOSE:
;       Check that keywords in a FITS header array match the associated data  
; EXPLANATION:
;       Given a FITS array IM, and a associated FITS header HDR, this
;       procedure will check that
;               (1) HDR is a string array, and IM is defined and numeric   
;               (2) The NAXISi values in HDR are appropriate to the dimensions 
;                   of IM
;               (3) The BITPIX value in HDR is appropriate to the datatype of IM
;       If the /UPDATE keyword is present, then the FITS header will be 
;       modified, if necessary, to force agreement with the image array
;
; CALLING SEQUENCE:
;       check_FITS, im, hdr, [ dimen, idltype, /UPDATE, /NOTYPE, /SILENT
;                              ERRMSG = ]'
;
; INPUT PARAMETERS:
;       IM -  FITS array, e.g. as read by READFITS
;       HDR - FITS header (string array) associated with IM
;
; OPTIONAL OUTPUTS:
;       dimen - vector containing actual array dimensions
;       idltype- data type of the FITS array as specified in the IDL SIZE
;               function (1 for BYTE, 2 for 16 bit integer, 3 for 32 bit integer etc.)
;
; OPTIONAL KEYWORD INPUTS:
;       /NOTYPE - If this keyword is set, then only agreement of the array
;               dimensions with the FITS header are checked, and not the 
;               data type.
;       /UPDATE - If this keyword is set then the BITPIX, NAXIS and NAXISi
;               FITS keywords will be updated to agree with the array
;       /FITS, /SDAS -  these are obsolete keywords that now do nothing 
;       /SILENT - If keyword is set and nonzero, the informational messages 
;               will not be printed
;       /ALLOW_DEGEN - Don't check for degenerate axes.
; OPTIONAL KEYWORD OUTPUT:
;       ERRMSG  = If this keyword is present, then any error messages will be
;                 returned to the user in this parameter rather than
;                 depending on the MESSAGE routine in IDL.  If no errors are
;                 encountered, then a null string is returned.  
;
; PROCEDURE:
;       Program checks the NAXIS and NAXISi keywords in the header to
;       see if they match the image array dimensions, and checks whether
;       the BITPIX keyword agrees with the array type.
;
; PROCEDURE CALLS:
;       FXADDPAR, FXPAR(), SXDELPAR
; MODIFICATION HISTORY:
;       Written, December 1991  W. Landsman Hughes/STX to replace CHKIMHD
;       No error returned if NAXIS=0 and IM is a scalar   W. Landsman  Feb 93
;       Fixed bug for REAL*8 STSDAS data W. Landsman July 93
;       Make sure NAXIS agrees with NAXISi  W. Landsman  October 93
;        Converted to IDL V5.0   W. Landsman   September 1997
;       Allow unsigned data types   W. Landsman December 1999
;       Allow BZERO = 0 for unsigned data types   W. Landsman January 2000
;       Added ERRMSG keyword, W. Landsman February 2000
;       Use FXADDPAR to put NAXISi in proper order   W. Landsman August 2000
;       Improper FXADDPAR call for DATATYPE keyword  W. Landsman December 2000
;       Remove explicit setting of obsolete !err W. Landsman February 2004
;       Remove SDAS support   W. Landsman       November 2006
;       Fix dimension errors introduced Nov 2006
;       Work again for null arrays W. Landsman/E. Hivon May 2007
;       Use V6.0 notation  W.L.  Feb. 2011
;       Add /ALLOW_DEGEN, William Thompson, 26-Jun-2019
;- 
 compile_opt idl2
 On_error,2

 if N_params() LT 2 then begin
    print,'Syntax - CHECK_FITS, im, hdr, dimen, idltype, '
    print,'            [ /UPDATE, /NOTYPE, ERRMSG=, /SILENT ]'
    return
 endif

 if arg_present(errmsg) then errmsg = ''       

 if size(hdr,/TNAME) NE 'STRING' then begin        ;Is hdr of string type?
        message= 'FITS header is not a string array'
        if  N_elements(ERRMSG) GT 0 then errmsg = message else $
             message, 'ERROR - ' + message, /CON
             return 
 endif

 im_info = size(im,/struc)
 ndimen = im_info.n_dimensions
 if ndimen GT 0 then dimen = im_info.dimensions[0:ndimen-1]
 idltype = im_info.type

 
 nax = fxpar( hdr, 'NAXIS', Count = N_naxis ) 
 if N_naxis EQ 0 then begin
        message = 'FITS header missing NAXIS keyword'
        if  N_elements(errmsg) GT 0 then errmsg = message else $
             message,'ERROR - ' + message,/CON 
             return 
 endif
        
 if ndimen EQ 0  then $             ;Null primary array
     if nax EQ 0 then return else begin
         message = 'FITS array is not defined'
         if  N_elements(errmsg) GT 0 then errmsg = message else $
             message,'ERROR - ' +message,/con 
             return 
     endelse

 
 naxis = fxpar( hdr, 'NAXIS*')
 naxi = N_elements( naxis )
 if nax GT naxi then begin                 ;Does NAXIS agree with # of NAXISi?
        if keyword_set( UPDATE) then begin
                fxaddpar, hdr, 'NAXIS', naxi
                if ~keyword_set(SILENT) then message, /INF, $
        'NAXIS changed from ' + strtrim(nax,2) + ' to ' + strtrim(naxi,2)
        endif else begin 
                message =  'FITS header has NAXIS = ' + strtrim(nax,2) + $
                ', but only ' + strtrim(naxi, 2) + ' axes defined'
                if  N_elements(ERRMSG) GT 0 then errmsg = message else $
                    message, 'ERROR - ' + message
                return
        endelse
 endif

 last = naxi-1                  ;Remove degenerate dimensions
 if ~keyword_set(allow_degen) then begin
     while ( (naxis[last] EQ 1) && (last GE 1) ) do last--
     if last NE nax-1 then begin
         naxis = naxis[ 0:last]
     endif
 endif

 if ( ndimen NE last + 1 ) then begin
    if ~keyword_set( UPDATE) THEN begin
        message = $
        '# of NAXISi keywords does not match # of array dimensions'
        if  N_elements(ERRMSG) GT 0 then errmsg = message else $
                                     message,'ERROR - ' + message,/CON 
        return 
 
     endif else goto, DIMEN_ERROR
 endif

 for i = 0,last do begin
      if naxis[i] NE dimen[i] then begin
      if ~keyword_set( UPDATE ) then begin
          message =  'Invalid NAXIS' + strtrim( i+1,2 ) + $
	             ' keyword value in header'
          if  N_elements(ERRMSG) GT 0 then errmsg = message else $ 
                                       message,'ERROR - ' + message,/CON
          return 
      endif else goto, DIMEN_ERROR
    endif
 endfor

BITPIX:     

 if ~keyword_set( NOTYPE ) then begin

 
  bitpix = fxpar( hdr, 'BITPIX')
  
    case idltype of

     1: if bitpix NE 8 then goto, BITPIX_ERROR
     2: if bitpix NE 16 then goto, BITPIX_ERROR  
     4: if bitpix NE -32 then goto, BITPIX_ERROR       
     3: if bitpix NE 32 then goto, BITPIX_ERROR 
     5: if bitpix NE -64 then goto, BITPIX_ERROR 
     12:if bitpix NE 16 then goto, BITPIX_ERROR
     13: if bitpix NE 32 then goto, BITPIX_ERROR
     
     else: begin
              message = 'Data array is not a valid FITS datatype'
             if  N_elements(ERRMSG) GT 0 then errmsg = message else $
                                          message,'ERROR - ' + message,/CON
             return 
      end

   endcase

 endif

 return

BITPIX_ERROR:
    if keyword_set( UPDATE ) then begin
    bpix = [0, 8, 16, 32, -32, -64, 32, 0, 0, 0, 0, 0, 16,32 ]
    comm = ['',' Character or unsigned binary integer', $
               ' 16-bit twos complement binary integer', $
               ' 32-bit twos complement binary integer', $
               ' IEEE single precision floating point', $
               ' IEEE double precision floating point', $
               ' 32-bit twos complement binary integer','','','','','', $
               ' 16-bit unsigned binary integer', $
               ' 32-bit unsigned binary integer' ]
    bitpix = bpix[idltype]
    comment = comm[idltype]
    if ~keyword_set(SILENT) then message, /INF, $
        'BITPIX value of ' + strtrim(bitpix,2) +  ' added to FITS header'
    fxaddpar, hdr, 'BITPIX', bitpix, comment
    return

  endif else begin 
       message = 'BITPIX value of ' + strtrim(bitpix,2) + $
                 ' in FITS header does not match array'
      if  N_elements(ERRMSG) GT 0  then errmsg = message else  $
          message,'ERROR - ' + message,/CON
      return
 endelse

DIMEN_ERROR:
   if keyword_set( UPDATE ) then begin
        fxaddpar, hdr, 'NAXIS', ndimen, before = 'NAXIS1'
	naxis = 'NAXIS' + strtrim(indgen(ndimen)+1,2)
        for i = 1, ndimen do fxaddpar, hdr, naxis[i-1], dimen[i-1], $
                'Number of positions along axis ' + strtrim(i,2), $
                after = 'NAXIS' + strtrim(i-1,2)          
        if naxi GT ndimen then begin
                for i = ndimen+1, naxi do sxdelpar, hdr, 'NAXIS'+strtrim(i,2)
        endif
        if ~keyword_set(SILENT) then message, /INF, $
                'NAXIS keywords in FITS header have been updated'
        goto, BITPIX
   endif

 end
 
         FUNCTION FXPAR, HDR, NAME, ABORT, COUNT=MATCHES, COMMENT=COMMENTS, $
                        START=START, PRECHECK=PRECHECK, POSTCHECK=POSTCHECK, $
                        NOCONTINUE = NOCONTINUE, DATATYPE=DATATYPE, $
                        NULL=K_NULL, NAN=NAN, MISSING=MISSING, $
                        MULTIVALUE=MULTIVALUE
;+
; NAME: 
;        FXPAR()
; PURPOSE: 
;       Obtain the value of a parameter in a FITS header.
; EXPLANATION: 
;       The first 8 chacters of each element of HDR are searched for a match to
;       NAME.  If the keyword is one of those allowed to take multiple values
;       ("HISTORY", "COMMENT", or "        " (blank)), then the value is taken
;       as the next 72 characters.  Otherwise, it is assumed that the next
;       character is "=", and the value (and optional comment) is then parsed
;       from the last 71 characters.  An error occurs if there is no parameter
;       with the given name.
;      
;       If the value is too long for one line, it may be continued on to the
;       the next input card, using the CONTINUE Long String Keyword convention.
;       For more info, http://fits.gsfc.nasa.gov/registry/continue_keyword.html
;       
;
;       Complex numbers are recognized as two numbers separated by one or more
;       space characters.
;
;       If a numeric value has no decimal point (or E or D) it is returned as
;       type LONG.  If it contains more than 8 numerals, or contains the
;       character 'D', then it is returned as type DOUBLE.  Otherwise it is
;       returned as type FLOAT.    If an integer is too large to be stored as
;       type LONG, then it is returned as DOUBLE.
;
;       If a keyword is in the header and has no value, then the default
;       missing value is returned as explained below.  This can be
;       distinguished from the case where the keyword is not found by the fact
;       that COUNT=0 in that case, while existing keywords without a value will
;       be returned with COUNT=1 or more.
;
; CALLING SEQUENCE: 
;       Result = FXPAR( HDR, NAME  [, ABORT, COUNT=, COMMENT=, /NOCONTINUE ] )
;
;       Result = FXPAR(HEADER,'DATE')           ;Finds the value of DATE
;       Result = FXPAR(HEADER,'NAXIS*')         ;Returns array dimensions as
;                                               ;vector
; REQUIRED INPUTS: 
;       HDR     = FITS header string array (e.g. as returned by FXREAD).  Each
;                 element should have a length of 80 characters
;       NAME    = String name of the parameter to return.  If NAME is of the
;                 form 'keyword*' then an array is returned containing values
;                 of keywordN where N is an integer.  The value of keywordN
;                 will be placed in RESULT(N-1).  The data type of RESULT will
;                 be the type of the first valid match of keywordN
;                 found, unless DATATYPE is given.
; OPTIONAL INPUT: 
;       ABORT   = String specifying that FXPAR should do a RETALL if a
;                 parameter is not found.  ABORT should contain a string to be
;                 printed if the keyword parameter is not found.  If not
;                 supplied, FXPAR will return with a negative !err if a keyword
;                 is not found.
; OUTPUT: 
;       The returned value of the function is the value(s) associated with the
;       requested keyword in the header array.
;
;       If the parameter is complex, double precision, floating point, long or
;       string, then the result is of that type.  Apostrophes are stripped from
;       strings.  If the parameter is logical, 1 is returned for T, and 0 is
;       returned for F.
;
;       If NAME was of form 'keyword*' then a vector of values are returned.
;
; OPTIONAL INPUT KEYWORDS: 
;       DATATYPE = A scalar value, indicating the type of vector
;                  data.  All keywords will be cast to this type.
;                  Default: based on first keyword.
;                  Example: DATATYPE=0.0D (cast data to double precision)
;       START   = A best-guess starting position of the sought-after
;                 keyword in the header.  If specified, then FXPAR
;                 first searches for scalar keywords in the header in
;                 the index range bounded by START-PRECHECK and
;                 START+POSTCHECK.  This can speed up keyword searches
;                 in large headers.  If the keyword is not found, then
;                 FXPAR searches the entire header.  
;
;                 If not specified then the entire header is searched.
;                 Searches of the form 'keyword*' also search the
;                 entire header and ignore START.
;
;                 Upon return START is changed to be the position of
;                 the newly found keyword.  Thus the best way to
;                 search for a series of keywords is to search for
;                 them in the order they appear in the header like
;                 this:
;
;                       START = 0L
;                       P1 = FXPAR('P1', START=START)
;                       P2 = FXPAR('P2', START=START)
;
;       PRECHECK = If START is specified, then PRECHECK is the number
;                  of keywords preceding START to be searched.
;                  Default: 5
;       POSTCHECK = If START is specified, then POSTCHECK is the number
;                   of keywords after START to be searched.
;                   Default: 20
;       /NOCONTINUE = If set, then continuation lines will not be read, even
;                 if present in the header
;       MISSING = By default, this routine returns 0 when keyword values are
;                 not found.  This can be overridden by using the MISSING
;                 keyword, e.g. MISSING=-1.
;       /NAN    = If set, then return Not-a-Number (!values.f_nan) for missing
;                 values.  Ignored if keyword MISSING is present.
;       /NULL   = If set, then return !NULL (undefined) for missing values.
;                 Ignored if MISSING or /NAN is present, or if earlier than IDL
;                 version 8.0.  If multiple values would be returned, then
;                 MISSING= or /NAN should be used instead of /NULL, making sure
;                 that the datatype is consistent with the non-missing values,
;                 e.g. MISSING='' for strings, MISSING=-1 for integers, or
;                 MISSING=-1.0 or /NAN for floating point.  /NAN should not be
;                 used if the datatype would otherwise be integer.
;       /MULTIVALUE = Allow multiple values to be returned, if found in the
;                     header.
; OPTIONAL OUTPUT KEYWORD:
;       COUNT   = Optional keyword to return a value equal to the number of
;                 parameters found by FXPAR.
;       COMMENTS= Array of comments associated with the returned values.
;
; PROCEDURE CALLS: 
;       GETTOK(), VALID_NUM
; SIDE EFFECTS: 
;
;       The system variable !err is set to -1 if parameter not found, 0 for a
;       scalar value returned.  If a vector is returned it is set to the number
;       of keyword matches found.    This use of !ERR is deprecated.
;
;       If a keyword occurs more than once in a header, a warning is given,
;       and the first occurence is used.  However, if the keyword is "HISTORY",
;       "COMMENT", or "        " (blank), then multiple values are returned.
;
; NOTES:
;	The functions SXPAR() and FXPAR() are nearly identical, although
;	FXPAR() has slightly more sophisticated parsing.   There is no
;	particular reason for having two nearly identical procedures, but
;	both are too widely used to drop either one.
;
; REVISION HISTORY: 
;       Version 1, William Thompson, GSFC, 12 April 1993.
;               Adapted from SXPAR
;       Version 2, William Thompson, GSFC, 14 October 1994
;               Modified to use VALID_NUM instead of STRNUMBER.  Inserted
;               additional call to VALID_NUM to trap cases where character
;               strings did not contain quotation marks.
;       Version 3, William Thompson, GSFC, 22 December 1994
;               Fixed bug with blank keywords, following suggestion by Wayne
;               Landsman.
;       Version 4, Mons Morrison, LMSAL, 9-Jan-98
;               Made non-trailing ' for string tag just be a warning (not
;               a fatal error).  It was needed because "sxaddpar" had an
;               error which did not write tags properly for long strings
;               (over 68 characters)
;       Version 5, Wayne Landsman GSFC, 29 May 1998
;               Fixed potential problem with overflow of LONG values
;       Version 6, Craig Markwardt, GSFC, 28 Jan 1998, 
;               Added CONTINUE parsing         
;       Version 7, Craig Markwardt, GSFC, 18 Nov 1999,
;               Added START, PRE/POSTCHECK keywords for better
;               performance
;       Version 8, Craig Markwardt, GSFC, 08 Oct 2003,
;               Added DATATYPE keyword to cast vector keywords type
;       Version 9, Paul Hick, 22 Oct 2003, Corrected bug (NHEADER-1)
;       Version 10, W. Landsman, GSFC  2 May 2012
;               Keywords of form "name_0" could confuse vector extractions
;       Version 11 W. Landsman, GSFC 24 Apr 2014
;               Don't convert LONG64 numbers to to double precision
;       Version 12, William Thompson, 13-Aug-2014
;               Add keywords MISSING, /NAN, and /NULL
;	Version 13, W. Landsman 25-Jan-2018
;		Return ULONG64 integer if LONG64 would overflow
;       Version 14, William Thompson, 03-Jun-2019
;               Add /MULTIVALUE keyword
;       Version 15, Mats LÃ¶fdahl, 11-Sep-2019
;               Read CONTINUE mechanism multi-line comments.
;-
;------------------------------------------------------------------------------
;
;  Check the number of parameters.
;
        IF N_PARAMS() LT 2 THEN BEGIN
            PRINT,'Syntax:  result =  FXPAR( HDR, NAME  [, ABORT ])'
            RETURN, -1
        ENDIF
;
;  Determine the default value for missing data.
;
        CASE 1 OF 
            N_ELEMENTS(MISSING) EQ 1: MISSING_VALUE = MISSING
            KEYWORD_SET(NAN): MISSING_VALUE = !VALUES.F_NAN
            KEYWORD_SET(K_NULL) AND !VERSION.RELEASE GE '8.': $
              DUMMY = EXECUTE('MISSING_VALUE = !NULL')
            ELSE: MISSING_VALUE = 0
        ENDCASE
        VALUE = MISSING_VALUE
;
;  Determine the abort condition.
;
        IF N_PARAMS() LE 2 THEN BEGIN
            ABORT_RETURN = 0
            ABORT = 'FITS Header'
        END ELSE ABORT_RETURN = 1
        IF ABORT_RETURN THEN ON_ERROR,1 ELSE ON_ERROR,2
;
;  Check for valid header.  Check header for proper attributes.
;
        S = SIZE(HDR)
        IF ( S[0] NE 1 ) OR ( S[2] NE 7 ) THEN $
            MESSAGE,'FITS Header (first parameter) must be a string array'
;
;  Convert the selected keyword NAME to uppercase.
;
        NAM = STRTRIM( STRUPCASE(NAME) )
;
;  Determine if NAME is of form 'keyword*'.  If so, then strip off the '*', and
;  set the VECTOR flag.  One must consider the possibility that NAM is an empty
;  string.
;
        NAMELENGTH1 = (STRLEN(NAM) - 1) > 1
        IF STRPOS( NAM, '*' ) EQ NAMELENGTH1 THEN BEGIN    
            NAM = STRMID( NAM, 0, NAMELENGTH1)  
            VECTOR = 1                          ;Flag for vector output  
            NAME_LENGTH = STRLEN(NAM)           ;Length of name 
            NUM_LENGTH = 8 - NAME_LENGTH        ;Max length of number portion  
            IF NUM_LENGTH LE 0 THEN MESSAGE,    $
                'Keyword length must be 8 characters or less'
;
;  Otherwise, extend NAME with blanks to eight characters.
;
        ENDIF ELSE BEGIN
            WHILE STRLEN(NAM) LT 8 DO NAM = NAM + ' '
            VECTOR = 0
        ENDELSE
;
;  If of the form 'keyword*', then find all instances of 'keyword' followed by
;  a number.  Store the positions of the located keywords in NFOUND, and the
;  value of the number field in NUMBER.
;
        IF N_ELEMENTS(START)     EQ 0 THEN START = -1L
        START = LONG(START[0])
        IF NOT VECTOR AND START GE 0 THEN BEGIN
            IF N_ELEMENTS(PRECHECK)  EQ 0 THEN PRECHECK = 5
            IF N_ELEMENTS(POSTCHECK) EQ 0 THEN POSTCHECK = 20
            NHEADER = N_ELEMENTS(HDR)
            MN = (START - PRECHECK)  > 0
            MX = (START + POSTCHECK) < (NHEADER-1)      ;Corrected bug
            KEYWORD = STRMID(HDR[MN:MX], 0, 8)
        ENDIF ELSE BEGIN
            RESTART:
            START   = -1L
            KEYWORD = STRMID( HDR, 0, 8)
        ENDELSE

        IF VECTOR THEN BEGIN
            NFOUND = WHERE(STRPOS(KEYWORD,NAM) GE 0, MATCHES)
            IF ( MATCHES GT 0 ) THEN BEGIN
                NUMST= STRMID(HDR[NFOUND], NAME_LENGTH, NUM_LENGTH)
                NUMBER = INTARR(MATCHES)-1
                FOR I = 0, MATCHES-1 DO         $
                    IF VALID_NUM( NUMST[I], NUM) THEN NUMBER[I] = NUM
                IGOOD = WHERE(NUMBER GE 0, MATCHES)
                IF MATCHES GT 0 THEN BEGIN
                    NFOUND = NFOUND[IGOOD]
                    NUMBER = NUMBER[IGOOD]
 		    G = WHERE(NUMBER GT 0, MATCHES)
 		    IF MATCHES GT 0 THEN NUMBER = NUMBER[G]     
		ENDIF
            ENDIF
;
;  Otherwise, find all the instances of the requested keyword.  If more than
;  one is found, and NAME is not one of the special cases, then print an error
;  message.
;
        ENDIF ELSE BEGIN
            NFOUND = WHERE(KEYWORD EQ NAM, MATCHES)
            IF MATCHES EQ 0 AND START GE 0 THEN GOTO, RESTART
            IF START GE 0 THEN NFOUND = NFOUND + MN
            IF (MATCHES GT 1) THEN BEGIN
                IF KEYWORD_SET(MULTIVALUE) THEN BEGIN
                    VECTOR = 1
                    NUMBER = INDGEN(MATCHES) + 1
                END ELSE IF (NAM NE 'HISTORY ') AND (NAM NE 'COMMENT ') $
                  AND (NAM NE '') THEN BEGIN
                    MESSAGE, /INFORMATIONAL, 'WARNING- Keyword ' +   $
                             NAM + 'located more than once in ' + ABORT
                ENDIF
            ENDIF
            IF (MATCHES GT 0) THEN START = NFOUND[MATCHES-1]
        ENDELSE
;
;  Extract the parameter field from the specified header lines.  If one of the
;  special cases, then done.
;
        IF MATCHES GT 0 THEN BEGIN
            VALUE = MISSING_VALUE
            LINE = HDR[NFOUND]
            SVALUE = STRTRIM( STRMID(LINE,9,71),2)
            IF (NAM EQ 'HISTORY ') OR (NAM EQ 'COMMENT ') OR    $
                    (NAM EQ '        ') THEN BEGIN
                VALUE = STRTRIM( STRMID(LINE,8,72),2)
                COMMENTS = STRARR(N_ELEMENTS(VALUE))
;
;  Otherwise, test to see if the parameter contains a string, signalled by
;  beginning with a single quote character (') (apostrophe).
;
            END ELSE FOR I = 0,MATCHES-1 DO BEGIN
                IF ( STRMID(SVALUE[I],0,1) EQ "'" ) THEN BEGIN
                    TEST = STRMID( SVALUE[I],1,STRLEN( SVALUE[I] )-1)
                    NEXT_CHAR = 0
                    OFF = 0
                    VALUE = ''
                    COMMENT = ''
;
;  Find the next apostrophe.
;
NEXT_APOST:
                    ENDAP = STRPOS(TEST, "'", NEXT_CHAR)
                    IF ENDAP LT 0 THEN MESSAGE,         $
                        'WARNING: Value of '+NAME+' invalid in '+ABORT+ " (no trailing ')", /info
                    VALUE = VALUE + STRMID( TEST, NEXT_CHAR, ENDAP-NEXT_CHAR )
;
;  Test to see if the next character is also an apostrophe.  If so, then the
;  string isn't completed yet.  Apostrophes in the text string are signalled as
;  two apostrophes in a row.
;
                    IF STRMID( TEST, ENDAP+1, 1) EQ "'" THEN BEGIN    
                        VALUE = VALUE + "'"
                        NEXT_CHAR = ENDAP+2      
                        GOTO, NEXT_APOST
                    ENDIF
;
;  Extract the comment, if any.
;
                    SLASH = STRPOS(TEST, "/", ENDAP)
                    IF SLASH GE 0 THEN COMMENT += STRMID(TEST, SLASH+1, STRLEN(TEST)-SLASH-1)
;
; CM 19 Sep 1997
; This is a string that could be continued on the next line.  Check this
; possibility with the following four criteria: *1) Ends with '&'
; (2) Next line is CONTINUE  (3) LONGSTRN keyword is present (recursive call to
;  FXPAR) 4. /NOCONTINUE is not set

                    IF NOT KEYWORD_SET(NOCONTINUE) THEN BEGIN
                       OFF = OFF + 1
                       VAL = STRTRIM(VALUE,2)

                       IF (STRLEN(VAL) GT 0) AND $
                      (STRMID(VAL, STRLEN(VAL)-1, 1) EQ '&') AND $
                      (STRMID(HDR[NFOUND[I]+OFF],0,8) EQ 'CONTINUE') THEN BEGIN
                       IF (SIZE(FXPAR(HDR, 'LONGSTRN',/NOCONTINUE)))[1] EQ 7 THEN BEGIN                    
                      VALUE = STRMID(VAL, 0, STRLEN(VAL)-1)
                      TEST = HDR[NFOUND[I]+OFF]
                      TEST = STRMID(TEST, 8, STRLEN(TEST)-8)
                      TEST = STRTRIM(TEST, 2)
                      IF STRMID(TEST, 0, 1) NE "'" THEN MESSAGE, $
                        'ERROR: Invalidly CONTINUEd string in '+ABORT
                      NEXT_CHAR = 1
                      GOTO, NEXT_APOST
                    ENDIF
                   ENDIF
    ENDIF

;
;  If not a string, then separate the parameter field from the comment field.
;  If there is no value field, then use the default "missing" value.
;
                ENDIF ELSE BEGIN
                    VALUE = MISSING_VALUE
                    TEST = SVALUE[I]
                    IF TEST EQ '' THEN BEGIN
                        COMMENT = ''
                        GOTO, GOT_VALUE
                    ENDIF
                    SLASH = STRPOS(TEST, "/")
                    IF SLASH GE 0 THEN BEGIN
                        COMMENT = STRMID(TEST, SLASH+1, STRLEN(TEST)-SLASH-1)
                        IF SLASH GT 0 THEN TEST = STRMID(TEST, 0, SLASH) ELSE $
                            GOTO, GOT_VALUE
                    END ELSE COMMENT = ''
;
;  Find the first word in TEST.  Is it a logical value ('T' or 'F')?
;
                    TEST2 = TEST
                    VALUE = GETTOK(TEST2,' ')
                    TEST2 = STRTRIM(TEST2,2)
                    IF ( VALUE EQ 'T' ) THEN BEGIN
                        VALUE = 1
                    END ELSE IF ( VALUE EQ 'F' ) THEN BEGIN
                        VALUE = 0
                    END ELSE BEGIN
;
;  Test to see if a complex number.  It's a complex number if the value and the
;  next word, if any, both are valid numbers.
;
                        IF STRLEN(TEST2) EQ 0 THEN GOTO, NOT_COMPLEX
                        VALUE2 = GETTOK(TEST2,' ')
                        IF VALID_NUM(VALUE,VAL1) AND VALID_NUM(VALUE2,VAL2) $
                                THEN BEGIN
                            VALUE = COMPLEX(VAL1,VAL2)
                            GOTO, GOT_VALUE
                        ENDIF
;
;  Not a complex number.  Decide if it is a floating point, double precision,
;  or integer number.  If an error occurs, then a string value is returned.
;  If the integer is not within the range of a valid long value, then it will 
;  be converted to a double.  
;
NOT_COMPLEX:
                        ON_IOERROR, GOT_VALUE
                        VALUE = TEST
                        IF NOT VALID_NUM(VALUE) THEN GOTO, GOT_VALUE
                        IF (STRPOS(VALUE,'.') GE 0) OR (STRPOS(VALUE,'E') $
                                GE 0) OR (STRPOS(VALUE,'D') GE 0) THEN BEGIN
                            IF ( STRPOS(VALUE,'D') GT 0 ) OR $
                                    ( STRLEN(VALUE) GE 8 ) THEN BEGIN
                            	VALUE = DOUBLE(VALUE)
                                END ELSE VALUE = FLOAT(VALUE)
                        ENDIF ELSE BEGIN
                            LMAX = 2.0D^31 - 1.0D
                            LMIN = -2.0D^31       ;Typo fixed Feb 2010
                            IF STRMID(VALUE,0,1) NE '-' THEN BEGIN
                            	VALUE = ULONG64(VALUE)
                            	IF VALUE LT ULONG64(2)^63-1 THEN VALUE = LONG64(VALUE)
                            ENDIF ELSE VALUE = LONG64(VALUE)
                            if (VALUE GE LMIN) and (VALUE LE LMAX) THEN $
                                VALUE = LONG(VALUE)
                        ENDELSE
                            
;
GOT_VALUE:
                        ON_IOERROR, NULL
                    ENDELSE
                ENDELSE         ; if string
;
;  Add to vector if required.
;
                IF VECTOR THEN BEGIN
                    MAXNUM = MAX(NUMBER)
                    IF ( I EQ 0 ) THEN BEGIN
                        IF N_ELEMENTS(DATATYPE) EQ 0 THEN BEGIN
                            ;; Data type determined from keyword
                            SZ_VALUE = SIZE(VALUE)
                        ENDIF ELSE BEGIN
                            ;; Data type requested by user
                            SZ_VALUE = SIZE(DATATYPE[0])
                        ENDELSE
                        RESULT = MAKE_ARRAY( MAXNUM, TYPE=SZ_VALUE[1])
                        COMMENTS = STRARR(MAXNUM)
                    ENDIF 
                    RESULT[   NUMBER[I]-1 ] =  VALUE
                    COMMENTS[ NUMBER[I]-1 ] =  COMMENT
                ENDIF ELSE BEGIN
                    COMMENTS = COMMENT
                ENDELSE
            ENDFOR
;
;  Set the value of !ERR for the number of matches for vectors, or simply 0
;  otherwise.
;
            IF VECTOR THEN BEGIN
                !ERR = MATCHES
                RETURN, RESULT
            ENDIF ELSE !ERR = 0
;
;  Error point for keyword not found.
;
        ENDIF ELSE BEGIN
            IF ABORT_RETURN THEN MESSAGE,'Keyword '+NAM+' not found in '+ABORT
            !ERR = -1
        ENDELSE
;
        RETURN, VALUE
        END
        
        function gettok,st,char, exact=exact, notrim=notrim
;+
; NAME:
;	GETTOK                                    
; PURPOSE:
;	Retrieve the first part of a (vector) string up to a specified character
; EXPLANATION:
;	GET TOKen - Retrieve first part of string until the character char 
;	is encountered.   
;
; CALLING SEQUENCE:
;	token = gettok( st, char, [ /EXACT, /NOTRIM ] )
;
; INPUT:
;	char - character separating tokens, scalar string
;
; INPUT-OUTPUT:
;	st - string to get token from (on output token is removed unless
;            /NOTRIM is set), scalar or vector
;
; OUTPUT:
;	token - extracted string value is returned, same dimensions as st
; OPTIONAL INPUT KEYWORD:
;       /EXACT -  The default behaviour of GETTOK is to remove any leading 
;              blanks and (if the token is a blank) convert tabs to blanks.    
;              Set the /EXACT keyword to skip these steps and leave the 
;              input string unchanged before searching for the  character 
;              tokens. 
;
;      /NOTRIM - if set, then the input string is left unaltered 
; EXAMPLE:
;	If ST is ['abc=999','x=3.4234'] then gettok(ST,'=') would return
;	['abc','x'] and ST would be left as ['999','3.4234'] 
;
; PROCEDURE CALLS:
;       REPCHR()
; HISTORY
;	version 1  by D. Lindler APR,86
;	Remove leading blanks    W. Landsman (from JKF)    Aug. 1991
;       V5.3 version, accept vector input   W. Landsman February 2000
;       Slightly faster implementation  W. Landsman   February 2001
;       Added EXACT keyword  W. Landsman March 2004
;       Assume since V5.4, Use COMPLEMENT keyword to WHERE W. Landsman Apr 2006
;       Added NOTRIM keyword W. L. March 2011
;-
;----------------------------------------------------------------------
  On_error,2                           ;Return to caller
  compile_opt idl2

   if N_params() LT 2 then begin
       print,'Syntax - token = gettok( st, char, [ /EXACT, /NOTRIM] )'
       return,-1
   endif

; if char is a blank treat tabs as blanks

 if ~keyword_set(exact) then begin
    st = strtrim(st,1)              ;Remove leading blanks and tabs
    if char EQ ' ' then begin 
       tab = string(9b)                 
       if max(strpos(st,tab)) GE 0 then st = repchr(st,tab,' ')
    endif
  endif
  token = st

; find character in string

  pos = strpos(st,char)
  test = pos EQ -1
  bad = where(test, Nbad, Complement = good, Ncomplement=Ngood)
  if Nbad GT 0 && ~keyword_set(notrim) then st[bad] = ''
 
; extract token
 if Ngood GT 0 then begin
    stg = st[good]
    pos = reform( pos[good], 1, Ngood )
    token[good] = strmid(stg,0,pos)
    if ~keyword_set(notrim) then st[good] = strmid(stg,pos+1)
 endif

;  Return the result.

 return,token
 end
 
 ;+
; NAME: 
;     VALID_NUM()
; PURPOSE:               
;     Check if a string is a valid number representation.
; EXPLANATION:              
;     The input string is parsed for characters that may possibly
;     form a valid number.  It is more robust than simply checking
;     for an IDL conversion error because that allows strings such
;     as '22.3qwert' to be returned as the valid number 22.3
;
;     This function had a major rewrite in August 2008 to use STREGEX
;     and allow vector input.    It should be backwards compatible.
; CALLING SEQUENCE: 
;     IDL> status = valid_num(string  [,value]  [,/integer])
;    
; INPUTS:
;     string  -  the string to be tested, scalar or array
;               
; RETURNS
;     status - byte scalar or array, same size as the input string
;              set to 1 where the string is a  valid number, 0 for invalid
; OPTIONAL OUTPUT:               
;     value     - The value the string decodes to, same size as input string.
;           This will be returned as a double precision number unless 
;           /INTEGER is present, in which case a long integer is returned.
;           
; OPTIONAL INPUT KEYWORD:          
;    /INTEGER   -  if present code checks specifically for an integer.
; EXAMPLES:
;     (1) IDL> print,valid_num(3.2,/integer) 
;        --> 0     ;Since 3.2 is not an integer 
;     (2) IDL> str =['-0.03','2.3g', '3.2e12']
;         IDL> test = valid_num(str,val)
;              test = [1,0,1]    &  val =  [-0.030000000 ,NaN ,3.2000000e+12]
; REVISION HISTORY:
;          Version 1, C D Pike, RAL, 24-May-93
;          Version 2, William Thompson, GSFC, 14 October 1994
;                       Added optional output parameter VALUE to allow
;                       VALID_NUM to replace STRNUMBER in FITS routines.
;          Version 3 Wayne Landsman rewrite to use STREGEX, vectorize
;          Version 4 W.L. (fix from C. Markwardt) Better Stregex expression, 
;                    was missing numbers like '134.' before Jan 1 2010
;-            

FUNCTION valid_num, string, value, INTEGER=integer
 On_error,2
 compile_opt idl2 
 
; A derivation of the regular expressions below can be found on 
; http://wiki.tcl.tk/989

   if keyword_set(INTEGER) then $ 
    st = '^[-+]?[0-9][0-9]*$'  else $                    ;Integer
     st = '^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eEdD][-+]?[0-9]+)?$' ;F.P.
   
;Simple return if we just need a boolean test.
    if N_params() EQ 1 then return, stregex(strtrim(string,2),st,/boolean)

   
      vv = stregex(strtrim(string,2),st,/boolean)      
      if size(string,/N_dimen) EQ 0 then begin     ;Scalar
         if vv then $
            value= keyword_set(integer) ? long(string) : double(string) 
      endif else begin                             ;Array 
         
      g = where(vv,Ng)
      if Ng GT 0 then begin      ;Need to create output vector
        if keyword_set(integer) then begin 
              value = vv*0L 
              value[g] = long(string[g])
        endif else begin 
                value = replicate(!VALUES.D_NAN,N_elements(vv))
                value[g] = double(string[g])
        endelse 
        endif   
        endelse 
     
       return,vv
      end
      
      pro sxdelpar, h, parname
;+
; NAME:
;	SXDELPAR
; PURPOSE:
;	Procedure to delete a keyword parameter(s) from a FITS header
;
; CALLING SEQUENCE:
;	sxdelpar, h, parname
;
; INPUTS:
;	h - FITS header, string array
;	parname - string or string array of keyword name(s) to delete
;
; OUTPUTS:
;	h - updated FITS header, If all lines are deleted from 
;		the header, then h is returned with a value of 0
;
; EXAMPLE:
;	Delete the astrometry keywords CDn_n from a FITS header, h
;
;	IDL> sxdelpar, h, ['CD1_1','CD1_2','CD2_1','CD2_2']
;
; NOTES:
;	(1)  No message is returned if the keyword to be deleted is not found
;	(2)  All appearances of a keyword in the header will be deleted
; HISTORY:
;   version 1  D. Lindler Feb. 1987
;	Test for case where all keywords are deleted    W. Landsman Aug 1995 
;   Allow for headers with more than 32767 lines W. Landsman Jan. 2003
;   Use ARRAY_EQUAL, cleaner syntax  W. L.  July 2009
;------------------------------------------------------------------
 On_error,2
 compile_opt idl2

 if N_Params() LT 2 then begin
      print,'Syntax - SXDELPAR, h, parname'
      return
 endif

; convert parameters to string array of upper case names of length 8 char


 if size(parname,/type) NE 7 then $
         message,'Keyword name(s) must be a string or string array'
 par = strtrim( strupcase(parname),2 ) 

 sz = size(h,/structure)
 if (sz.N_dimensions NE 1) || (sz.type NE 7) then $
	message,'FITS header (1st parameter) must be a string array'

 nlines = sz.N_elements		;number of lines in header array
 pos = 0L		;position in compressed header with keywords removed

; loop on header lines

 keyword = strtrim( strmid(h,0,8), 2 )
 for i = 0L, nlines-1 do begin
        if array_equal(keyword[i] NE par, 1b) then begin   
 	   h[pos] = h[i]		;keep it
	   pos++		;increment number of lines kept
	   if keyword[i] eq 'END' then break  	;end of header
        endif 
 endfor

 if pos GT 0 then h = h[0:pos-1] else h = 0	      ;truncate

 return
 end
 
 function SXPAR, hdr, name, abort, COUNT=matches, COMMENT = comments, $
                IFound = number, NoContinue = NoContinue, SILENT = silent, $
                DUP = dup, NULL = K_Null, NAN = NaN, MISSING = Missing
;+
; NAME:
;      SXPAR
; PURPOSE:
;      Obtain the value of a parameter in a FITS header
;
; CALLING SEQUENCE:
;      result = SXPAR( Hdr, Name, [ Abort, COUNT=, COMMENT =, /NoCONTINUE, 
;                                        DUP=,   /SILENT  ])   
;
; INPUTS:
;      Hdr =  FITS header array, (e.g. as returned by READFITS) 
;             string array, each element should have a length of 80 characters      
;
;      Name = String name of the parameter to return.   If Name is of the
;             form 'keyword*' then an array is returned containing values of
;             keywordN where N is a positive (non-zero) integer.  The value of 
;             keywordN will be placed in RESULT[N-1].  The data type of RESULT 
;             will be the type of the first valid match of keywordN found.
;
; OPTIONAL INPUTS:
;       ABORT - string specifying that SXPAR should do a RETALL
;               if a parameter is not found.  ABORT should contain
;               a string to be printed if the keyword parameter is not found.
;               The default string is 'FITS Header' if abort is a integer or 
;               single character.     If ABORT is not supplied, SXPAR will return
;               quietly with COUNT = 0  if a keyword is not found.
;
; OPTIONAL INPUT KEYWORDS: 
;       DUP =  When the FITS keyword exists more than once in a header, set 
;                 DUP to a positive integer to specify which value is to be 
;                 read.   For example, set DUP = 2 to read the value of the 
;                 second appearance of a keyword in the header.    If DUP is not
;                 supplied, then by default, SXPAR() reads the *last* appearance and
;                 issues a warning.
;       MISSING = By default, this routine returns 0 when keyword values are
;                 not found.  This can be overridden by using the MISSING
;                 keyword, e.g. MISSING=-1.
;       /NAN    = If set, then return Not-a-Number (!values.f_nan) for missing
;                 values.  Ignored if keyword MISSING is present.
;       /NOCONTINUE = If set, then continuation lines will not be read, even
;                 if present in the header
;       /NULL   = If set, then return !NULL (undefined) for missing values.
;                 Ignored if MISSING or /NAN is present, or if earlier than IDL
;                 version 8.0.  If multiple values would be returned, then
;                 MISSING= or /NAN should be used instead of /NULL, making sure
;                 that the datatype is consistent with the non-missing values,
;                 e.g. MISSING='' for strings, MISSING=-1 for integers, or
;                 MISSING=-1.0 or /NAN for floating point.  /NAN should not be
;                 used if the datatype would otherwise be integer.
;       /SILENT - Set this keyword to suppress warning messages about duplicate
;                 keywords in the FITS header.
;
; OPTIONAL OUTPUT KEYWORDS:
;       COUNT - Optional keyword to return a value equal to the number of 
;               parameters found by SXPAR, integer scalar
;
;       COMMENT - Array of comments associated with the returned values
;       IFOUND - Array of found keyword indicies when Name is of the form keyword*
;              For example, one searches for 'TUNIT*' and the FITS header contains
;              TUNIT1, TUNIT2, TUNIT4, and TUNIT6 then IFOUND woud be returned as
;              [1,2,4,6].    Set to zero if Name is not of the form keyword*.
;
; OUTPUTS:
;       Function value = value of parameter in header.
;               If parameter is double precision, floating, long or string,
;               the result is of that type.  Apostrophes are stripped
;               from strings.  If the parameter is logical, 1b is
;               returned for T, and 0b is returned for F.
;               If Name was of form 'keyword*' then a vector of values
;               are returned.
;
; SIDE EFFECTS:
;       !ERR is set to -1 if parameter not found, 0 for a scalar
;       value returned.  If a vector is returned !ERR is set to the
;       number of keyword matches found.    The use of !ERR is deprecated, and
;       instead the COUNT keyword is preferred
;
;       If a keyword (except HISTORY or COMMENT) occurs more than once in a 
;       header, and the DUP keyword is not supplied, then a warning is given, 
;       and the *last* occurrence is used.
;
; EXAMPLES:
;       Given a FITS header, h, return the values of all the NAXISi values
;       into a vector.    Then place the history records into a string vector.
;
;       IDL> naxisi = sxpar( h ,'NAXIS*')         ; Extract NAXISi value
;       IDL> history = sxpar( h, 'HISTORY' )      ; Extract HISTORY records
;
; PROCEDURE:
;       The first 8 chacters of each element of Hdr are searched for a 
;       match to Name.  The value from the last 20 characters is returned.  
;       An error occurs if there is no parameter with the given name.
;
;       If a numeric value has no decimal point it is returned as type
;       LONG.   If it has a decimal point, and contains more than 8 numerals, or 
;       contains the character 'D', then it is returned as type DOUBLE.  Otherwise
;       it is returned as type FLOAT.    Very large integer values, outside
;       the range of valid LONG, are returned as LONG64.
;
;       If the value is too long for one line, it may be continued on to the
;       the next input card, using the CONTINUE convention.  For more info,
;       see http://fits.gsfc.nasa.gov/registry/continue_keyword.html
;
;       Complex numbers are recognized as two numbers separated by one or more
;       space characters.
;
; NOTES:
;       The functions SXPAR() and FXPAR() are nearly identical, although
;       FXPAR() has slightly more sophisticated parsing, and additional keywords
;       to specify positions in the header to search (for speed), and to force
;       the output to a specified data type..   There is no
;       particular reason for having two nearly identical functions, but
;       both are too widely used to drop either one.
;
; PROCEDURES CALLED:
;       cgErrorMsg(), GETTOK(), VALID_NUM()
; MODIFICATION HISTORY:
;       DMS, May, 1983, STPAR Written.
;       D. Lindler Jan 90 added ABORT input parameter
;       J. Isensee Jul,90 added COUNT keyword
;       W. Thompson, Feb. 1992, added support for FITS complex values.
;       W. Thompson, May 1992, corrected problem with HISTORY/COMMENT/blank
;               keywords, and complex value error correction.
;       W. Landsman, November 1994, fix case where NAME is an empty string 
;       W. Landsman, March 1995,  Added COMMENT keyword, ability to read
;               values longer than 20 character
;       W. Landsman, July 1995, Removed /NOZERO from MAKE_ARRAY call
;       T. Beck May 1998, Return logical as type BYTE
;       W. Landsman May 1998, Make sure integer values are within range of LONG
;       W. Landsman Feb 1998, Recognize CONTINUE convention 
;       W. Landsman Oct 1999, Recognize numbers such as 1E-10 as floating point
;       W. Landsman Jan 2000, Only accept integer N values when name = keywordN
;       W. Landsman Dec 2001, Optional /SILENT keyword to suppress warnings
;       W. Landsman/D. Finkbeiner  Mar 2002  Make sure extracted vectors 
;             of mixed data type are returned with the highest type.
;       W.Landsman Aug 2008  Use vector form of VALID_NUM()
;       W. Landsman Jul 2009  Eliminate internal recursive call
;       W. Landsman Apr 2012  Require vector numbers be greater than 0
;       W. Landsman Apr 2014  Don't convert Long64 numbers to double
;       W. Landsman Nov 2014  Use cgErrorMsg rather than On_error,2
;       W. Landsman Dec 2014  Return Logical as IDL Boolean in IDL 8.4 or later
;       W. Landsman May 2015  Added IFound output keyword
;       J. Slavin Aug 2015 Allow for 72 character par values (fixed from 71)
;       W. Landsman Sep 2015  Added Missing, /NULL and /NaN keywords
;       W. Landsman Oct 2017   Added DUP keyword,  Needed to support distortion
;			table lookup parameters
;       W. Landsman Jan 2018 Return ULONG64 integer if LONG64 will overflow
;       W. Landsman/Y. Yang May 2018 MISSING keyword was always returning 0
;       W. Landsman/M. Knight Aug 2018 Clean up use of Abort parameter
;-
;----------------------------------------------------------------------
 compile_opt idl2

 if N_params() LT 2 then begin
     print,'Syntax -  result =  sxpar( hdr, name, [abort])'
     print,'   Input Keywords:    /NOCONTINUE, /SILENT, MISSING=, /NAN, /NULL'
     print,'   Output Keywords:   COUNT=,  COMMENT= '
     return, -1
 endif 
 
 ;
;  Determine the default value for missing data.
;
        CASE 1 OF 
            N_ELEMENTS(MISSING) EQ 1: MISSING_VALUE = MISSING
            KEYWORD_SET(NAN): MISSING_VALUE = !VALUES.F_NAN
            KEYWORD_SET(K_NULL) AND !VERSION.RELEASE GE '8.': $
              DUMMY = EXECUTE('MISSING_VALUE = !NULL')
            ELSE: MISSING_VALUE = 0
        ENDCASE
        VALUE = MISSING_VALUE
;
 
 if N_elements(abort) EQ 0 then begin
      abort_return = 0
      abort = 'FITS Header'
 endif else begin 
       if strlen(strtrim(abort,2)) LE 1 then abort = 'FITS header'
       abort_return = 1
 endelse
 
 if abort_return then On_error,1 else begin
      Catch, theError
      if theError NE 0 then begin
           Catch,/Cancel
	   void = cgErrorMsg(/quiet)
	   return,-1
	   endif
   endelse
;       Check for valid header

;Check header for proper attributes.
  if ( size(hdr,/N_dimen) NE 1 ) || ( size(hdr,/type) NE 7 ) then $
           message,'FITS Header (first parameter) must be a string array'

  nam = strtrim( strupcase(name) )      ;Copy name, make upper case     


;  Determine if NAME is of form 'keyword*'.  If so, then strip off the '*', and
;  set the VECTOR flag.  One must consider the possibility that NAM is an empty
;  string.

   namelength1 = (strlen(nam) - 1 ) > 1         
   if strpos( nam, '*' ) EQ namelength1 then begin    
            nam = strmid( nam, 0, namelength1)  
            vector = 1                  ;Flag for vector output  
            name_length = strlen(nam)   ;Length of name 
            num_length = 8 - name_length        ;Max length of number portion  
            if num_length LE 0 then  $ 
                  message, 'Keyword length must be 8 characters or less'

;  Otherwise, extend NAME with blanks to eight characters.

    endif else begin  
                while strlen(nam) LT 8 do nam += ' ' ;Make 8 chars long
                vector = 0      
    endelse


;  If of the form 'keyword*', then find all instances of 'keyword' followed by
;  a number.  Store the positions of the located keywords in IFOUND, and the
;  value of the number field in NUMBER.

        histnam = (nam eq 'HISTORY ') || (nam eq 'COMMENT ') || (nam eq '') 
        keyword = strmid( hdr, 0, 8)
	    number = 0
 
        if vector then begin
            ifound = where(strpos(keyword,nam) GE 0, matches)
            if  matches GT 0  then begin
                numst= strmid( hdr[ifound], name_length, num_length)
	        	igood = where(VALID_NUM(numst,/INTEGER), matches)

				if matches GT 0 then begin 
		     			ifound = ifound[igood]
             			number = long(numst[igood])
		     			g = where(number GT 0, matches)
 		     			if matches GT 0 then number = number[g]
				endif 
           endif

;  Otherwise, find all the instances of the requested keyword.  If more than
;  one is found, and NAME is not one of the special cases, then check if DUP keyword 
;  supplied to determine which one to use.   If DUP not supplied then issue a
;  warning and use the *last* appearance.
;  

        endif else begin
            ifound = where(keyword EQ nam, matches)
             if (matches GT 1) && ~histnam then begin
                if N_elements(dup) EQ 1 then begin
                if dup LE matches then begin 
                        ifound = ifound[dup-1] 
                        matches = 1
                endif else begin
                   message,/inf,'Warning - keyword ' + strtrim(nam,2) + $
                   ' located ' + strtrim(matches,2) + ' times in FITS header'
                   message,/inf,'But DUP specified as ' + strtrim(dup,2)
                endelse 
                endif else begin  
                if ~keyword_set(silent) then $
                message,/informational, 'Warning - keyword ' +   $
                strtrim(nam,2) + ' located more than once in ' + abort
                endelse
        endif
        endelse


; Process string parameter 

 if matches GT 0 then begin
  line = hdr[ifound]
  svalue = strtrim( strmid(line,9,71),2)
  if histnam then $
       value = strtrim(strmid(line,8,72),2) else for i = 0,matches-1 do begin
      if ( strmid(svalue[i],0,1) EQ "'" ) then begin   ;Is it a string?
                  test = strmid( svalue[i],1,strlen( svalue[i] )-1)
                  next_char = 0
                  off = 0
                  value = '' 
          NEXT_APOST:
                  endap = strpos(test, "'", next_char)      ;Ending apostrophe  
                  if endap LT 0 then $ 
                            MESSAGE,'Value of '+name+' invalid in '+abort
                  value += strmid( test, next_char, endap-next_char )  

;  Test to see if the next character is also an apostrophe.  If so, then the
;  string isn't completed yet.  Apostrophes in the text string are signalled as
;  two apostrophes in a row.

                 if strmid( test, endap+1, 1) EQ "'" then begin    
                    value += "'"
                    next_char = endap+2         
                    goto, NEXT_APOST
                 endif      

; Extract the comment, if any
                
                slash = strpos( test, "/", endap )
                if slash LT 0 then comment = '' else    $
                        comment = strmid( test, slash+1, strlen(test)-slash-1 )

; This is a string that could be continued on the next line.  Check this
; possibility with the following four criteria: *1) Ends with '&'
; (2) Next line is CONTINUE  (3) LONGSTRN keyword is present (recursive call to
; SXPAR) 4. /NOCONTINE is not set

    if ~keyword_set(nocontinue) then begin
                off++
                val = strtrim(value,2)

                if (strlen(val) gt 0) && $
                  (strmid(val, strlen(val)-1, 1) EQ '&') && $
                  (strmid(hdr[ifound[i]+off],0,8) EQ 'CONTINUE') then $
		      if ~array_equal(keyword EQ 'LONGSTRN',0b) then begin 
                  value = strmid(val, 0, strlen(val)-1)
                  test = hdr[ifound[i]+off]
                  test = strmid(test, 8, strlen(test)-8)
                  test = strtrim(test, 2)
                  if strmid(test, 0, 1) NE "'" then message, $
                    'ERROR: Invalidly CONTINUEd string in '+ abort
                  next_char = 1
                  GOTO, NEXT_APOST
                ENDIF
    ENDIF


; Process non-string value  

          endif else begin
               value = missing_value
               test = svalue[i]
               if test EQ '' then begin
                        comment = ''
                        GOTO, got_value
                endif
                slash = strpos( test, "/" )
                if slash GE 0 then begin
                        comment = strmid( test, slash+1, strlen(test)-slash-1 )
                        if slash GT 0 then test = strmid(test, 0, slash) else $
                            GOTO, got_value
                endif else comment = ''

; Find the first word in TEST.  Is it a logical value ('T' or 'F') ?

                test2 = test
                value = gettok(test2,' ')
                true = 1b
                false = 0b
                if !VERSION.RELEASE GE 8.4 then begin
                	true =  boolean(true) 
                	false = boolean(false)
               endif                      

               if ( value EQ 'T' ) then value = true else $
               if ( value EQ 'F' ) then value = false else begin

;  Test to see if a complex number.  It's  a complex number if the value and
;  the next word, if any, are both valid values.

                if strlen(test2) EQ 0 then goto, NOT_COMPLEX
                value2 = gettok( test2, ' ') 
                if value2 EQ '' then goto, NOT_COMPLEX
                On_ioerror, NOT_COMPLEX
                value2 = float(value2)
                value = complex(value,value2)
                goto, GOT_VALUE

;  Not a complex number.  Decide if it is a floating point, double precision,
;  or integer number.

NOT_COMPLEX:
                On_IOerror, GOT_VALUE
                 
                  if (strpos(value,'.') GE 0) || (strpos(value,'E') GT 0) $
                  || (strpos(value,'D') GE 0) then begin  ;Floating or double?
                      if ( strpos(value,'D') GT 0 ) || $  ;Double?
                         ( strlen(value) GE 8 ) then value = double(value) $
                                                else value = float(value)
                       endif else begin                   ;Long integer
                            lmax = 2.0d^31 - 1.0d
                            lmin = -2.0d^31      ;Typo fixed Feb 2010

                            if strmid(value,0,1) NE '-' then begin 
                            		value =ulong64(value) 
                            		if value lt ulong64(2)^63-1 then value = long64(value)
                            endif else value = long64(value)
                            if (value GE lmin) && (value LE lmax) then $
                                value = long(value) 
                       endelse

GOT_VALUE:
                On_IOerror, NULL
                endelse
             endelse; if c eq apost

;  Add to vector if required

         if vector then begin
               if ( i EQ 0 ) then begin
                     maxnum = max(number)
                     dtype = size(value,/type)
                     result = make_array( maxnum, TYPE = dtype )
                     comments = strarr( maxnum )
               endif 
               if size(value,/type) GT dtype then begin   ;Do we need to recast?
                    result = result + 0*value
                    dtype = size(value,/type)
               endif
               result[ number[i]-1 ] =  value
               comments[ number[i]-1 ] = comment
          endif else $
                comments = comment
  endfor

  if vector then begin
         !ERR = matches     
         return, result
  endif else !ERR = 0

endif  else  begin  

     if abort_return then message,/CON, $
        'Keyword '+ strtrim(nam,2) + ' not found in '+abort
     !ERR = -1
endelse     

return, value       

END                 
