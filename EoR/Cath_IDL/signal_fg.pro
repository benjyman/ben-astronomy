pro signal_fg,nu,a,sig,pder

common wblock,ww

print,'A: ',A

; procedure to define the fitting function and its derivatives

c = 2.997e8   ; m/s
lambda = c/(nu*1.e6)
boltz = 1380.

; nu is in MHz

; [A, B, C, D, E, nu_c, delt_nu]
; [Jy, unitless, unitless, Jy, MHz, MHz]

nu_0 = 50.   ; MHz

sig1 = a[0]*(nu/nu_0)^(a[1]+a[2]*(alog(nu)-a[7]))/lambda^2*ww*2.*boltz
;sig2 = a[3]*(tanh((nu-a[4])/a[5]))*nu^2*ww
sig2 = a[3]*nu/lambda^2*ww*2.*boltz
sig3 = -a[4]*exp(-(nu-a[5])^2/2./a[6]^2)/lambda^2*ww*2.*boltz

;sechbase = (nu-a[4])/a[5]

;print,'sechbase ',sechbase

sig = sig1 + sig2 + sig3

;pder = [ [sig1/a[0]] , [alog(nu)*sig1] , [(alog(nu))*(alog(nu)-a[6])*sig1] , [nu^2*ww*tanh((nu-a[4])/a[5])] , [-nu^2*ww*a[3]/a[5]/(cosh(sechbase))^2] , [-nu^2*ww*a[3]*sechbase/a[5]/(cosh(sechbase))^2] , [alog(nu)*sig1*(-a[2])] ]

pder = [ [sig1/a[0]] , [alog(nu)*sig1] , [(alog(nu))*(alog(nu)-a[6])*sig1] , [sig2/a[3]] , [sig3/a[4]] , [sig3*(nu-a[5])/a[6]^2] , [sig3*(nu-a[5])^2/a[6]^3] , [alog(nu)*sig1*(-a[2])] ]


;print,'Signal'
;print,sig
plot,nu,sig
print,'a:'
print,a
;wait,0.1

end
