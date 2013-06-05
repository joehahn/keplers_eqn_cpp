;keplers_eqn.pro
;by Joe Hahn, jhahn@spacescience.org, 4 June 2013.

;This IDL script reads the output generated by test_keplers_eqn.cpp,
;and makes some plots to show that the numerical errors in the call to
;keplers_eqn(), which solves Kepler's equation numerically, are within
;expected limits. All results are plotted
;against x=1-e where log10(x) is uniformly distributed over -6<log10(x)<0
;so that eccentricity e ranges over 10^(-6)<e<0.999999. The figure shows that
;keplers_eqn() has to iterate more when e is close to one. In summary,
;when N=10^7 (which is the number of times that keplers_eqn() is called)
;and max_error=1e-15, all numerical errors are<max_error as expected and
;the execution time on my 6 year old Dell Precision T3400 is 5.82 seconds.
;To execute, start IDL and enter @keplers_eqn

;read input file that was generated by test_keplers_eqn.cpp
file = 'keplers_eqn.dat'
openr, lun, file, get_lun=1
xtime_sec = read_binary(lun, data_type=5, data_dims=0)
N = read_binary(lun, data_type=3, data_dims=0)
e = read_binary(lun, data_type=5, data_dims=N)
M = read_binary(lun, data_type=5, data_dims=N)
EA = read_binary(lun, data_type=5, data_dims=N)
its = read_binary(lun, data_type=3, data_dims=N)
close, lun

;plot the error in the solution for kepler's equation versus x=1-e. Also show
;that the mean anomalies M are uniformly distributed over -Pi<M<Pi,
;and plot the number of iterations that keplers_eqn() used to solve Kepler's equation.
!p.multi = [0, 1, 3]
window, xs=550, ys=850, retain=2
x = 1 - e
err = abs(EA - e*sin(EA) - M)
plot_oo, x, err, psym=3, xrange=[1d-6, 1], yrange=[1d-18, 1d-12], xstyle=1, ystyle=1, $
    xtitle='x = 1 - e', ytitle='error', charsize=3.0
plot_oi, x, M/!dpi, psym=3, xrange=[1d-6, 1], yrange=[-1, 1], xstyle=1, ystyle=1, $
    xtitle='x = 1 - e', ytitle='mean anomaly M/!7p!3', charsize=3.0
plot_oi, x, its, psym=3, xrange=[1d-6, 1], yrange=[0, 7], xstyle=1, ystyle=1, $
    xtitle='x = 1 - e', ytitle='number of iterations', charsize=3.0
oplot, x, its, psym=6, symsize=0.4, color=128
!p.multi = 0

print, 'maximum error in solution    = ', max(err)
print, 'maximum number of iterations = ', max(its)
print, 'execution time (sec)         = ', xtime_sec

