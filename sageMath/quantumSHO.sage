# https://doxdrum.wordpress.com/2010/02/19/playing-around-with-sagemath/

reset()
def Harm_Osc(n, x):
 return sqrt(1/(2^n*factorial(n)))/(pi)^(0.25)*e^(-x^2/2)*hermite(n,x)
p = plot(x^2/2, (x, -3.2, 3.2), color='black') + text('$x^2/2$', (3.9, 5), color='black', fontsize=13)
for n in range(5):
 p += plot(Harm_Osc(n, x)/2 + n + 0.5, (x, -5, 5), color=hue(n/5.0), fill=n + 0.5) + text('$\Psi_%s$' %n, (5.2, n+0.5), color=hue(n/5.0), fontsize=13)
show(p, axes_labels=['$x$','$E$'])
