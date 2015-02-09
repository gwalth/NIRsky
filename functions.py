
from numpy import *
from types import *
from scipy.optimize import leastsq
from scipy import interpolate

def poly_n(p, data, x):
    p = array(p)
    if type(x) == FloatType: y = 0
    else: y = zeros((x.shape[0],))
    for n in arange(p.shape[0]):
        y += p[n]*x**n
    return y - data

def gauss_fn(param, data, x):
    # From Dan's gaussfit
    """ Parameters: a0 = height of exp, a1 = center of exp, a2 = sigma
    (the width), a3 = constant term, a4 = linear term"""
    a0,a1,a2,a3,a4 = param
    return a0*exp(-(((x-a1)/a2)**2)/2) + a3 + a4*x - data

# 2D gaussian on a plane
def gauss2d_rot_fn(param,x,y):
    """amp = amplitude, mu_x = mean in x, mu_y = mean in y,
       sig_x = standard deviation in x, sig_y = standard deviation in y,
       c = constant term, c_x = linear term in x, c_y = linear term in y"""

    amp,mu_x,mu_y,sig_x,sig_y,c,cx,cy,rot = param
    #amp,mu_x,mu_y,sig_x,sig_y,c,cx,cy,cx2,cy2,rot = param

    # rotating the ellipse
    xp = (x-mu_x) * cos(rot) - (y-mu_y) * sin(rot)
    yp = (x-mu_x) * sin(rot) + (y-mu_y) * cos(rot)

    poly = c + cx*x + cy*y
    #poly = c + cx*x + cy*y + cx2*x**2 + cy2*y**2

    return amp*exp(-(((xp/sig_x)**2)/2 + ((yp/sig_y)**2)/2)) + poly

# 2D gaussian on a plane
def gauss2d_fn(param,x,y):
    """amp = amplitude, mu_x = mean in x, mu_y = mean in y,
       sig_x = standard deviation in x, sig_y = standard deviation in y,
       c = constant term, c_x = linear term in x, c_y = linear term in y"""

    amp,mu_x,mu_y,sig_x,sig_y,c,cx,cy = param
    xp = x-mu_x
    yp = y-mu_y
    poly = c + cx*x + cy*y
    return amp*exp(-(((xp/sig_x)**2)/2 + ((yp/sig_y)**2)/2)) + poly

# from http://wiki.scipy.org/Cookbook/SignalSmooth
def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

#def eparams(data,x):
#    #slope = 0
#    intercept = min(data)
#    slope = (data[-1] - data[0])/(x[-1] - x[0])
#    fdata = data - poly_n(array([intercept,slope]), 0, x)
#    return [intercept,slope]

def N_gauss_fn(param,data,x):
    # Basic N-dim gaussian function to be used with leastsq
    # Determines the number of dimensions from param
    """ Parameters: a = amplitude, xm = center of exp, sigma (the width),
        c0 = constant term, c1 = linear term"""
    # example:
    #          a1, xm1, sigma1, a2, xm2, sigma2,      c0, c1
    #   param = [0.045,95,3.5] + [0.015,80,3.5] + [0.003,0.0]

    coef  = param[:-2] # a1, xm1, sigma1, a2, xm2, sigma2,...
    c0,c1 = param[-2:]

    N = len(coef)/3

    v = 0
    for i in range(N):
        a     = coef[3*i]
        xm    = coef[3*i+1]
        sigma = coef[3*i+2]
        v += a*exp(-((x - xm)/sigma)**2/2)
    v += c0 + c1*x

    return v - data

def g_params(data,x,sigma):
    slope = 0
    intercept = min(data)
    fdata = data - poly_n(array([intercept,slope]), 0, x)
    center = argmax(fdata)
    height = fdata[center]
    width = abs(sigma*(x[1]-x[0]))
    center = x[center]
    return [height,center,width,intercept,slope]

def poly_x(x,p):
    return poly_n(p, 0, x)

def gauss_x(x,param):
    return gauss_fn(param, 0, x)

def sum2i_1(i):
    x = lambda n: 2**(n-1)

    s = 0
    for n in arange(i)+1:
        s += x(n)
        print n,x(n)
    print
    print s

def fit_gauss(y,xguess,window,sigma=2.0):

    xsize = y.shape[0]

    x0 = int(xguess-window/2.)
    x1 = int(xguess+window/2.)

    if x0 < 0: x0 = 0
    if x1 > xsize-1: x1 = xsize-1

    xfit = arange(x0,x1)
    yfit = y[x0:x1]

    param = g_params(yfit,xfit,sigma)
    sol = leastsq(gauss_fn,param,args=(yfit,xfit),full_output=1)
    return sol,[x0,x1]

def continuum_fit(spec,wav,weight,fit="spline",full_output=0,N2=25):
    """
    input:
      spec       - input spectrum
      wav        - wavelengths (length of spectrum)
      weigbts    - weights (length of spectrum) 
                     smaller values have less weight
                     larger values hvae greater weight
    optional:
      fit         - type of fit [spline,legendre], default is spline 
      full_output - return residuals, knots
    """
    if fit == "spline":
        print 
        #N2 = 25
        #N2 = 27 
        #N2 = 43
        ssize = (max(wav) - min(wav))/N2
        wav_knot = arange(min(wav),max(wav),ssize)[1:]
        C = interpolate.LSQUnivariateSpline(wav,spec,wav_knot,w=weight)
        spec_fn = C.__call__(wav)
        spec_knot = C.__call__(wav_knot)
        #print C.get_knots()
        print "Rchi2 =",C.get_residual()

        fit_dict = {"spec_knot":spec_knot,"wav_knot":wav_knot}

          
    elif fit == "legendre":
    
        Xw = (wav[-1] - wav[0])/2.0
        Xc = (wav[-1] + wav[0])/2.0
        X = (wav-Xc)/Xw
        #print Xw
        #print Xc
        A = legendre(X,order)
    
        coefficients = linalg.lstsq(A.T*weight[::,newaxis],spec*weight)[0]
        #print "Iteration=",niter
        print "Coefficients=",coefficients
        spec_fn = dot(coefficients,A)
        resid = spec - spec_fn
        rms_scatter = sqrt(sum(weight*power(resid,2))/sum(weight))
        print "RMS=",rms_scatter
 
        # probably return resid...?
        fit_dict = {}

    else:
        print "not a supported fitting method"
        sys.exit(1)


    if full_output:
        return spec_fn,fit_dict
    else:
        return spec_fn
    
