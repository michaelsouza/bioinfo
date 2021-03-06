# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import time
from scipy.optimize import minimize

def print_problem(G):
    print('Problem')
    print('  nnodes: ', G['nnodes'])
    print('  nedges: ', G['nedges'])
    print('     G.I: ', G['I'])
    print('     G.J: ', G['J'])
    print('     G.L: ', G['L'])
    print('     G.U: ', G['U'])

def load_problem(fname):
    fid_problem = '%s%s' % (fname, '.csv')
    print('Reading problem : %s' % fid_problem)
    G = pd.read_csv(fid_problem)
    G = G.to_dict()
    G['nedges'] = len(G['I'].keys())
    nid_max = 0
    for k in range(G['nedges']):
        i = G['I'][k]
        j = G['J'][k]
        nid_max = max([nid_max, i, j])
        lij = G['L'][k]
        uij = G['U'][k]
        if i > j:
            print('   Changing edge index order from (%d,%d) to (%d,%d)' % (i,j,j,i))
            G['I'][k] = j
            G['J'][k] = i
        if lij > uij:
            raise Exception('Inconsistent data (lij > uij): lij = %g and uij %g' % (lij, uij))
    G['nnodes'] = nid_max
    G['I'] = np.array([G['I'][k] - 1 for k in range(G['nedges'])], dtype=int)
    G['J'] = np.array([G['J'][k] - 1 for k in range(G['nedges'])], dtype=int)
    G['L'] = np.array([G['L'][k] for k in range(G['nedges'])], dtype=float)
    G['U'] = np.array([G['U'][k] for k in range(G['nedges'])], dtype=float)
    fid_xsol = '%s%s' % (fname, '_xsol.csv')
    print('Reading solution: %s' % fid_xsol)
    xsol = pd.read_csv(fid_xsol)
    xsol = np.array([[xi.x, xi.y, xi.z] for index, xi in xsol.iterrows()],dtype=float)
    print_problem(G)
    return G, xsol

def array2matrix(x):
    if len(x.shape) == 1:
        # convert array to matrix
        return x.reshape((int(len(x) / 3), 3))
    return x

def matrix2array(x):
    if len(x.shape) > 1:
        return x.reshape((int(len(x) * 3),))
    return x
    
def xinit_jjmore(G, seed=None):
    nnodes = G['nnodes']
    nedges = G['nedges']
    x = np.zeros((nnodes,3), dtype=float)
    D = {} # sparse matrix (CSR)
    for i in range(nnodes):
        D[i] = {}
    for k in range(nedges):
        i = G['I'][k]
        j = G['J'][k]
        lij = G['L'][k]
        uij = G['U'][k]
        D[i][j] = (lij + uij) / 2.0
    if seed is not None: np.random.seed(None)
    for i in range(nnodes):
        xi = x[i]
        for j in D[i].keys():
            # random spherical coords centered on x[i]
            dij   = D[i][j]
            phi   = np.random.random_sample() * (2 * np.pi)
            theta = np.random.random_sample() * np.pi
            x[j]  = xi + dij * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    return x


def xinit_random(G):
    nnodes = G['nnodes']
    max_dij = np.max(G['L'] + G['U']) / 2.0
    x = max_dij * nnodes * np.random.random_sample((nnodes, 3))
    return x

    
def norm(x):
    s = 0
    for xi in x:
        s += xi**2
    return np.sqrt(s)


def fobj(G,x,dij_tol=0.0):
    f = 0.0
    for k in range(G['nedges']):
        lij = (1-dij_tol) * G['L'][k]
        uij = (1+dij_tol) * G['U'][k]
        xi = x[G['I'][k]]
        xj = x[G['J'][k]]
        dij = norm(xi - xj)
        f_lij = np.max([lij - dij, 0])
        f_uij = np.max([dij - uij, 0])
        f += f_lij + f_uij
    return f


def theta(tau, x, grad=False):
    val = tau**2
    for xk in x:
        val += xk**2
    val = np.sqrt(val)
    if grad:
        return val, x / val
    else:
        return val


def phi(alpha, tau, x, diff=False):
    val = alpha * x + np.sqrt((alpha * x)**2 + tau**2)
    if diff:
        return val, alpha * val / (val - alpha * x)
    else:
        return val


def fobj_smooth(alpha, tau, G, y, grad=False):
    x = array2matrix(y)
        
    # initializations
    f = 0.0
    g = np.zeros((len(x), 3))
    for k in range(G['nedges']):
        lij = G['L'][k]
        uij = G['U'][k]
        i = G['I'][k]
        j = G['J'][k]
        xi = x[i]
        xj = x[j]
        theta_ij, g_theta  = theta(tau, xi - xj, grad=True)
        phi_lij, g_phi_lij = phi(alpha, tau, lij - theta_ij, diff=True)
        phi_uij, g_phi_uij = phi(alpha, tau, theta_ij - uij, diff=True)
        f +=  phi_lij + phi_uij
        g[i] += (g_phi_uij - g_phi_lij) * g_theta
        g[j] -= g[i]
    if grad:
        return f, matrix2array(g)
    else:
        return f

        
def fobj_smooth_rot(alpha, tau, G, y, wi, wx, grad=False):
    if len(y.shape) == 1:
        # convert array to matrix
        x = y.reshape((int(len(y) / 3), 3))
    else:
        x = y

    # applying rotates
    for j in range(len(wi)):
        k = wi[j]
        a = wx[j]
        for i in range(k,G['nnodes']):
            x[i] = rotate(x[k-1],x[k-1]-x[k-2],a,x[i],diff=False)
            
    # diff rotate
    g_x = np.zeros((len(x),len(wi),3),dtype=float)
    for j in range(len(wi)):
        k = wi[j]
        a = 0
        for i in range(k,G['nnodes']):
            x[i],g_x[i][j] = rotate(x[k-1],x[k-1]-x[k-2],0,x[i],diff=True)
        
    # eval smooth function
    f = 0.0
    g = np.zeros((len(wi),), dtype=float)
    for k in range(G['nedges']):
        lij = G['L'][k]
        uij = G['U'][k]
        i = G['I'][k]
        j = G['J'][k]
        xi = x[i]
        xj = x[j]
        theta_ij, g_theta  = theta(tau, xi - xj, grad=True)
        phi_lij, g_phi_lij = phi(alpha, tau, lij - theta_ij, diff=True)
        phi_uij, g_phi_uij = phi(alpha, tau, theta_ij - uij, diff=True)
        f +=  phi_lij + phi_uij
        gi = (g_phi_uij - g_phi_lij) * g_theta
        gj = -gi
        g += np.dot(g_x[i], gi) + np.dot(g_x[j], gj)
    
    if grad:
        return f, g
    else:
        return f
        
        
def numdiff(f, x):
    fx = f(x)
    if np.isscalar(x):
        x = np.array([x],dtype=float)
    if np.isscalar(fx):
        fx = np.array([fx],dtype=float)
    g = np.zeros((len(fx),len(x)))
    d = np.zeros(len(x))
    h = 0.001
    for i in range(len(fx)):
        for j in range(len(x)):
            d[j] = h * max([np.abs(x[j]), 1.0])
            fx = f(x + d) - f(x - d)
            if np.isscalar(fx):
                g[i,j] = fx / (2 * d[j])
            else:
                g[i,j] = fx[i] / (2*d[j])
            d[j] = 0
    return g
      
    
def rotate(vec_p,vec_d,theta,vec_x):
    # -*- coding: utf-8 -*-
    """Rotates the vec_x around the axis vec_d anchored at the point vec_p.
        
    References:
        [1] http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisrotate/
    """
    (a,b,c) = vec_p # base point
    (u,v,w) = vec_d/norm(vec_d) # line direction
    (x,y,z) = vec_x # point to be rotated
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    vec_d_dot_vec_x = u*x+v*y+w*z
    omt = 1 - cos_theta
    C1 = (a*(v**2+w**2) - u*(b*v+c*w-vec_d_dot_vec_x))
    C2 = (b*(u**2+w**2) - v*(a*u+c*w-vec_d_dot_vec_x))
    C3 = (c*(u**2+v**2) - w*(a*u+b*v-vec_d_dot_vec_x))
    S1 = (-c*v + b*w - w*y + v*z)
    S2 = ( c*u - a*w + w*x - u*z)
    S3 = (-b*u + a*v - v*x + u*y)
    return np.array([C1*omt + x*cos_theta + S1*sin_theta,
                     C2*omt + y*cos_theta + S2*sin_theta,
                     C3*omt + z*cos_theta + S3*sin_theta],
                     dtype=float)

    
def rotate_diff(vec_p,vec_d,theta,vec_x):
    (a,b,c) = vec_p # base point
    (u,v,w) = vec_d/norm(vec_d) # line direction
    (x,y,z) = vec_x # point to be rotated
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    vec_d_dot_vec_x = u*x+v*y+w*z
    C1 = (a*(v**2+w**2) - u*(b*v+c*w-vec_d_dot_vec_x))
    C2 = (b*(u**2+w**2) - v*(a*u+c*w-vec_d_dot_vec_x))
    C3 = (c*(u**2+v**2) - w*(a*u+b*v-vec_d_dot_vec_x))
    S1 = (-c*v + b*w - w*y + v*z)
    S2 = ( c*u - a*w + w*x - u*z)
    S3 = (-b*u + a*v - v*x + u*y)
    return np.array([(C1-x)*sin_theta+S1*cos_theta,
                     (C2-y)*sin_theta+S2*cos_theta,
                     (C3-z)*sin_theta+S3*cos_theta],dtype=float)
    

def rotors_apply(x,wi,wx):
    y = array2matrix(x)
    y = y.copy()
    # applying rotates
    k_old=-1
    for j in range(len(wi)):
        k = wi[j]
        a = wx[j]
        if k <= k_old:
            raise Exception('The wi array must be sorted in ascending order.')
        for i in range(k,len(y)):
            y[i] = rotate(y[k-1],y[k-1]-y[k-2],a,y[i])
    return y.reshape(x.shape)
    
    
def rotors_diff(x,wi,wx):
    y = array2matrix(x)
     # diff rotate
    g = np.zeros((len(y),len(wi),3),dtype=float)
    for j in range(len(wi)):
        k = wi[j]
        a = wx[j]
        for i in range(k,len(y)):
            g[i][j] = rotate_diff(y[k-1],y[k-1]-y[k-2],a,y[i])
    return g
            
def check_theta():
    print('Checking function theta')
    tau = 0.001
    xij = np.array([1,2,3])
    dij = norm(xij)
    print('dij = ', dij)
    tij, gij = theta(tau, xij, grad=True)
    print('tij = ', tij)
    print('gij     = ', gij)
    f = lambda y: theta(tau, y, grad=False)
    gij_num = numdiff(f, xij)
    print('gij_num = ', gij_num)
    

def check_phi():
    print('Checking function phi')
    alpha = 0.5
    tau   = 0.001
    y     = np.pi
    fun   = max([y, 0])
    f,g   = phi(alpha, tau, y, diff=True)
    print('fun = ', fun)
    print('phi = ', f)
    f = lambda y: phi(alpha, tau, y, diff=False)
    print('g     = ', g)
    print('g_num = ', numdiff(f,y))
    

def check_xsol(G, x, detailed=False):
    print('Checking solution:')
    max_eij = 0.0
    for k in range(G['nedges']):
        i = G['I'][k]
        j = G['J'][k]
        lij = G['L'][k]
        uij = G['U'][k]
        xi = x[i]
        xj = x[j]
        dij = norm(xi - xj)
        f_lij = np.max([lij - dij, 0])
        f_uij = np.max([dij - uij, 0])
        eij = np.max([f_lij, f_uij])
        if eij > max_eij: max_eij = eij
        if detailed: print('   (%3d,%3d) [%3.2f,%3.2f]: %3.2f (error: %3.2E)' % (i,j,lij,uij,dij, eij))
    print('   max error: ', max_eij)
        

def check_fobj_smooth(G,xsol):
    fun = fobj(G,xsol)
    tau   = 0.001
    alpha = 0.5
    f,g = fobj_smooth(alpha, tau, G, xsol, grad=True)
    print('fun=',fun)
    print('f  =',f)
    f = lambda y: fobj_smooth(alpha, tau, G, y, grad=False)
    x = xsol.reshape((3 * len(xsol),))
    print('g    =', g)
    print('g_num=', numdiff(f,x))
    
def check_fobj_smooth_rot(G,xsol):
    print('Checking fobj_smooth_rot function')
#    y = rotate(xsol[])
#    fun = fobj(G,xsol)
#    tau   = 0.001
#    alpha = 0.5
#    wi = [4,5]
#    wx = [np.pi/8, np.pi/10]
#    f,g = fobj_smooth_rot(alpha, tau, G, xsol, wi, wx, grad=True)
#    print('fun=',fun)
#    print('f  =',f)
#    f = lambda wx: fobj_smooth_rot(alpha, tau, G, xsol, wi, wx, grad=False)
#    x = xsol.reshape((3 * len(xsol),))
#    print('g    =', g)
#    print('g_num=', numdiff(f,x))


def check_xinit(G):
    print('Generating solution using xinit_random')
    x = xinit_random(G)
    check_xsol(G, x)
    print('Generating solution using xinit_jjmore')
    x = xinit_jjmore(G)
    check_xsol(G, x)

    
def check_rotate():
    G,x = load_problem('../../instances/mdpg_1GPV_N008_EPS0.16')
    x = array2matrix(x)
    p = x[2]
    d = x[2] - x[1]
    theta = np.pi / 3
    print('y=',rotate(p,d,theta,x[3]))
    f = lambda t: rotate(p,d,t,x[3])
    print('g=',rotate_diff(p,d,theta,x[3]))
    print('g_num=',numdiff(f,theta))
    

def check_rotors_apply():
    print('Checking rotors apply')
    G,x = load_problem('../../instances/mdpg_1GPV_N008_EPS0.16')
    wi = np.array([3,4],dtype=int)
    wx = np.array([np.pi / 8, np.pi / 10], dtype=float)
    print('x=', rotors_apply(x,wi,wx))
    
    
def check_rotors_diff():    
    print('Checking rotors apply')
    G,x = load_problem('../../instances/mdpg_1GPV_N008_EPS0.16')
    wi = np.array([3,4], dtype=int)
    wx = np.array([np.pi / 8, np.pi / 10], dtype=float)
    x = matrix2array(x)
    f = lambda wx: rotors_apply(x,wi,wx)
    print('g=',rotors_diff(x,wi,wx))
    print('g_num = ', numdiff(f,wx))
    
def sph(G):
    nedges = G['nedges']
    D = (G['L'] + G['U']) / 2.0
    alpha = 0.5
    tau   = np.sort(D)[int(0.75 * nedges)]
    rho   = 0.99
    maxit = 1000
    ftol  = 1E-8
    dtol  = 1E-2
    print('SPH')
    print('   alpha = ', alpha)
    print('   tau   = ', tau)
    print('   rho   = ', rho)
    print('   maxit = ', maxit)
    print('   ftol  = ', ftol)
    print('   dtol  = ', dtol)
    
    x  = xinit_jjmore(G)
    fx = fobj(G,x)
    print('   fx    = ', fx)
    # check_xsol(G,x, detailed=True)
    for niter in range(maxit):
        tic = time.time()
        fx_old = fx
        x_old = x
        # defining and solving unconstrained problem
        f = lambda y: fobj_smooth(alpha, tau, G, y, grad=True)
        x = x.reshape((3 * len(x),)) # matrix to vector
        res = minimize(f, x, method='BFGS',jac=True)
        toc = time.time() - tic
        # update x and fx
        x = res.x
        
        x = x.reshape((int(len(x) / 3), 3)) # vector to matrix
        fx = fobj(G,x,dij_tol=dtol)
        df = (fx - fx_old) / fx_old
        dx = max(norm(x - x_old)) / max(norm(x))
        if (niter+1) % 100 == 1:
            print(' iter    fx        df       dx      nit   rho  time')
        if (niter+1) % 20 == 1:
            print('%5d %5.2E % 5.2E %5.2E %4d  %5.2E %3.2f' % (niter+1, fx, df, dx, res.nit, tau, toc))
        
        # stop criteria    
        if fx < ftol:
            break
        tau *= rho
    print('%5d %5.2E % 5.2E %5.2E %4d  %5.2E %3.2f' % (niter+1, fx, df, dx, res.nit, tau, toc))
    # check_xsol(G,x,detailed=True)
    return x

        
def sph_rot(G):
    # set the initial solution satisfying the dmdgp 
    # exact distance constraints
    x = xinit_dmdgp(G)
    w = np.zeros(len(x), dtype=float)
    nnodes = G['nnodes']
    
    # calculating derivatives
    z = x.copy()
    for k in range(4, nnodes):
        z[k],gk = rotate()
    
if __name__ == "__main__":
#    fname = '../../instances/mdpg_1GPV_N008_EPS0.16'
#    G, xsol = load_problem(fname)
    # print('xsol=', xsol)
    # check_theta()
    # check_phi()
#    check_fobj_smooth(G,xsol)
    # check_xinit(G)
#    check_rotate()
#    check_rotors_apply()
    check_rotors_diff()
#    check_fobj_smooth_rot(G,xsol)
#    sph(G)
