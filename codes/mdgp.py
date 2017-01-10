# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.optimize import minimize

def load_problem(fid_problem, fid_xsol):
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
    G['nnodes'] = nid_max + 1
    G['L'] = np.array([G['L'][k] for k in range(G['nedges'])], dtype=float)
    G['U'] = np.array([G['U'][k] for k in range(G['nedges'])], dtype=float)
    print('Reading solution: %s' % fid_xsol)
    xsol = pd.read_csv(fid_xsol)
    xsol = np.array([[xi.x, xi.y, xi.z] for index, xi in xsol.iterrows()],dtype=float)
    return G, xsol

def xinit_jjmore(G):
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
    
    np.random.seed(0)
    for i in range(nnodes):
        M = []
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

def fobj(G,x):
    f = 0.0
    for k in range(G['nedges']):
        lij = G['L'][k]
        uij = G['U'][k]
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
    if len(y.shape) == 1:
        # convert array to matrix
        x = y.reshape((int(len(y) / 3), 3))
    else:
        x = y
        
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
        # print('phi_lij[%d] = % g' % (k, phi_lij))
        # print('phi_uij[%d] = % g' % (k, phi_uij))
        f +=  phi_lij + phi_uij
        g[i] += (g_phi_uij - g_phi_lij) * g_theta
        g[j] -= g[i]
    if grad:
        return f, g.reshape((3 * len(x),))
    else:
        return f

def numdiff(f, x):
    if np.isscalar(x):
        x = np.array([x],dtype=float)
    g = np.zeros(x.shape)
    fx = f(x)
    for k in range(len(x)):
        dxk  = 0.001 * max([np.abs(x[k]), 1.0])
        x[k] = x[k] + dxk
        g[k] = (f(x) - fx) / dxk
        x[k] = x[k] - dxk
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
        if detailed: print('   (%d,%d) [%3.2g,%3.2g]: %3.2g (error: %3.2g)' % (i,j,lij,uij,dij, eij))
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

def check_xinit(G):
    print('Generating solution using xinit_random')
    x = xinit_random(G)
    check_xsol(G, x)
    print('Generating solution using xinit_jjmore')
    x = xinit_jjmore(G)
    check_xsol(G, x)

def sph(G):
    nedges = G['nedges']
    D = (G['L'] + G['U']) / 2.0
    alpha = 0.5
    tau   = np.sort(D)[int(0.75 * nedges)]
    rho   = 0.99
    print('SPH')
    print('   alpha = ', alpha)
    print('   tau   = ', tau)
    print('   rho   = ', rho)
    maxit = 1000
    ftol  = 1E-8
    
    x  = xinit_jjmore(G)
    fx = fobj(G,x)
    print('   fx    = ', fx)
    check_xsol(G,x, detailed=True)
    for niter in range(maxit):
        fx_old = fx
        x_old  = x
        # defining and solving unconstrained problem
        f = lambda y: fobj_smooth(alpha, tau, G, y, grad=True)
        x = x.reshape((3 * len(x),)) # matrix to vector
        res = minimize(f, x, method='BFGS',jac=True)
        # update x and fx
        x = res.x
        
        x = x.reshape((int(len(x) / 3), 3)) # vector to matrix
        fx = fobj(G,x)
        df = (fx - fx_old) / fx_old
        dx = max(norm(x - x_old)) / max(norm(x))
        if (niter+1) % 20 == 1:
            print(' iter    fx        df       dx      nit   rho')
        print('%5d %5.2E % 5.2E %5.2E %4d  %5.2E' % (niter+1, fx, df, dx, res.nit, tau))
        # stop criteria
        if fx < ftol:
            break
        tau *= rho
    check_xsol(G,x,detailed=True)
    return x
    
G, xsol = load_problem('mdgp_08A.csv', 'mdgp_08A_xsol.csv')
print('G=\n', G)
print('xsol=', xsol)
# check_theta()
# check_phi()
# check_fobj_smooth(G,xsol)
# check_xinit(G)
sph(G)