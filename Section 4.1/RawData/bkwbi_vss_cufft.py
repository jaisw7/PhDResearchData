#/usr/bin/python

"""
Author: Shashank Jaiswal

This is a Python-3 script that you can run by typing 
$ python bkwbi_vss_cufft 1.0 16 6 4 gj 4.0 64

The script needs seven inputs:
- mr mass ratio of the two species
- N number of points in each direction of the velocity mesh
- M number of points on sphere
- Nrho Number of points in the radial direction
- quadType Quadrature type used for integration
- t0 Initial time for Krook-Wu solution
- blkSize CUDA block size

Refer to the PhDThesis to understand the meaning of these parameters.

The inputs are purely meant for benchmarking. You can write a simple bash 
loop to test the effect of input on the accuracy. To run the script, 
you will need to install the following Python modules installed on your system:
- mako
- pycuda
- numpy

"""

import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda import compiler, gpuarray
from math import gamma, cos, pi
import sys
import itertools as it
from mako.template import Template
#import skcuda.fft as cufft
from cufft import (cufftPlan3d, cufftPlanMany, 
                            cufftExecD2Z, cufftExecZ2Z,
                            CUFFT_D2Z, CUFFT_Z2Z, CUFFT_FORWARD, CUFFT_INVERSE
                        )

def get_grid_for_block(block, nrow, ncol=1):
    return (int((nrow + (-nrow % block[0])) // block[0]),
            int((ncol + (-ncol % block[1])) // block[1]))

def jacobi(n, a, b, z):
    j = [1]

    if n >= 1:
        j.append(((a + b + 2)*z + a - b) / 2)
    if n >= 2:
        apb, bbmaa = a + b, b*b - a*a

        for q in range(2, n + 1):
            qapbpq, apbp2q = q*(apb + q), apb + 2*q
            apbp2qm1, apbp2qm2 = apbp2q - 1, apbp2q - 2

            aq = apbp2q*apbp2qm1/(2*qapbpq)
            bq = apbp2qm1*bbmaa/(2*qapbpq*apbp2qm2)
            cq = apbp2q*(a + q - 1)*(b + q - 1)/(qapbpq*apbp2qm2)

            # Update
            j.append((aq*z - bq)*j[-1] - cq*j[-2])

    return j


def jacobi_diff(n, a, b, z):
    dj = [0]

    if n >= 1:
        dj.extend(jp*(i + a + b + 2)/2
                  for i, jp in enumerate(jacobi(n - 1, a + 1, b + 1, z)))

    return dj

def jacobz(n, a, b):
    if(not n): return [] 

    z = np.zeros(n)
    dth = pi/(2.0*n)
    rlast = 0.0

    for k in range(n):
        r = -cos((2.0*k + 1.0)*dth);

        if k>=1: 
            r = 0.5*(r + rlast)

        for j in range(1, 5000):
            poly = jacobi(n, a, b, r)[-1];
            pder = jacobi_diff(n, a, b, r)[-1];

            tsum = 0.0
            for i in range(k): 
                tsum += 1.0/(r - z[i]);

            delr = -poly / (pder - tsum * poly);
            r   += delr;

            if( abs(delr) < 1e-20 ): break;

        z[k]  = r;
        rlast = r;
    return z


# Compute Gauss-Lobatto-Jacobi points and weights
def zwglj(Np, a, b):
    z= np.zeros(Np)
    w = np.ones(Np)*2.0

    if Np>=1:
        z[0], z[Np-1] = -1.0, 1.0;

        z[1:-1] = jacobz(Np-2, a+1.0, b+1.0); 
        w = jacobi(Np-1, a, b, z)[-1];

        fac  = pow(2.0, a+b+1.0)*gamma(a+Np)*gamma(b+Np);
        fac /= (Np-1.0)*gamma(Np)*gamma(a+b+Np+1.0);

        w = fac/(w*w)
        w[0], w[Np-1] = w[0]*(b  + 1.0), w[Np-1]*(a + 1.0);

    return z, w

# Compute Gauss-Jacobi points and weights
def zwgj(Np, a, b):
    z= jacobz(Np, a, b)
    w = jacobi_diff(Np, a, b, z)[-1]

    fac  = pow(2.0, a+b+1.0)*gamma(a+Np+1.0)*gamma(b+Np+1.0);
    fac /= gamma(Np+1.0)*gamma(a+b+Np+1.0);

    for i in range(Np): 
        w[i] = fac/(w[i]*w[i]*(1-z[i]*z[i]))
    
    return z, w

# Compute Gauss-Radau-Jacobi points and weights (z=-1)
def zwgrjm(Np, a, b):
    z= np.zeros(Np)
    w = np.ones(Np)*2.0

    if Np>=1:
        z[0] = -1.0;

        z[1:] = jacobz(Np-1, a, b+1.0); 
        w = jacobi(Np, a, b, z)[-1];

        fac  = pow(2.0, a+b)*gamma(a+Np)*gamma(b+Np);
        fac /= (b+Np)*gamma(Np)*gamma(a+b+Np+1.0);

        w = fac*(1-z)/(w*w)
        w[0] = w[0]*(b  + 1.0)

    return z, w

# Compute Gauss-Radau-Jacobi points and weights (z=+1)
def zwgrjp(Np, a, b):
    z= np.zeros(Np)
    w = np.ones(Np)*2.0

    if Np>=1:
        z[Np-1] = 1.0

        z[:-1] = jacobz(Np-1, a+1, b); 
        w = jacobi(Np, a, b, z)[-1];

        fac  = pow(2.0, a+b)*gamma(a+Np)*gamma(b+Np);
        fac /= (a+Np)*gamma(Np)*gamma(a+b+Np+1.0);

        w = fac*(1+z)/(w*w)
        w[Np-1] = w[Np-1]*(a  + 1.0)

    return z, w


def bench(mr, Nv, Nrho, M, quad, t0, blksize):

    # check if power of 2
    assert ((blksize & (blksize - 1)) == 0) and blksize!=0, "Need power of 2"

    # without loss of generality
    m1 = 1.0
    m2 = mr*m1

    # cases
    cases=list(it.product([0,1], repeat=2))
    # list of masses
    masses=[m1,m2]

    # the viscosity exponent: gamma = 2*(1-omega)
    gamma00 = gamma01 = gamma10 = gamma11 = 0.
    gamma = {'00': gamma00, '01': gamma01, '10': gamma10, '11': gamma11}

    # the angular scattering parameter: eta = alpha-1
    eta00 = eta01 = eta10 = eta11 = 0.
    eta = {'00': eta00, '01': eta01, '10': eta10, '11': eta11}

    # number of points in the single direction of velocity domain
    #Nv = 4
    vsize = Nv**3

    # size of compact support
    S = 3.0

    # Radius of the sphere
    R = 2*S

    # Length of the velocity domain
    #L = (3.0+2**0.5)*S/2.0
    L1 = (max(4.*mr/(1.+mr),2.0)+np.sqrt(1.0+mr)+1.0)/2.0*S;
    L2 = (max(4.*1./(1.+mr),2.)+np.sqrt(1.0+1.0/mr)+1.0)/2.0*S;
    L = max(L1,L2);
    L = 12;
    print("velocityMesh: (%s %s)"%(-L,L))

    # Number of quadrature points
    #Nrho = int(Nv/2)

    # Number of points on full sphere 
    validSph = {
          2 : "womersley-ss001-m00002.txt",
          6 : "womersley-ss003-m00006.txt",
         12 : "womersley-ss005-m00012.txt",
         32 : "womersley-ss007-m00032.txt",
         48 : "womersley-ss009-m00048.txt",
         70 : "womersley-ss011-m00070.txt",
         94 : "womersley-ss013-m00094.txt",
        120 : "womersley-ss015-m00120.txt",
        156 : "womersley-ss017-m00156.txt",
        192 : "womersley-ss019-m00192.txt"
    }
    mpath = "./"
    #M = 6
    sz = np.loadtxt(mpath + validSph[M])
    assert sz.shape == (M, 3), "Issue with spherical design points"
    sw = 4*np.pi/M
    d_sz_x = gpuarray.to_gpu(np.ascontiguousarray(sz[:,0]))
    d_sz_y = gpuarray.to_gpu(np.ascontiguousarray(sz[:,1]))
    d_sz_z = gpuarray.to_gpu(np.ascontiguousarray(sz[:,2]))

    # number of pre-integration points on full sphere
    Mpre = 192
    szpre = np.loadtxt(mpath + validSph[Mpre])
    assert szpre.shape == (Mpre, 3), "Issue with spherical design points"
    swpre = 4*np.pi/Mpre

    # gauss quadrature points for integration
    validQuads = {'gj':zwgj, 'gll':zwglj, 'grm':zwgrjm, 'grp':zwgrjp }
    assert quad in validQuads, "Valid Quads: "+str(validQuads.keys())
    qz, qw = validQuads[quad](Nrho, 0.0, 0.0)
    # scale the quadrature from [-1, 1] to [0, R] 
    qz = (R/2.0)*(1.0+qz)
    qw = ((R-0.0)/2.0)*qw
    # load quadrature zeros onto the memory
    #d_qw = gpuarray.to_gpu(qw)

    # velocity mesh weight and points 
    cw = (2.0*L/Nv)**3
    c0 = np.linspace(-L+L/Nv, L-L/Nv, Nv)
    cv = c0[np.mgrid[0:Nv, 0:Nv, 0:Nv]]
    cv = cv.reshape((3,vsize))
    
    # Krook-Wu solution
    #t0 = 6.5  # time
    tmax = t0

    # See Krook-Wu paper, $\mu$ in Eq. 14
    muMass = 4*mr/((1+mr)*(1+mr));

    # Number-densities
    n1 = 1.0;
    n2 = 1.0;
        
    # collision prefactors
    b11 = 1.0/(4.0*np.pi);
    b22 = 1.0/(4.0*np.pi);
    b12 = 1.0/(8.0*np.pi);
    b21 = b12;
    prefac = {'00': b11, '01': b12, '10': b21, '11': b22}

    lambda11 = b11*4*np.pi*n1;
    lambda12 = b21*4*np.pi*n1;
    lambda21 = b12*4*np.pi*n2;
    lambda22 = b22*4*np.pi*n2;

    p1 = lambda22 - lambda21*muMass*(3-2*muMass);
    p2 = lambda11 - lambda12*muMass*(3-2*muMass);

    A = 1./6.*(lambda11+lambda21*muMass*(3-2.0*muMass*p2/p1));
    B = 1./3.*(lambda11*p1+lambda21*muMass*(3-2.0*muMass)*p2);

    RR = A/(A*np.exp(A*tmax)-B);
    R1 = p1*RR;
    R2 = p2*RR;
    K = (n1+n2)/((n1+n2)+2*(n1*p1+n2*p2)*RR);

    P1 = 1-3*R1;
    P2 = 1-3*R2;

    dRR = -pow(A,3)*np.exp(A*tmax)/pow(A*np.exp(A*tmax)-B, 2.0);
    dR1 = p1*dRR;
    dR2 = p2*dRR;
    dP1 = -3*p1*dRR;
    dP2 = -3*p2*dRR;
    dK = -2*(n1+n2)*(n1*p1+n2*p2)/pow((n1+n2)+2*(n1*p1+n2*p2)*RR,2)*dRR;

    # KW soln at t=tmax
    h_f1 = np.empty(vsize)   # init 1st PDF
    h_f2 = np.empty(vsize)   # init 2nd PDF
    h_df1 = np.empty(vsize)  # RHS of 1st equation
    h_df2 = np.empty(vsize)  # RHS of 2nd equation
    for l in range(vsize):
        cSqr = np.dot(cv[:,l], cv[:,l])
            
        h_f1[l] = (n1*pow(m1/(2*np.pi*K),1.5)
                    *np.exp(-m1*(cSqr)/(2*K))*(P1+m1/K*R1*(cSqr)));
        h_f2[l] = (n2*pow(m2/(2*np.pi*K),1.5)
                    *np.exp(-m2*(cSqr)/(2*K))*(P2+m2/K*R2*(cSqr)));

        h_df1[l] = (h_f1[l]*(-3/(2*K)*dK+m1*(cSqr)/(2*K*K)*dK)
                    + n1*pow(m1/(2*np.pi*K),1.5)*np.exp(-m1*(cSqr)/(2*K))
                    *(dP1+(m1/K*dR1-m1*dK/(K*K)*R1)*(cSqr)));
        h_df2[l] = (h_f2[l]*(-3/(2*K)*dK+m2*(cSqr)/(2*K*K)*dK)
                    + n2*pow(m2/(2*np.pi*K),1.5)*np.exp(-m2*(cSqr)/(2*K))
                    *(dP2+(m2/K*dR2-m2*dK/(K*K)*R2)*(cSqr)));

    # precomputation
    l0 = np.concatenate((np.arange(0,Nv/2), np.arange(-Nv/2, 0)))
    l = l0[np.mgrid[0:Nv, 0:Nv, 0:Nv]]
    l = l.reshape((3,vsize)).astype(np.int32)
    d_lx = gpuarray.to_gpu(np.ascontiguousarray(l[0,:]))
    d_ly = gpuarray.to_gpu(np.ascontiguousarray(l[1,:]))
    d_lz = gpuarray.to_gpu(np.ascontiguousarray(l[2,:]))

    """
    # the vars
    #h_bb1 = np.zeros(Nrho*vsize)
    #h_bb2_00 = np.zeros(vsize)
    #h_bb2_01 = np.zeros(vsize)
    #h_bb2_10 = np.zeros(vsize)
    #h_bb2_11 = np.zeros(vsize)

    # Final precomputation
    # TODO: transform this to matrix computations
    for p in range(Nrho):
        #for q in range(M):
        #    for r in range(vsize):
        #        h_aa[ (p*M+q)*vsize + r] = (
        #            np.pi/L*qz[p]*np.dot(l[:,r], sz[q,:])
        #        );

        # sinc function has different defn in numpy
        for r in range(vsize):
            h_bb1[p*vsize+r] = np.pi/L*qz[p]*np.sqrt(np.dot(l[:,r],l[:,r]))
            h_bb2_00[r] += qw[p]*pow(qz[p],(gamma['00']+2))*16*np.pi*np.pi*(
                      np.sinc(1./L*qz[p]*np.sqrt(np.dot(l[:,r],l[:,r]))))
            h_bb2_01[r] += qw[p]*pow(qz[p],(gamma['01']+2))*16*np.pi*np.pi*(
                      np.sinc(1./L*qz[p]*np.sqrt(np.dot(l[:,r],l[:,r]))))
            h_bb2_10[r] += qw[p]*pow(qz[p],(gamma['10']+2))*16*np.pi*np.pi*(
                      np.sinc(1./L*qz[p]*np.sqrt(np.dot(l[:,r],l[:,r]))))
            h_bb2_11[r] += qw[p]*pow(qz[p],(gamma['11']+2))*16*np.pi*np.pi*(
                      np.sinc(1./L*qz[p]*np.sqrt(np.dot(l[:,r],l[:,r]))))
    print("Finished precomputation")
    """

    # transfer to device
    d_f1 = gpuarray.to_gpu(h_f1)
    d_f2 = gpuarray.to_gpu(h_f2)
    d_df1 = gpuarray.to_gpu(h_df1)
    d_df2 = gpuarray.to_gpu(h_df2)
    h_Q = np.empty_like(h_df1)
    d_Q = gpuarray.empty_like(d_df1)

    # define scratch  spaces
    d_FTf = gpuarray.empty(vsize, dtype=np.complex128)
    d_FTg = gpuarray.empty_like(d_FTf)
    d_f1C = gpuarray.empty_like(d_FTf)
    d_f2C = gpuarray.empty_like(d_FTf)
    d_QG = gpuarray.empty_like(d_FTf)
    d_t1 = gpuarray.empty(M*Nrho*vsize, dtype=np.complex128)
    d_t2 = gpuarray.empty_like(d_t1)
    d_t3 = gpuarray.empty_like(d_t1)
    #d_t4 = gpuarray.empty_like(d_t1)
    #d_t5 = gpuarray.empty_like(d_t1)

    # define complex to complex plan
    rank = 3
    n = np.array([Nv, Nv, Nv], dtype=np.int32)

    #planD2Z = cufftPlan3d(Nv, Nv, Nv, CUFFT_D2Z)
    planZ2Z_MNrho = cufftPlanMany(rank, n.ctypes.data,
        None, 1, vsize, 
        None, 1, vsize, 
        CUFFT_Z2Z, M*Nrho)
    planZ2Z = cufftPlan3d(Nv, Nv, Nv, CUFFT_Z2Z)

    # Determine the grid/block (for pre-computation)
    block = (blksize, 1, 1)
    grid = get_grid_for_block(block, vsize)
    print(block); print(grid)

    # read the template
    makosrc = "".join(open('dgfsBatchedBi_VSS_cufft.mako', 'r+').readlines())    
    src = Template(makosrc).render(
        Nrho=Nrho, M=M, vsize=vsize, sw=sw, prefac=prefac,
        cases=cases, masses=masses,
        qw=qw, qz=qz, 
        L=L, sz=sz, 
        gamma=gamma, eta=eta,
        Mpre=Mpre, szpre=szpre, swpre=swpre
    )
    #makowrite = open('dgfsBatchedBi_VSS_cufft.cu', 'w').write(src)

    # Compile the source code and retrieve the kernel
    module = compiler.SourceModule(src)

    print("Starting precomputation, this may take some time ...")
    start, end = cuda.Event(), cuda.Event()
    cuda.Context.synchronize()
    start.record()
    start.synchronize()

    d_aa = gpuarray.empty(Nrho*M*vsize, dtype=h_f1.dtype)
    precompute_aa = module.get_function("precompute_a")
    precompute_aa.prepare('PPPP')
    precompute_aa.set_cache_config(cuda.func_cache.PREFER_L1)
    precompute_aa.prepared_call(grid, block, d_lx.ptr, d_ly.ptr, d_lz.ptr, 
        d_aa.ptr)

    d_bb1 = {}; d_bb2 = {}
    precompute_bb = {}
    for cp, cq in cases:
        cpcq = str(cp)+str(cq)
        d_bb1[cpcq] = gpuarray.empty(Nrho*M*vsize, dtype=d_FTf.dtype)
        d_bb2[cpcq] = gpuarray.zeros(vsize, dtype=d_FTf.dtype)
        precompute_bb[cpcq] = module.get_function("precompute_bc_"+cpcq)
        precompute_bb[cpcq].prepare('IIdddPPPPPPPP')
        precompute_bb[cpcq].set_cache_config(cuda.func_cache.PREFER_L1)

        for p in range(Nrho):
            fac = np.pi/L*qz[p]
            fac_b = swpre*pow(qz[p],(gamma[cpcq]+2));
            fac_c = qw[p]*sw*fac_b
            for q in range(M):
                precompute_bb[cpcq].prepared_call(grid, block, 
                    np.int32(p), np.int32(q), np.float64(fac), 
                    np.float64(fac_b), np.float64(fac_c),
                    d_lx.ptr, d_ly.ptr, d_lz.ptr, 
                    d_sz_x.ptr, d_sz_y.ptr, d_sz_z.ptr, 
                    d_bb1[cpcq].ptr, d_bb2[cpcq].ptr
                )

    end.record()
    end.synchronize()
    secs = start.time_till(end)*1e-3
    print("Finished precomputation in: %fs" % (secs))

    # transform scalar to complex
    r2zKern = module.get_function("r2z")
    r2zKern.prepare('PP')
    r2zKern.set_cache_config(cuda.func_cache.PREFER_L1)

    # Prepare the cosSinMul kernel for execution
    cosSinMultKern = {}
    #computeQGKern = {}
    outKern = {}
    for cp, cq in cases:
        idx = str(cp) + str(cq)
        cosSinMultKern[idx] = module.get_function("cosSinMul_"+idx)
        cosSinMultKern[idx].prepare('PPPPP')
        cosSinMultKern[idx].set_cache_config(cuda.func_cache.PREFER_L1)

        #computeQGKern[idx] = module.get_function("computeQG_"+idx)
        #computeQGKern[idx].prepare('PPP')
        #computeQGKern[idx].set_cache_config(cuda.func_cache.PREFER_L1)

        outKern[idx] = module.get_function("output_"+idx)
        outKern[idx].prepare('PPPP')
        outKern[idx].set_cache_config(cuda.func_cache.PREFER_L1)

    # Prepare the computeQGKern kernel for execution
    computeQGKern = module.get_function("computeQG")
    computeQGKern.prepare('PPP')
    computeQGKern.set_cache_config(cuda.func_cache.PREFER_L1)

    # Prepare the prodKern kernel for execution
    prodKern = module.get_function("prod")
    prodKern.prepare('PPP')
    prodKern.set_cache_config(cuda.func_cache.PREFER_L1)

    # Prepare the ax kernel for execution
    ax2Kern = module.get_function("ax2")
    ax2Kern.prepare('PPP')
    ax2Kern.set_cache_config(cuda.func_cache.PREFER_L1)

    # Prepare the out kernel for execution
    selfSubKern = module.get_function("selfSub")
    selfSubKern.prepare('PP')
    selfSubKern.set_cache_config(cuda.func_cache.PREFER_L1)

    warmup = 10
    for i in range(warmup):
        r2zKern.prepared_call(grid, block, d_bb2['00'].ptr, d_f1C.ptr)

    # timers
    start, end = cuda.Event(), cuda.Event()
    cuda.Context.synchronize()
    start.record()
    start.synchronize()


    # single call
    def collisionSolve(idx, d_f1, d_f2, d_Q):
        # construct d_fC from d_f0
        r2zKern.prepared_call(grid, block, d_f1.ptr, d_f1C.ptr)
        r2zKern.prepared_call(grid, block, d_f2.ptr, d_f2C.ptr)

        # compute forward FFT of f | FTf = fft(f1)
        cufftExecZ2Z(planZ2Z, d_f1C.ptr, d_FTf.ptr, CUFFT_FORWARD)

        # compute forward FFT of g | FTg = fft(f2)
        cufftExecZ2Z(planZ2Z, d_f2C.ptr, d_FTg.ptr, CUFFT_FORWARD)
        
        # compute t1_{pqr} = exp(1i*m_q/(m_p+m_q)*a_{pqr})*FTf_r; 
        # t2_{pqr} = exp(-1i*m_p/(m_p+m_q)*a_{pqr})*FTg_r
        # includes scaling of FTf and FTg
        cosSinMultKern[idx].prepared_call(grid, block, d_aa.ptr, d_FTf.ptr, 
            d_FTg.ptr, d_t1.ptr, d_t2.ptr)

        # compute inverse fft 
        cufftExecZ2Z(planZ2Z_MNrho, d_t1.ptr, d_t3.ptr, CUFFT_INVERSE)
        cufftExecZ2Z(planZ2Z_MNrho, d_t2.ptr, d_t1.ptr, CUFFT_INVERSE)

        # compute t2 = t3*t1
        prodKern.prepared_call(grid, block, d_t3.ptr, d_t1.ptr, d_t2.ptr)

        # compute t1 = fft(t2)
        cufftExecZ2Z(planZ2Z_MNrho, d_t2.ptr, d_t1.ptr, CUFFT_FORWARD)

        # compute f1C_r = wrho_p*ws*b1_p*t1_r
        # scaling for forward FFT is included here
        #computeQGKern[idx].prepared_call(grid, block, d_bb1[idx].ptr, d_t1.ptr,
        #    d_f1C.ptr)
        computeQGKern.prepared_call(grid, block, d_bb1[idx].ptr, d_t1.ptr,
            d_f1C.ptr)

        # inverse fft| QG = iff(fC)  [Gain computed]
        cufftExecZ2Z(planZ2Z, d_f1C.ptr, d_QG.ptr, CUFFT_INVERSE)

        # compute FTg_r = b2_r*FTg_r
        #axKern.prepared_call(grid, block, d_bb2[idx].ptr, d_FTg.ptr)
        ax2Kern.prepared_call(grid, block, d_bb2[idx].ptr, d_FTg.ptr, 
            d_FTf.ptr)

        # inverse fft| f2C = iff(FTf)
        cufftExecZ2Z(planZ2Z, d_FTf.ptr, d_f2C.ptr, CUFFT_INVERSE)
        
        # compute Q = real(d_QG - d_f_r*d_fC)
        outKern[idx].prepared_call(grid, block, d_QG.ptr, d_f2C.ptr, d_f1.ptr, 
            d_Q.ptr)


    collisionSolve('00', d_f1, d_f1, d_Q)
    selfSubKern.prepared_call(grid, block, d_Q.ptr, d_df1.ptr)

    collisionSolve('01', d_f1, d_f2, d_Q)
    selfSubKern.prepared_call(grid, block, d_Q.ptr, d_df1.ptr)

    collisionSolve('10', d_f2, d_f1, d_Q)
    selfSubKern.prepared_call(grid, block, d_Q.ptr, d_df2.ptr)

    collisionSolve('11', d_f2, d_f2, d_Q)
    selfSubKern.prepared_call(grid, block, d_Q.ptr, d_df2.ptr)

    end.record()
    end.synchronize()
    secs = start.time_till(end)*1e-3/4.
    print("Time elapsed per kernel: %fs" % (secs))

    h_df1 = d_df1.get()
    print("Error norm (1): L_infty", np.max(np.abs(h_df1)))
    print("Error norm (1): L_1", np.sum(np.abs(h_df1)))
    print("Error norm (1): L_2", np.sqrt(np.sum(np.abs(h_df1)**2)))

    h_df2 = d_df2.get()
    print("Error norm (2): L_infty", np.max(np.abs(h_df2)))
    print("Error norm (2): L_1", np.sum(np.abs(h_df2)))
    print("Error norm (2): L_2", np.sqrt(np.sum(np.abs(h_df2)**2)))


if __name__ == "__main__":
    assert len(sys.argv)==8, "Need 7 arguments: Read the file introduction"
    try:
        maps = [float, int, int, int, str, float, int]
        mr, Nv, M, Nrho, quad, t0, blksize = (maps[i](sys.argv[i+1]) for i in range(7))
    except:
        raise ValueError("Need: mr N M Nrho quadType t0 blksize")
    
    bench(mr=mr, Nv=Nv, M=M, Nrho=Nrho, quad=quad, t0=t0, blksize=blksize)
    
