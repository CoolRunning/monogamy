#!/usr/bin/env python
"""
    Monogamy of quantum correlations
    A python module that lets you explore what is possible and what is not

    @author Marcus Märtens  (mmarcusx@gmail.com)
"""
import numpy as np
import cvxpy as cvx
from itertools import product
from math import log

#def idx2bool(l):
#    """ Takes a list of indices ("""
#    return [i in l for i


def is_density_operator(A, ftol=1.0e-4):
    """ returns True if A is a valid density operator
        @param  A   numpy-array
        @param  ftol    tolerance value for numerical instability"""
    # A needs to be quadratic
    if (A.shape[0] !=  A.shape[1]):
        return False

    # The trace of A needs to be one
    if (np.fabs(A.trace() - 1.0) > ftol):
        return False

    # A needs to be hermitian
    if not(np.allclose(A, np.matrix.getH(A), atol=ftol)):
        return False

    # A needs to be positive semi-definite
    return is_psd(A, ftol)


def werner(w):
    """ returns a Werner-state defined as
        rho_AB = w|\phi_+><\phi_+| + 0.5*(1-w)*(I_A \dot I_B)
        @parameter  w   interpolates between maximally entangled (w=1) and no correlations (w=0)
        @return numpy-Array
    """
    if not(0 <= w <= 1):
        raise Exception(ValueError, "Parameter w needs to be between 0 and 1!")

    rho_AB = np.zeros((4,4))
    rho_AB[0,0] = 0.25 + w*0.25
    rho_AB[3,3] = 0.25 + w*0.25
    rho_AB[1,1] = 0.25 - w*0.25
    rho_AB[2,2] = 0.25 - w*0.25
    rho_AB[0,3] = w*0.5
    rho_AB[3,0] = w*0.5

    return rho_AB


def is_psd(A, ftol=1.0e-4):
    """ checks whether A is a posi-semidefinite matrix by
        computing all eigenvalues. """
    for (i,j),_ in np.ndenumerate(A):
        if (np.linalg.det(minor(A,i,j))) < -1.0 * ftol:
            return False
    return True


def minor(A,i,j):
    """ computes the minor of a numpy array A.
        @parameters A   numpy array
        @parameters i   row index
        @parameters j   col index

        @see    http://stackoverflow.com/questions/3858213"""
    return A[np.array(list(range(i))+list(range(i+1,A.shape[0])))[:,np.newaxis],
             np.array(list(range(j))+list(range(j+1,A.shape[1])))]


def partial_trace(A, sel, d):
    """ Partial trace operation.

        @parameter  A   numpy-array
        @parameter  d   vector describing the dimension of each subsystem
        @parameter  sel list of booleans, True marks registers that are gonna be traced out
                        can also be given as a list of indices

        @return numpy-array
    """
    # transform sel if given as a list of indices:
    try:
        if type(sel[0]) == int:
            sel = [i in sel for i in list(range(len(d)))]
    except IndexError:
        print('Warning: no subsystem to be traced out.')
        sel = [False] * len(d)

    # compute the new dimension of the reduced state
    dr = [dim for (dim, b) in zip(d,sel) if not b]  # dimensions not traced out
    dr_inv = [dim for (dim, b) in zip(d,sel) if b]  # dimensions traced out
    msize = np.prod(dr)
    R = np.zeros( [msize]*2 )

    # we produce the set of all valid indices in the reduced system and loop over it
    mbase = list(product(*[range(x) for x in dr]))
    for v in product(mbase, mbase):
        val = 0.0
        r_idx = list(v[0])
        c_idx = list(v[1])
        for z in list(product(*[list(range(x)) for x in dr_inv])):
            # we walk through all dimensions: if the dimension needs to be
            # traced out, we force the same subindex in the row and column selection
            r_idx2, c_idx2 = [], []
            i, j = 0, 0
            for x in sel:
                if x:
                    r_idx2.append(z[j])
                    c_idx2.append(z[j])
                    j += 1
                else:
                    r_idx2.append(r_idx[i])
                    c_idx2.append(c_idx[i])
                    i += 1
            val += A[get_coeff_idx( (r_idx2, c_idx2), d )]
        R[get_coeff_idx(v,dr)] = val

    return R


def partial_trace_constraints(rho_e, rho_c, sel, d):
    """ This function generates a set of constraints that reflect the
        conditions given by a partial-trace relation. For example, the
        condition:
        tr_A(rho_AB) = rho_B
        will be translated into 4 equality constraints, each covering
        up one entry of the matrix rho_B

        @parameter  rho_e cvxpy expression for the matrix that we optimize
        @parameter  sel bitmask that marks register to be traced out by True
                        can also be given as a list of indices
        @parameter  rho_c resulting density matrix as a numpy-array
        @parameter  d   vector describing the dimension of each subsystem

        @return list of cvxpy expression
    """
    # transform sel if given as a list of indices:
    try:
        if type(sel[0]) == int:
            sel = [i in sel for i in list(range(len(d)))]
    except IndexError:
        print('Warning: no subsystem to be traced out.')
        sel = [False] * len(d)

    constr = []

    # compute the new dimension of the reduced state
    dr = [dim for (dim, b) in zip(d,sel) if not b]  # dimensions not traced out
    dr_inv = [dim for (dim, b) in zip(d,sel) if b] # dimensions traced out

    # we produce the set of all valid indices in the reduced system and loop over it
    mbase = list(product(*[range(x) for x in dr]))
    for v in product(mbase, mbase):
        r_idx = list(v[0])
        c_idx = list(v[1])
        expr = []
        for z in list(product(*[list(range(x)) for x in dr_inv])):
            # we walk through the dimension: if the dimension needs to be
            # traced out, we force the same subindex in the row and column selection
            r_idx2, c_idx2 = [], []
            i, j = 0, 0
            for x in sel:
                if x:
                    r_idx2.append(z[j])
                    c_idx2.append(z[j])
                    j += 1
                else:
                    r_idx2.append(r_idx[i])
                    c_idx2.append(c_idx[i])
                    i += 1
            if expr == []:
                expr = [rho_e[get_coeff_idx( (r_idx2, c_idx2), d)]]
            else:
                expr = [expr[0] + rho_e[get_coeff_idx( (r_idx2, c_idx2), d)]]
        constr += [ expr[0] == rho_c[get_coeff_idx(v,dr)] ]

    return constr


def get_coeff_idx(v, d):
    """ Used to translate from a multi-index into a flat index that
        can be used to address a numpy-matrix. For example, if
        system A and system B have a binary basis, the joined system
        AB has 4 basis vectors: (0,0),(0,1),(1,0),(1,1). These
        vectors are mapped to integer:
        (0,0) ->  0
        (0,1) ->  1
        (1,0) ->  2
        (1,1) ->  3
    """
    row = sum([v[0][i] * np.prod(d[i+1:]) for i in range(len(d))])
    col = sum([v[1][i] * np.prod(d[i+1:]) for i in range(len(d))])
    return (int(row), int(col))


def construct_global_state(n, d, margs, ftol=10e-4, use_scs=False):
    """ This function takes a set of marginal states and tries to construct a global state in which
        all parties share the respective correlations.
        @parameter  n   number of registers of the global state
        @parameter  d   list which describes the dimension of each register
        @parameter margs    list of tuples (rho, m) where rho is a numpy-array describing a marginal state
                                and m is bitmask (list of booleans) that mark which register is active in the marginal
        @parameter  use_scs if True, the SCS solver will be used instead of CVXOPT. SCS is faster, but less accurate.
        @return (state, rhoG)   where state is 'optimal' if a global state was found and 'infeasible' if not. If

                                a state was found, it is returned by rhoG as a numpy-array.
    """
    # transform sel if given as a list of indices:
    marginals = []
    for (rho, m) in margs:
        try:
            if type(m[0]) == int:
                if max(m) >= len(d):
                    raise Exception(ValueError, "Error in Indexing of subsystem (indexing starts with 0)!")
                m = [i in m for i in list(range(len(d)))]
        except IndexError:
            raise Exception(IndexError, "No system to be traced out!")
            m = [False] * len(d)
        marginals.append( (rho, m) )


    # some sanity checks
    if len(d) != n:
        raise Exception(ValueError, "The number of registers does not equal the length of the given vector of dimensions!")
    for (rho, m) in marginals:
        if not is_density_operator(rho, ftol):
            raise Exception(ValueError, "One of the given marginals is not a valid density operator!")

    for (rho, m) in [marginal for marginal in marginals]:
        # check whether the marginal is consistent with its bitmask
        if (np.prod([x for x,b in zip(d,m) if not b])) != rho.shape[0]:
            raise Exception(ValueError, "The dimensions of the marginal with the traced out parameters do not match")

    # create expression for global state
    rho_global = cvx.Semidef(int(np.prod(d)))

    # adding constraints
    constr = [cvx.trace(rho_global) == 1]
    for marginal in marginals:
        constr += partial_trace_constraints(rho_global, marginal[0], marginal[1], d)

    # there is no objective function, so we define a dummy
    obj = cvx.Minimize(1)

    # define problem
    prob = cvx.Problem(obj, constr)

    # apply solver
    if use_scs:
        prob.solve(solver=cvx.SCS)
    else:
        prob.solve(solver=cvx.CVXOPT)

    # extract solution (if found) and return
    return (prob.status, np.array(rho_global.value) if prob.status == cvx.OPTIMAL else None)


def entropy(A):
    """ Computes the Von-Neumann-Entropy of a density-operator A
        @parameter  A   numpy-array
        @return Von-Neumann-Entropy
    """
    return -1 * sum([e * log(e,2) for e in np.linalg.eigvalsh(A) if e != 0.0])
