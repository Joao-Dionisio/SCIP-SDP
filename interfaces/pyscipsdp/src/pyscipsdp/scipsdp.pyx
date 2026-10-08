# cython: language_level=3
from libc.stdlib cimport malloc, free

import numpy as np

from pyscipopt.scip cimport Constraint, Variable
from pyscipopt.scip import PY_SCIP_CALL, str_conversion, Expr

SDP_CONSHDLRS = ("SDP", "SDPrank1")


cdef SCIP_CONS* _getSdpCons(Constraint cons) except NULL:
    if cons.getConshdlrName() not in SDP_CONSHDLRS:
        raise ValueError("constraint %s is not an SDP constraint (handler %s)" % (cons.name, cons.getConshdlrName()))
    return cons.scip_cons

def _denseSymmetric(mat, n=None, what="matrix"):
    if hasattr(mat, "toarray"):
        mat = mat.toarray()
    arr = np.asarray(mat, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        raise ValueError("%s must be square, got shape %s" % (what, arr.shape))
    if n is not None and arr.shape[0] != n:
        raise ValueError("%s has size %d, expected %d" % (what, arr.shape[0], n))
    if not np.allclose(arr, arr.T):
        raise ValueError("%s must be symmetric" % what)
    return arr

def _lowerTriangularNonzeros(matrix):
    # nonzeros of the lower triangular part as (number, rows, cols, vals); the arrays get one
    # spare slot so that the address of their first element is valid even without nonzeros
    rows, cols = np.nonzero(np.tril(matrix))
    n = len(rows)
    row = np.zeros(n + 1, dtype=np.intc)
    col = np.zeros(n + 1, dtype=np.intc)
    val = np.zeros(n + 1, dtype=np.float64)
    row[:n] = rows
    col[:n] = cols
    val[:n] = matrix[rows, cols]
    return n, row, col, val

def _linearTerms(entry):
    if isinstance(entry, Expr):
        terms = []
        for term, coef in entry.terms.items():
            if len(term) == 0:
                terms.append((None, coef))
            elif len(term) == 1:
                terms.append((term[0], coef))
            else:
                raise ValueError("SDP matrix entries must be linear, got %s" % entry)
        return terms
    try:
        return [(None, float(entry))]
    except TypeError:
        raise TypeError("unsupported SDP matrix entry %r" % (entry,)) from None


cdef class SDPModel(Model):
    """
    A PySCIPOpt Model with the SCIP-SDP plugins (SDP constraint handler, SDP relaxator,
    CBF/SDPA readers, SDP heuristics, ...) included. The constructor is the one of Model.

    """

    def includeDefaultPlugins(self):
        """Includes the default SCIP plugins, the SCIP-SDP plugins and SCIP-SDP's default parameter settings."""
        PY_SCIP_CALL(SCIPSDPincludeDefaultPlugins(self._scip))

    def getSDPSolverName(self):
        """
        Gets the name of the SDP solver SCIP-SDP was built with ("none" if there is none).

        Returns
        -------
        str

        """
        return bytes(SCIPsdpiGetSolverName()).decode('utf-8')

    def addConsSdp(self, vars, matrices, constant=None, name="", rank1=False, removeduplicates=True):
        """
        Add an SDP constraint: sum_j matrices[j] * vars[j] - constant must be positive semidefinite.

        Parameters
        ----------
        vars : list of Variable
            variables of the constraint
        matrices : list of matrices
            one square symmetric matrix per variable (nested lists, numpy arrays or scipy.sparse)
        constant : matrix or None, optional
            square symmetric matrix that is subtracted; None means zero (Default value = None)
        name : str, optional
            name of the constraint (Default value = "")
        rank1 : bool, optional
            should the matrix additionally have rank one? (Default value = False)
        removeduplicates : bool, optional
            should duplicate matrix entries be removed? (Default value = True)

        Returns
        -------
        Constraint
            The newly created SDP constraint

        """
        cdef SCIP_CONS* scip_cons
        cdef int nvars
        cdef int blocksize
        cdef int i
        cdef int[::1] rowview
        cdef int[::1] colview
        cdef SCIP_Real[::1] valview
        cdef int** row = NULL
        cdef int** col = NULL
        cdef SCIP_Real** val = NULL
        cdef SCIP_VAR** _vars = NULL
        cdef int nnonz
        cdef int constnnonz
        cdef int[::1] nvarnonz
        cdef int[::1] constrow
        cdef int[::1] constcol
        cdef SCIP_Real[::1] constval

        vars = list(vars)
        matrices = list(matrices)
        if len(vars) != len(matrices):
            raise ValueError("got %d variables but %d matrices" % (len(vars), len(matrices)))
        if not vars and constant is None:
            raise ValueError("an SDP constraint needs at least one variable or a constant matrix")

        matrices = [_denseSymmetric(m, what="matrices[%d]" % i) for i, m in enumerate(matrices)]
        blocksize = matrices[0].shape[0] if matrices else _denseSymmetric(constant, what="constant").shape[0]
        for i, m in enumerate(matrices):
            _denseSymmetric(m, blocksize, "matrices[%d]" % i)
        if constant is None:
            constant = np.zeros((blocksize, blocksize))
        constant = _denseSymmetric(constant, blocksize, "constant")

        # SCIP-SDP takes the nonzeros of the lower triangular parts as (row, col, val) arrays;
        # the numpy arrays own the memory, the C arrays below only point into them
        nvars = len(vars)
        nonzeros = [_lowerTriangularNonzeros(m) for m in matrices]
        nvarnonz = np.array([n for n, _, _, _ in nonzeros] + [0], dtype=np.intc)
        nnonz = sum(nvarnonz)
        constnnonz, constrow, constcol, constval = _lowerTriangularNonzeros(constant)

        if name == '':
            name = 'c'+str(SCIPgetNConss(self._scip)+1)

        try:
            row = <int**>malloc(max(nvars, 1) * sizeof(int*))
            col = <int**>malloc(max(nvars, 1) * sizeof(int*))
            val = <SCIP_Real**>malloc(max(nvars, 1) * sizeof(SCIP_Real*))
            _vars = <SCIP_VAR**>malloc(max(nvars, 1) * sizeof(SCIP_VAR*))
            if row == NULL or col == NULL or val == NULL or _vars == NULL:
                raise MemoryError()

            for i in range(nvars):
                _, rowview, colview, valview = nonzeros[i]
                row[i] = &rowview[0]
                col[i] = &colview[0]
                val[i] = &valview[0]
                _vars[i] = (<Variable?>vars[i]).scip_var

            if rank1:
                PY_SCIP_CALL(SCIPcreateConsSdpRank1(self._scip, &scip_cons, str_conversion(name), nvars,
                    nnonz, blocksize, &nvarnonz[0], col, row, val, _vars, constnnonz,
                    &constcol[0], &constrow[0], &constval[0], removeduplicates))
            else:
                PY_SCIP_CALL(SCIPcreateConsSdp(self._scip, &scip_cons, str_conversion(name), nvars,
                    nnonz, blocksize, &nvarnonz[0], col, row, val, _vars, constnnonz,
                    &constcol[0], &constrow[0], &constval[0], removeduplicates))
        finally:
            free(row)
            free(col)
            free(val)
            free(_vars)

        PY_SCIP_CALL(SCIPaddCons(self._scip, scip_cons))
        pyCons = self._getOrCreateCons(scip_cons)
        PY_SCIP_CALL(SCIPreleaseCons(self._scip, &scip_cons))

        return pyCons

    def addConsPSD(self, matrix, name="", rank1=False):
        """
        Add an SDP constraint: matrix must be positive semidefinite.

        Parameters
        ----------
        matrix : square matrix of linear expressions
            symmetric matrix (nested lists or numpy array, e.g. from addMatrixVar) whose
            entries are linear expressions, variables or numbers
        name : str, optional
            name of the constraint (Default value = "")
        rank1 : bool, optional
            should the matrix additionally have rank one? (Default value = False)

        Returns
        -------
        Constraint
            The newly created SDP constraint

        """
        entries = np.asarray(matrix, dtype=object)
        if entries.ndim != 2 or entries.shape[0] != entries.shape[1]:
            raise ValueError("matrix must be square, got shape %s" % (entries.shape,))
        n = entries.shape[0]

        coefs = {}
        vars = []
        constant = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                for var, coef in _linearTerms(entries[i, j]):
                    if var is None:
                        constant[i, j] += coef
                        continue
                    if var.ptr() not in coefs:
                        coefs[var.ptr()] = np.zeros((n, n))
                        vars.append(var)
                    coefs[var.ptr()][i, j] += coef

        matrices = [coefs[var.ptr()] for var in vars]
        for m in matrices + [constant]:
            if not np.allclose(m, m.T):
                raise ValueError("matrix must be symmetric")

        # matrix = sum_j A_j x_j + C  is psd  <=>  sum_j A_j x_j - (-C)  is psd
        return self.addConsSdp(vars, matrices, -constant, name=name, rank1=rank1)

    def getNVarsSdp(self, Constraint sdp_cons):
        """
        Gets number of variables in SDP constraint.

        Parameters
        ----------
        sdp_cons : Constraint
            SDP constraint to get the number of variables from.

        Returns
        -------
        int

        """
        return SCIPconsSdpGetNVars(self._scip, _getSdpCons(sdp_cons))

    def getVarsSdp(self, Constraint sdp_cons):
        """
        Gets variables in SDP constraint.

        Parameters
        ----------
        sdp_cons : Constraint
            SDP constraint to get the variables from.

        Returns
        -------
        list of Variable

        """
        cdef SCIP_CONS* scip_cons = _getSdpCons(sdp_cons)
        cdef int nvars = SCIPconsSdpGetNVars(self._scip, scip_cons)
        cdef SCIP_VAR** _vars = SCIPconsSdpGetVars(self._scip, scip_cons)
        cdef int i

        return [self._getOrCreateVar(_vars[i]) for i in range(nvars)]

    def getBlocksizeSdp(self, Constraint sdp_cons):
        """
        Gets the size of the matrix of an SDP constraint.

        Parameters
        ----------
        sdp_cons : Constraint
            SDP constraint to get the blocksize from.

        Returns
        -------
        int

        """
        return SCIPconsSdpGetBlocksize(self._scip, _getSdpCons(sdp_cons))

    def getNNonzSdp(self, Constraint sdp_cons):
        """
        Gets number of nonzeros (lower triangular part) in SDP constraint.

        Parameters
        ----------
        sdp_cons : Constraint
            SDP constraint to get the number of nonzeros from.

        Returns
        -------
        tuple of int
            number of nonzeros of the variable matrices and of the constant matrix

        """
        cdef int nnonz
        cdef int constnnonz

        PY_SCIP_CALL(SCIPconsSdpGetNNonz(self._scip, _getSdpCons(sdp_cons), &nnonz, &constnnonz))

        return nnonz, constnnonz

    def getMatricesSdp(self, Constraint sdp_cons):
        """
        Gets the matrices of the SDP constraint sum_j A_j x_j - A_0 psd, ordered as getVarsSdp().

        Parameters
        ----------
        sdp_cons : Constraint
            SDP constraint to get the matrices from.

        Returns
        -------
        list of numpy.ndarray
            the dense matrices A_j

        """
        cdef SCIP* scip = self._scip
        cdef SCIP_CONS* scip_cons = _getSdpCons(sdp_cons)
        cdef int blocksize = SCIPconsSdpGetBlocksize(scip, scip_cons)
        cdef int nvars = SCIPconsSdpGetNVars(scip, scip_cons)
        cdef SCIP_Real[:, ::1] buf
        cdef int j

        matrices = []
        for j in range(nvars):
            matrix = np.zeros((blocksize, blocksize))
            buf = matrix
            PY_SCIP_CALL(SCIPconsSdpGetFullAj(scip, scip_cons, j, &buf[0, 0]))
            matrices.append(matrix)

        return matrices

    def getConstantSdp(self, Constraint sdp_cons):
        """
        Gets the constant matrix A_0 of the SDP constraint sum_j A_j x_j - A_0 psd.

        Parameters
        ----------
        sdp_cons : Constraint
            SDP constraint to get the constant matrix from.

        Returns
        -------
        numpy.ndarray

        """
        cdef SCIP* scip = self._scip
        cdef SCIP_CONS* scip_cons = _getSdpCons(sdp_cons)
        cdef int blocksize = SCIPconsSdpGetBlocksize(scip, scip_cons)
        cdef SCIP_Real[:, ::1] buf

        constant = np.zeros((blocksize, blocksize))
        buf = constant
        PY_SCIP_CALL(SCIPconsSdpGetFullConstMatrix(scip, scip_cons, &buf[0, 0]))

        return constant

    def isRankOneSdp(self, Constraint sdp_cons):
        """
        Returns whether the SDP constraint additionally requires rank one.

        Parameters
        ----------
        sdp_cons : Constraint
            SDP constraint to check.

        Returns
        -------
        bool

        """
        return bool(SCIPconsSdpShouldBeRankOne(_getSdpCons(sdp_cons)))

    def getMatrixValSdp(self, sdp_cons, sol=None):
        """
        Gets the value of the matrix sum_j A_j x_j - A_0 of an SDP constraint in a solution.

        Parameters
        ----------
        sdp_cons : Constraint
            SDP constraint to evaluate.
        sol : Solution or None, optional
            solution to evaluate in; None means the best solution (Default value = None)

        Returns
        -------
        numpy.ndarray

        """
        if sol is None:
            sol = self.getBestSol()

        matrix = -self.getConstantSdp(sdp_cons)
        for var, A in zip(self.getVarsSdp(sdp_cons), self.getMatricesSdp(sdp_cons)):
            matrix = matrix + self.getSolVal(sol, var) * A

        return matrix
