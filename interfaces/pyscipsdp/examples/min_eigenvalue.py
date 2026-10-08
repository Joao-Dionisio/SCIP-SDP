"""Maximise the smallest eigenvalue of A(x) = A0 + x1*A1 + x2*A2 with x1 + x2 = 1, x >= 0.

    max t  s.t.  A(x) - t*I  is positive semidefinite
"""
import numpy as np

from pyscipsdp import SDPModel

A0 = np.array([[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 1.0]])
A1 = np.diag([0.0, 0.0, 2.0])
A2 = np.array([[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

m = SDPModel()
x1 = m.addVar("x1", lb=0)
x2 = m.addVar("x2", lb=0)
t = m.addVar("t", lb=None)
m.addCons(x1 + x2 == 1)

# matrix of linear expressions; entries may be numbers, variables or expressions
M = [[A0[i, j] + A1[i, j] * x1 + A2[i, j] * x2 - (t if i == j else 0) for j in range(3)] for i in range(3)]
cons = m.addConsPSD(M, name="psd")
m.setObjective(t, "maximize")

m.optimize()

print("t  =", m.getVal(t))
print("x  =", m.getVal(x1), m.getVal(x2))
print("eigenvalues of A(x):", np.linalg.eigvalsh(m.getMatrixValSdp(cons) + m.getVal(t) * np.eye(3)))
