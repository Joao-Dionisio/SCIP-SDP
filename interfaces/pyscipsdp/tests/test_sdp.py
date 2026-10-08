import glob
import os

import numpy as np
import pytest
import pyscipopt

from pyscipsdp import SDPModel, SDP_CONSHDLRS


HAS_SDP_SOLVER = SDPModel().getSDPSolverName() != "none"


@pytest.fixture(params=[1, 0], ids=["sdp", "lp"])
def model(request):
    """SDPModel solving SDP relaxations (needs an SDP solver) or LPs with eigenvector cuts."""
    if request.param == 1 and not HAS_SDP_SOLVER:
        pytest.skip("PySCIPSDP was built without SDP solver")
    m = SDPModel()
    m.setParam("misc/solvesdps", request.param)
    m.hideOutput()
    return m


def test_is_pyscipopt_model():
    m = SDPModel()
    assert isinstance(m, pyscipopt.Model)
    # SCIP-SDP plugins are present
    assert m.getParam("misc/solvesdps") in (0, 1)
    m.free()


def test_sdp_solver_name():
    assert SDPModel().getSDPSolverName() in ("none", "Clarabel")


def test_without_default_plugins():
    m = SDPModel(defaultPlugins=False)
    with pytest.raises(KeyError):
        m.getParam("misc/solvesdps")
    m.includeDefaultPlugins()
    m.getParam("misc/solvesdps")


def test_addConsSdp_matrices(model):
    # [[x, 1], [1, y]] psd, y <= 4  =>  min x = 1/4
    m = model
    x = m.addVar("x", lb=0, obj=1)
    y = m.addVar("y", lb=0, ub=4)
    cons = m.addConsSdp([x, y], [[[1, 0], [0, 0]], [[0, 0], [0, 1]]], constant=[[0, -1], [-1, 0]])
    m.optimize()
    assert m.getStatus() == "optimal"
    assert m.getObjVal() == pytest.approx(0.25, abs=1e-4)
    assert cons.getConshdlrName() == "SDP"
    assert np.linalg.eigvalsh(m.getMatrixValSdp(cons)).min() >= -1e-5


def test_addConsPSD_expressions_misdp(model):
    m = model
    x = m.addVar("x", lb=0, obj=1)
    y = m.addVar("y", vtype="I", lb=1, ub=3)
    m.addConsPSD([[x, 1], [1, y]])
    m.optimize()
    assert m.getObjVal() == pytest.approx(1 / 3, abs=1e-4)
    assert m.getVal(y) == pytest.approx(3)


def test_addConsPSD_matrix_variable(model):
    # X psd, X_00 = 1, X_11 = 1, minimize X_01  =>  -1
    m = model
    X = m.addMatrixVar((2, 2), lb=-10, ub=10, name="X")
    m.addCons(X[0, 0] == 1)
    m.addCons(X[1, 1] == 1)
    m.addCons(X[0, 1] == X[1, 0])
    m.setObjective(X[0, 1] + X[1, 0])
    # SDP matrices must be symmetric, so build the matrix from the lower triangle
    M = [[X[0, 0], X[1, 0]], [X[1, 0], X[1, 1]]]
    m.addConsPSD(M)
    m.optimize()
    assert m.getObjVal() == pytest.approx(-2, abs=1e-4)


def test_getters():
    m = SDPModel()
    x = m.addVar("x")
    y = m.addVar("y")
    A = np.array([[2.0, 1.0], [1.0, 0.0]])
    B = np.array([[0.0, 0.0], [0.0, 3.0]])
    C = np.array([[1.0, 0.5], [0.5, 1.0]])
    cons = m.addConsSdp([x, y], [A, B], C, name="mysdp")
    assert cons.name == "mysdp"
    assert m.getBlocksizeSdp(cons) == 2
    assert m.getNNonzSdp(cons) == (3, 3)
    assert m.getNVarsSdp(cons) == 2
    vars, mats, const = m.getVarsSdp(cons), m.getMatricesSdp(cons), m.getConstantSdp(cons)
    assert [v.name for v in vars] == ["x", "y"]
    assert np.allclose(mats[0], A) and np.allclose(mats[1], B) and np.allclose(const, C)
    assert vars[0] is x
    assert not m.isRankOneSdp(cons)


def test_rank1():
    m = SDPModel()
    x = m.addVar("x")
    cons = m.addConsSdp([x], [np.eye(2)], rank1=True)
    assert cons.getConshdlrName() == "SDPrank1"
    assert m.isRankOneSdp(cons)


def test_input_validation():
    m = SDPModel()
    x = m.addVar("x")
    with pytest.raises(ValueError):
        m.addConsSdp([x], [[[0, 1], [2, 0]]])  # not symmetric
    with pytest.raises(ValueError):
        m.addConsSdp([x], [])
    with pytest.raises(ValueError):
        m.addConsPSD([[x * x, 0], [0, 1]])  # nonlinear
    lin = m.addCons(x >= 0)
    with pytest.raises(ValueError):
        m.getMatricesSdp(lin)


SDPDIR = os.environ.get("SCIPSDP_SOURCE", os.path.join(os.path.dirname(__file__), "..", "..", ".."))
# SCIP-SDP's CBF reader does not support second-order cones
CBF = sorted(p for p in glob.glob(os.path.join(SDPDIR, "instances", "*.cbf")) if "soc" not in os.path.basename(p))


@pytest.mark.skipif(not CBF, reason="no SCIP-SDP instances found")
@pytest.mark.parametrize("path", CBF, ids=os.path.basename)
def test_read_and_solve_cbf(path):
    m = SDPModel()
    m.hideOutput()
    m.setParam("limits/time", 60)
    if not HAS_SDP_SOLVER:
        m.setParam("misc/solvesdps", 0)
    m.readProblem(path)
    assert any(c.getConshdlrName() in SDP_CONSHDLRS for c in m.getConss())
    m.optimize()
    assert m.getStatus() in ("optimal", "infeasible", "unbounded", "inforunbd")


def test_empty_matrix_and_default_name():
    m = SDPModel()
    x = m.addVar("x")
    y = m.addVar("y")
    cons = m.addConsSdp([x, y], [np.eye(2), np.zeros((2, 2))])
    assert cons.name == "c1"
    assert m.getNNonzSdp(cons) == (2, 0)
    assert np.allclose(m.getMatricesSdp(cons)[1], 0)


def test_constant_only():
    m = SDPModel()
    cons = m.addConsSdp([], [], constant=-np.eye(3))
    assert m.getNVarsSdp(cons) == 0
    assert np.allclose(m.getConstantSdp(cons), -np.eye(3))
