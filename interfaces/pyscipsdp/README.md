# PySCIPSDP

Python interface for [SCIP-SDP](https://github.com/scipopt/SCIP-SDP), built on top of
[PySCIPOpt](https://github.com/scipopt/PySCIPOpt).

```python
from pyscipsdp import SDPModel

m = SDPModel()
x = m.addVar("x", obj=1)
y = m.addVar("y", vtype="I", lb=1, ub=3)
m.addConsPSD([[x, 1], [1, y]])        # [[x, 1], [1, y]] is positive semidefinite
m.optimize()
```

`SDPModel` is a subclass of `pyscipopt.Model` (`cdef class SDPModel(Model)`), so everything
PySCIPOpt can do works unchanged.

## Installation

```bash
pip install pyscipsdp
```

(Wheels are not published yet.)

## Building from source

This is for modifying SCIP-SDP. It needs a SCIP installation of the same major.minor version as the
SCIP of PySCIPOpt (for its headers and CMake files), and CMake.

```bash
git clone https://github.com/scipopt/SCIP-SDP
cd SCIP-SDP/interfaces/pyscipsdp
pip install .
```

The SCIP installation is looked for in the conda environment, `/usr/local`, `/opt/homebrew` and
`/usr`; set `SCIPOPTDIR` to use another one.

`pip install .` builds SCIP-SDP with its CMake into `SCIP-SDP/build` and then PySCIPSDP against it.
The SDP solver is chosen as in SCIP-SDP's build: `SDPS=msk|sdpa|dsdp|cbl|none` (default `none`), with
the solver found as described in SCIP-SDP's `INSTALL`. After changing SCIP-SDP, run `pip install .`
again; to use an already installed SCIP-SDP instead, set `SCIPSDPDIR`.

PySCIPSDP and SCIP-SDP use the SCIP that comes with the PySCIPOpt from PyPI, so that a single SCIP is
used. PySCIPSDP is compiled against the `.pxd` files of PySCIPOpt; after upgrading PySCIPOpt, rebuild
PySCIPSDP. Building from source is currently tested on macOS only.

## API (in addition to `pyscipopt.Model`)

| method | |
|---|---|
| `addConsSdp(vars, matrices, constant=None, name="", rank1=False)` | `Σ_j matrices[j]·vars[j] − constant` is psd. Matrices may be nested lists, numpy arrays or scipy.sparse. |
| `addConsPSD(matrix, name="", rank1=False)` | `matrix` is psd, where the entries are linear expressions, variables or numbers |
| `getNVarsSdp(sdp_cons)`, `getVarsSdp(sdp_cons)` | |
| `getBlocksizeSdp(sdp_cons)`, `getNNonzSdp(sdp_cons)` | |
| `getMatricesSdp(sdp_cons)`, `getConstantSdp(sdp_cons)` | the dense matrices `A_j` and `A_0` |
| `isRankOneSdp(sdp_cons)` | |
| `getMatrixValSdp(sdp_cons, sol=None)` | value of `Σ A_j x_j − A_0` in a solution |
| `getSDPSolverName()` | the SDP solver SCIP-SDP was built with |

The CBF and SDPA readers are included, so `m.readProblem("instance.cbf")` works as well.

## Tests

```bash
pytest tests          # also solves the instances in SCIP-SDP's instances/ directory
```
