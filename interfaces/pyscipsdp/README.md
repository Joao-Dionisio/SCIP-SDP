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

This is for modifying SCIP-SDP. Build SCIP-SDP with its CMake, with the SDP solver of your choice
(see SCIP-SDP's `INSTALL`), then build PySCIPSDP against it. PySCIPSDP and SCIP-SDP use the SCIP that
comes with the PySCIPOpt from PyPI, so that a single SCIP is used.

```bash
git clone https://github.com/scipopt/SCIP-SDP
cd SCIP-SDP
cmake -B build -DSCIP_DIR=/path/to/scip -DSDPS=...    # msk, sdpa, dsdp, cbl or none
cmake --build build

cd interfaces/pyscipsdp
export SCIPOPTDIR=/path/to/scip    # SCIP installation for its headers; same major.minor
                                   # version as the SCIP of PySCIPOpt
pip install .
```

By default PySCIPSDP uses SCIP-SDP's `build/` directory; set `SCIPSDPDIR` to use a SCIP-SDP
installation (`cmake --install`) instead. After changing SCIP-SDP, rebuild it and run
`pip install .` again.

PySCIPSDP is compiled against the `.pxd` files of PySCIPOpt. After upgrading PySCIPOpt, rebuild
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
