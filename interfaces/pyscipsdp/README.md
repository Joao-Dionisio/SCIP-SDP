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

This is for modifying SCIP-SDP. PySCIPSDP lives in SCIP-SDP's repository and compiles the
SCIP-SDP sources of that checkout into its extension. It links against the SCIP that comes with
the PySCIPOpt from PyPI, so that a single SCIP is used.

```bash
git clone https://github.com/scipopt/SCIP-SDP
cd SCIP-SDP/interfaces/pyscipsdp
export SCIPOPTDIR=/path/to/scip          # SCIP installation for its headers; same major.minor
                                         # version as the SCIP of PySCIPOpt
export SDPS=cbl                          # SDP solver: none (default) or cbl (Clarabel)
export CLARABEL_DIR=/path/to/Clarabel.cpp
pip install .
```

After changing SCIP-SDP, run `pip install .` again. `SCIPSDP_SOURCE` builds against a different
SCIP-SDP checkout.

With `SDPS=none`, SCIP-SDP can only solve with LP relaxations and eigenvector cuts (set
`misc/solvesdps` to 0). For `SDPS=cbl`, build Clarabel's C library first (needs Rust):

```bash
git clone --recurse-submodules https://github.com/oxfordcontrol/Clarabel.cpp
cd Clarabel.cpp/rust_wrapper
cargo build --release --features "sdp,clarabel/sdp-accelerate"   # Linux: clarabel/sdp-openblas
```

PySCIPSDP is compiled against the `.pxd` files of PySCIPOpt. After upgrading PySCIPOpt, rebuild
PySCIPSDP. Building from source is supported on Linux and macOS.

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
