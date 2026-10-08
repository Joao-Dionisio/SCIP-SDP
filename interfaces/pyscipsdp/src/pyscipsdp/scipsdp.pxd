from pyscipopt.scip cimport SCIP, SCIP_CONS, SCIP_VAR, SCIP_Real, SCIP_Bool, SCIP_RETCODE
from pyscipopt.scip cimport SCIPgetNConss, SCIPaddCons, SCIPreleaseCons
from pyscipopt.scip cimport Model

cdef extern from "scipsdpdefplugins.h":
    SCIP_RETCODE SCIPSDPincludeDefaultPlugins(SCIP* scip)

cdef extern from "sdpi/sdpi.h":
    const char* SCIPsdpiGetSolverName()

cdef extern from "cons_sdp.h":
    SCIP_RETCODE SCIPcreateConsSdp(SCIP* scip,
                                   SCIP_CONS** cons,
                                   const char* name,
                                   int nvars,
                                   int nnonz,
                                   int blocksize,
                                   int* nvarnonz,
                                   int** col,
                                   int** row,
                                   SCIP_Real** val,
                                   SCIP_VAR** vars,
                                   int constnnonz,
                                   int* constcol,
                                   int* constrow,
                                   SCIP_Real* constval,
                                   SCIP_Bool removeduplicates)
    SCIP_RETCODE SCIPcreateConsSdpRank1(SCIP* scip,
                                        SCIP_CONS** cons,
                                        const char* name,
                                        int nvars,
                                        int nnonz,
                                        int blocksize,
                                        int* nvarnonz,
                                        int** col,
                                        int** row,
                                        SCIP_Real** val,
                                        SCIP_VAR** vars,
                                        int constnnonz,
                                        int* constcol,
                                        int* constrow,
                                        SCIP_Real* constval,
                                        SCIP_Bool removeduplicates)
    int SCIPconsSdpGetNVars(SCIP* scip, SCIP_CONS* cons)
    SCIP_VAR** SCIPconsSdpGetVars(SCIP* scip, SCIP_CONS* cons)
    int SCIPconsSdpGetBlocksize(SCIP* scip, SCIP_CONS* cons)
    SCIP_RETCODE SCIPconsSdpGetNNonz(SCIP* scip, SCIP_CONS* cons, int* nnonz, int* constnnonz)
    SCIP_RETCODE SCIPconsSdpGetFullAj(SCIP* scip, SCIP_CONS* cons, int j, SCIP_Real* Aj)
    SCIP_RETCODE SCIPconsSdpGetFullConstMatrix(SCIP* scip, SCIP_CONS* cons, SCIP_Real* mat)
    SCIP_Bool SCIPconsSdpShouldBeRankOne(SCIP_CONS* cons)

cdef class SDPModel(Model):
    pass
