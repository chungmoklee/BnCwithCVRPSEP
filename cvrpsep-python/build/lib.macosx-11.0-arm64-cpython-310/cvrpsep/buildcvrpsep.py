import os
import glob
from cffi import FFI


ffibuilder = FFI()

FILE_PATH = os.path.dirname(__file__)

SRC_PATH = os.path.abspath(os.path.join(FILE_PATH, 'src'))
INC_PATH = os.path.abspath(os.path.join(FILE_PATH, 'src'))

ffibuilder.cdef(
"""
#define CMGR_CT_MIN_ROUTES        101
#define CMGR_CT_NODE_DEGREE       102  /* Degree = 2 for each customer. */
#define CMGR_CT_CAP               103  /* Capacity constraint. */
#define CMGR_CT_GENCAP            104  /* Generalized capacity constraint. */
#define CMGR_CT_FCI               104  /* For public version. */
#define CMGR_CT_TWOMATCHING       105
#define CMGR_CT_COMB              106
#define CMGR_CT_STR_COMB          107  /* Strengthened comb. */
#define CMGR_CT_HYPOTOUR          108
#define CMGR_CT_EXT_HYPOTOUR      109
#define CMGR_CT_MSTAR             110  /* Homogeneous multistar. */
#define CMGR_CT_WMSTAR            111  /* Weak multistar. */
#define CMGR_CT_DJCUT             112  /* Disjunctive cut. */
#define CMGR_CT_GOMORY            113  /* By variable numbers */
#define CMGR_CT_TWOEDGES_HYPOTOUR 114  /* 2EH inequality */
#define CMGR_BT_CLIQUE_DOWN       201  /* x(S:S) <= RHS */
#define CMGR_BT_CLIQUE_UP         202  /* x(S:S) >= RHS */
#define CMGR_BT_STAR_DOWN         301  /* x(i:F) <=RHS */
#define CMGR_BT_STAR_UP           302  /* x(i:F) >=RHS */

#define CMGR_CT_SLB               401  /* x(F) >= RHS. Simple lower bound */


typedef struct
{
    int CType; /* Constraint Type. */
    int Key;
    int IntListSize;
    int *IntList;
    int ExtListSize;
    int *ExtList;
    int CListSize;
    int *CList;
    double *CoeffList;
    int A,B,L; /* For MSTARs: Lambda=L/B, Sigma=A/B. */
    double RHS;
    int BranchLevel;
    int GlobalNr;
} CnstrRecord;
typedef CnstrRecord *CnstrPointer;


typedef CnstrPointer *CnstrPointerList;

typedef struct
{
    CnstrPointerList CPL;
    int Dim;
    int Size;
} CnstrMgrRecord;
typedef CnstrMgrRecord *CnstrMgrPointer;


void CAPSEP_SeparateCapCuts(
    int NoOfCustomers,
    int *Demand,
    int CAP,
    int NoOfEdges,
    int *EdgeTail,
    int *EdgeHead,
    double *EdgeX,
    CnstrMgrPointer CMPExistingCuts,
    int MaxNoOfCuts,
    double EpsForIntegrality,
    char *IntegerAndFeasible,
    double *MaxViolation,
    CnstrMgrPointer CutsCMP
);

void MSTARSEP_SeparateMultiStarCuts(
    int NoOfCustomers,
    int *Demand,
    int CAP,
    int NoOfEdges,
    int *EdgeTail,
    int *EdgeHead,
    double *EdgeX,
    CnstrMgrPointer CMPExistingCuts,
    int MaxNoOfCuts,
    double *MaxViolation,
    CnstrMgrPointer CutsCMP
);

void GLMSEP_SeparateGLM(
    int NoOfCustomers,
    int *Demand,
    int CAP,
    int NoOfEdges,
    int *EdgeTail,
    int *EdgeHead,
    double *EdgeX,
    int *CustList,
    int *CustListSize,
    double *Violation
);

void FCISEP_SeparateFCIs(
    int NoOfCustomers,
    int *Demand,
    int CAP,
    int NoOfEdges,
    int *EdgeTail,
    int *EdgeHead,
    double *EdgeX,
    CnstrMgrPointer CMPExistingCuts,
    int MaxNoOfTreeNodes,
    int MaxNoOfCuts,
    double *MaxViolation,
    CnstrMgrPointer CutsCMP
);

void COMBSEP_SeparateCombs(
    int NoOfCustomers,
    int *Demand,
    int CAP,
    int QMin,
    int NoOfEdges,
    int *EdgeTail,
    int *EdgeHead,
    double *EdgeX,
    int MaxNoOfCuts,
    double *MaxViolation,
    CnstrMgrPointer CutsCMP
);

void HTOURSEP_SeparateHTours(
    int NoOfCustomers,
    int *Demand,
    int CAP,
    int NoOfEdges,
    int *EdgeTail,
    int *EdgeHead,
    double *EdgeX,
    CnstrMgrPointer CMPExistingCuts,
    int MaxNoOfCuts,
    double *MaxViolation,
    CnstrMgrPointer CutsCMP
);

void BRNCHING_GetCandidateSets(
    int NoOfCustomers,
    int *Demand,
    int CAP,
    int NoOfEdges,
    int *EdgeTail,
    int *EdgeHead,
    double *EdgeX,
    CnstrMgrPointer CMPExistingCuts,
    double BoundaryTarget,
    int MaxNoOfSets,
    CnstrMgrPointer SetsCMP
);

void CMGR_CreateCMgr(CnstrMgrPointer *CMP, int Dim);

void CMGR_MoveCnstr(CnstrMgrPointer SourcePtr,
                    CnstrMgrPointer SinkPtr,
                    int SourceIndex,
                    int SinkIndex);
"""
)



src_files = [
    'basegrph.c', 
    'brnching.c', 
    'combsep.c', 
    'cutbase.c', 
    'fcits.c', 
    'hpmstar.c', 
    'memmod.c', 
    'newhtour.c', 
    'strngcmp.c', 
    'binpack.c', 
    'capsep.c', 
    'compcuts.c', 
    'fcapfix.c', 
    'glmsep.c', 
    'htoursep.c', 
    'mstarsep.c', 
    'sort.c', 
    'twomatch.c', 
    'blocks.c', 
    'cnstrmgr.c', 
    'compress.c', 
    'fcisep.c', 
    'grsearch.c', 
    'intap.c', 
    'mxf.c', 
    'strcomb.c',
]


# with open(os.path.join(PATH, 'src/capsep.c'), 'r') as f:
ffibuilder.set_source(
    "cvrpsep",
    """
    #include <stdlib.h>
    #include <stdio.h>
    #include "memmod.h"
    #include "basegrph.h"
    #include "sort.h"
    #include "cnstrmgr.h"
    #include "cutbase.h"
    #include "compcuts.h"
    #include "compress.h"
    #include "fcapfix.h"
    #include "grsearch.h"
    #include "capsep.h"
    #include "mstarsep.h"
    #include "glmsep.h"
    #include "fcisep.h"
    #include "combsep.h"
    #include "htoursep.h"
    #include "brnching.h"
    """,
    libraries=["c"],
    sources=[os.path.join(SRC_PATH, f) for f in src_files],
    # sources=src_files,
    include_dirs=[INC_PATH],
    extra_compile_args=['-O3']
)

# print(f'CWD: {os.getcwd()}')
# print(f'{INC_PATH=}')
# print(f'{src_files=}')
# print(f'{os.listdir(INC_PATH)}')
    

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
