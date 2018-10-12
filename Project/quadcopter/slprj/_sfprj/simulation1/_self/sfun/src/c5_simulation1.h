#ifndef __c5_simulation1_h__
#define __c5_simulation1_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc5_simulation1InstanceStruct
#define typedef_SFc5_simulation1InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c5_sfEvent;
  boolean_T c5_isStable;
  boolean_T c5_doneDoubleBufferReInit;
  uint8_T c5_is_active_c5_simulation1;
  real_T (*c5_u)[4];
  real_T (*c5_udist)[4];
  real_T *c5_dist;
} SFc5_simulation1InstanceStruct;

#endif                                 /*typedef_SFc5_simulation1InstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c5_simulation1_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c5_simulation1_get_check_sum(mxArray *plhs[]);
extern void c5_simulation1_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
