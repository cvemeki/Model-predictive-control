#ifndef __c5_simulation2_h__
#define __c5_simulation2_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc5_simulation2InstanceStruct
#define typedef_SFc5_simulation2InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c5_sfEvent;
  boolean_T c5_isStable;
  boolean_T c5_doneDoubleBufferReInit;
  uint8_T c5_is_active_c5_simulation2;
  uint32_T c5_method;
  boolean_T c5_method_not_empty;
  uint32_T c5_state;
  boolean_T c5_state_not_empty;
  uint32_T c5_b_state[2];
  boolean_T c5_b_state_not_empty;
  uint32_T c5_c_state[625];
  boolean_T c5_c_state_not_empty;
  real_T (*c5_d)[3];
  real_T *c5_dist;
} SFc5_simulation2InstanceStruct;

#endif                                 /*typedef_SFc5_simulation2InstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c5_simulation2_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c5_simulation2_get_check_sum(mxArray *plhs[]);
extern void c5_simulation2_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
