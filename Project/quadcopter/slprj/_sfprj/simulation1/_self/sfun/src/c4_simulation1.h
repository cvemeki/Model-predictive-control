#ifndef __c4_simulation1_h__
#define __c4_simulation1_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef struct_sdz1UxT0ktEtEZMZuZvzEsH
#define struct_sdz1UxT0ktEtEZMZuZvzEsH

struct sdz1UxT0ktEtEZMZuZvzEsH
{
  real_T thrustLimits[8];
  real_T thrustRateLimits[4];
  real_T thrustDir[12];
  real_T L[12];
  real_T nL;
  real_T rad;
  real_T bladeRad;
  real_T mass;
  real_T I[9];
  real_T km;
  real_T kf;
  real_T K[16];
  real_T velocityGain;
  real_T vGain[2];
  real_T zGain;
  real_T thetaGain[3];
  real_T omegaGain[3];
  real_T Ts2;
};

#endif                                 /*struct_sdz1UxT0ktEtEZMZuZvzEsH*/

#ifndef typedef_c4_sdz1UxT0ktEtEZMZuZvzEsH
#define typedef_c4_sdz1UxT0ktEtEZMZuZvzEsH

typedef struct sdz1UxT0ktEtEZMZuZvzEsH c4_sdz1UxT0ktEtEZMZuZvzEsH;

#endif                                 /*typedef_c4_sdz1UxT0ktEtEZMZuZvzEsH*/

#ifndef typedef_SFc4_simulation1InstanceStruct
#define typedef_SFc4_simulation1InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c4_sfEvent;
  boolean_T c4_isStable;
  boolean_T c4_doneDoubleBufferReInit;
  uint8_T c4_is_active_c4_simulation1;
  c4_sdz1UxT0ktEtEZMZuZvzEsH c4_quad;
  uint32_T c4_method;
  boolean_T c4_method_not_empty;
  uint32_T c4_state;
  boolean_T c4_state_not_empty;
  uint32_T c4_b_state[2];
  boolean_T c4_b_state_not_empty;
  uint32_T c4_c_state[625];
  boolean_T c4_c_state_not_empty;
  real_T (*c4_u)[3];
  real_T *c4_dist;
  real_T (*c4_omega)[3];
  real_T (*c4_omegaDot)[3];
} SFc4_simulation1InstanceStruct;

#endif                                 /*typedef_SFc4_simulation1InstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c4_simulation1_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c4_simulation1_get_check_sum(mxArray *plhs[]);
extern void c4_simulation1_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
