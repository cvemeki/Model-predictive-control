#ifndef __c3_simulation1_h__
#define __c3_simulation1_h__

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

#ifndef typedef_c3_sdz1UxT0ktEtEZMZuZvzEsH
#define typedef_c3_sdz1UxT0ktEtEZMZuZvzEsH

typedef struct sdz1UxT0ktEtEZMZuZvzEsH c3_sdz1UxT0ktEtEZMZuZvzEsH;

#endif                                 /*typedef_c3_sdz1UxT0ktEtEZMZuZvzEsH*/

#ifndef typedef_SFc3_simulation1InstanceStruct
#define typedef_SFc3_simulation1InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c3_sfEvent;
  boolean_T c3_isStable;
  boolean_T c3_doneDoubleBufferReInit;
  uint8_T c3_is_active_c3_simulation1;
  c3_sdz1UxT0ktEtEZMZuZvzEsH c3_quad;
  real_T *c3_uTot;
  real_T *c3_dist;
  real_T (*c3_vDot)[3];
  real_T (*c3_theta)[3];
} SFc3_simulation1InstanceStruct;

#endif                                 /*typedef_SFc3_simulation1InstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c3_simulation1_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c3_simulation1_get_check_sum(mxArray *plhs[]);
extern void c3_simulation1_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
