#ifndef __c1_simulation2_h__
#define __c1_simulation2_h__

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

#ifndef typedef_c1_sdz1UxT0ktEtEZMZuZvzEsH
#define typedef_c1_sdz1UxT0ktEtEZMZuZvzEsH

typedef struct sdz1UxT0ktEtEZMZuZvzEsH c1_sdz1UxT0ktEtEZMZuZvzEsH;

#endif                                 /*typedef_c1_sdz1UxT0ktEtEZMZuZvzEsH*/

#ifndef typedef_SFc1_simulation2InstanceStruct
#define typedef_SFc1_simulation2InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_simulation2;
  c1_sdz1UxT0ktEtEZMZuZvzEsH c1_quad;
  real_T *c1_uTot;
  real_T (*c1_vDot)[3];
  real_T (*c1_theta)[3];
} SFc1_simulation2InstanceStruct;

#endif                                 /*typedef_SFc1_simulation2InstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c1_simulation2_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_simulation2_get_check_sum(mxArray *plhs[]);
extern void c1_simulation2_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
