/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation2_sfun.h"
#include "c2_simulation2.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simulation2_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c2_debug_family_names[6] = { "quad", "nargin", "nargout",
  "u", "omega", "omegaDot" };

/* Function Declarations */
static void initialize_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance);
static void initialize_params_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance);
static void enable_c2_simulation2(SFc2_simulation2InstanceStruct *chartInstance);
static void disable_c2_simulation2(SFc2_simulation2InstanceStruct *chartInstance);
static void c2_update_debugger_state_c2_simulation2
  (SFc2_simulation2InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c2_simulation2
  (SFc2_simulation2InstanceStruct *chartInstance);
static void set_sim_state_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_st);
static void finalize_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance);
static void sf_gateway_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance);
static void mdl_start_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance);
static void initSimStructsc2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber);
static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData);
static void c2_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_omegaDot, const char_T *c2_identifier, real_T c2_y[3]);
static void c2_b_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[3]);
static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static real_T c2_c_emlrt_marshallIn(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_info_helper(const mxArray **c2_info);
static const mxArray *c2_emlrt_marshallOut(const char * c2_b_u);
static const mxArray *c2_b_emlrt_marshallOut(const uint32_T c2_b_u);
static void c2_eml_scalar_eg(SFc2_simulation2InstanceStruct *chartInstance);
static void c2_eml_xgemm(SFc2_simulation2InstanceStruct *chartInstance, real_T
  c2_A[9], real_T c2_B[3], real_T c2_C[3], real_T c2_b_C[3]);
static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static int32_T c2_d_emlrt_marshallIn(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_e_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId,
  c2_sdz1UxT0ktEtEZMZuZvzEsH *c2_y);
static void c2_f_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[8]);
static void c2_g_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[4]);
static void c2_h_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[12]);
static void c2_i_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[9]);
static void c2_j_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[16]);
static void c2_k_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[2]);
static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static uint8_T c2_l_emlrt_marshallIn(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_simulation2, const char_T
  *c2_identifier);
static uint8_T c2_m_emlrt_marshallIn(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_b_eml_xgemm(SFc2_simulation2InstanceStruct *chartInstance, real_T
  c2_A[9], real_T c2_B[3], real_T c2_C[3]);
static void init_dsm_address_info(SFc2_simulation2InstanceStruct *chartInstance);
static void init_simulink_io_address(SFc2_simulation2InstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance)
{
  chartInstance->c2_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c2_is_active_c2_simulation2 = 0U;
}

static void initialize_params_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance)
{
  const mxArray *c2_m0 = NULL;
  const mxArray *c2_mxField;
  c2_sdz1UxT0ktEtEZMZuZvzEsH c2_r0;
  c2_m0 = sf_mex_get_sfun_param(chartInstance->S, 0, 1);
  c2_mxField = sf_mex_getfield(c2_m0, "thrustLimits", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.thrustLimits, 1, 0,
                      0U, 1, 0U, 2, 2, 4);
  c2_mxField = sf_mex_getfield(c2_m0, "thrustRateLimits", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.thrustRateLimits, 1,
                      0, 0U, 1, 0U, 2, 1, 4);
  c2_mxField = sf_mex_getfield(c2_m0, "thrustDir", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.thrustDir, 1, 0, 0U,
                      1, 0U, 2, 3, 4);
  c2_mxField = sf_mex_getfield(c2_m0, "L", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.L, 1, 0, 0U, 1, 0U,
                      2, 3, 4);
  c2_mxField = sf_mex_getfield(c2_m0, "nL", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.nL, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "rad", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.rad, 1, 0, 0U, 0,
                      0U, 0);
  c2_mxField = sf_mex_getfield(c2_m0, "bladeRad", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.bladeRad, 1, 0, 0U,
                      0, 0U, 0);
  c2_mxField = sf_mex_getfield(c2_m0, "mass", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.mass, 1, 0, 0U, 0,
                      0U, 0);
  c2_mxField = sf_mex_getfield(c2_m0, "I", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.I, 1, 0, 0U, 1, 0U,
                      2, 3, 3);
  c2_mxField = sf_mex_getfield(c2_m0, "km", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.km, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "kf", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.kf, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "K", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.K, 1, 0, 0U, 1, 0U,
                      2, 4, 4);
  c2_mxField = sf_mex_getfield(c2_m0, "velocityGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.velocityGain, 1, 0,
                      0U, 0, 0U, 0);
  c2_mxField = sf_mex_getfield(c2_m0, "vGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.vGain, 1, 0, 0U, 1,
                      0U, 1, 2);
  c2_mxField = sf_mex_getfield(c2_m0, "zGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.zGain, 1, 0, 0U, 0,
                      0U, 0);
  c2_mxField = sf_mex_getfield(c2_m0, "thetaGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.thetaGain, 1, 0, 0U,
                      1, 0U, 1, 3);
  c2_mxField = sf_mex_getfield(c2_m0, "omegaGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), c2_r0.omegaGain, 1, 0, 0U,
                      1, 0U, 1, 3);
  c2_mxField = sf_mex_getfield(c2_m0, "Ts2", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c2_mxField), &c2_r0.Ts2, 1, 0, 0U, 0,
                      0U, 0);
  sf_mex_destroy(&c2_m0);
  chartInstance->c2_quad = c2_r0;
}

static void enable_c2_simulation2(SFc2_simulation2InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c2_simulation2(SFc2_simulation2InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c2_update_debugger_state_c2_simulation2
  (SFc2_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c2_simulation2
  (SFc2_simulation2InstanceStruct *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_y = NULL;
  int32_T c2_i0;
  real_T c2_b_u[3];
  const mxArray *c2_b_y = NULL;
  uint8_T c2_hoistedGlobal;
  uint8_T c2_c_u;
  const mxArray *c2_c_y = NULL;
  c2_st = NULL;
  c2_st = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createcellmatrix(2, 1), false);
  for (c2_i0 = 0; c2_i0 < 3; c2_i0++) {
    c2_b_u[c2_i0] = (*chartInstance->c2_omegaDot)[c2_i0];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_b_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c2_y, 0, c2_b_y);
  c2_hoistedGlobal = chartInstance->c2_is_active_c2_simulation2;
  c2_c_u = c2_hoistedGlobal;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_c_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 1, c2_c_y);
  sf_mex_assign(&c2_st, c2_y, false);
  return c2_st;
}

static void set_sim_state_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_st)
{
  const mxArray *c2_b_u;
  real_T c2_dv0[3];
  int32_T c2_i1;
  chartInstance->c2_doneDoubleBufferReInit = true;
  c2_b_u = sf_mex_dup(c2_st);
  c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_b_u, 0)),
                      "omegaDot", c2_dv0);
  for (c2_i1 = 0; c2_i1 < 3; c2_i1++) {
    (*chartInstance->c2_omegaDot)[c2_i1] = c2_dv0[c2_i1];
  }

  chartInstance->c2_is_active_c2_simulation2 = c2_l_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c2_b_u, 1)),
     "is_active_c2_simulation2");
  sf_mex_destroy(&c2_b_u);
  c2_update_debugger_state_c2_simulation2(chartInstance);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance)
{
  int32_T c2_i2;
  int32_T c2_i3;
  int32_T c2_i4;
  real_T c2_b_u[3];
  int32_T c2_i5;
  real_T c2_b_omega[3];
  uint32_T c2_debug_family_var_map[6];
  c2_sdz1UxT0ktEtEZMZuZvzEsH c2_b_quad;
  real_T c2_nargin = 3.0;
  real_T c2_nargout = 1.0;
  real_T c2_b_omegaDot[3];
  static c2_sdz1UxT0ktEtEZMZuZvzEsH c2_r1 = { { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
      0.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 1.0 }, { 0.0, 0.2, 0.0, -0.2, 0.0, 0.0, 0.0, -0.2, 0.0,
      0.2, 0.0, 0.0 }, 0.2, 0.040000000000000008, 0.080000000000000016, 8.0, {
      10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 15.0 }, 11.0, 28.0, { 28.0, 0.0,
      -5.6000000000000005, 11.0, 28.0, 5.6000000000000005, 0.0, -11.0, 28.0, 0.0,
      5.6000000000000005, 11.0, 28.0, -5.6000000000000005, 0.0, -11.0 }, 10.0, {
      0.1, 0.1 }, 1.0, { 1.0, 1.0, 1.0 }, { 0.2, 0.2, 0.2 }, 0.3 };

  int32_T c2_i6;
  real_T c2_b[3];
  int32_T c2_i7;
  real_T c2_y[3];
  int32_T c2_i8;
  static real_T c2_a[9] = { 10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 15.0 };

  real_T c2_b_a[9];
  int32_T c2_i9;
  real_T c2_b_b[3];
  int32_T c2_i10;
  real_T c2_c1;
  real_T c2_c2;
  real_T c2_c3;
  int32_T c2_i11;
  int32_T c2_i12;
  int32_T c2_i13;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  for (c2_i2 = 0; c2_i2 < 3; c2_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_u)[c2_i2], 0U);
  }

  for (c2_i3 = 0; c2_i3 < 3; c2_i3++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_omega)[c2_i3], 1U);
  }

  chartInstance->c2_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  for (c2_i4 = 0; c2_i4 < 3; c2_i4++) {
    c2_b_u[c2_i4] = (*chartInstance->c2_u)[c2_i4];
  }

  for (c2_i5 = 0; c2_i5 < 3; c2_i5++) {
    c2_b_omega[c2_i5] = (*chartInstance->c2_omega)[c2_i5];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c2_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_b_quad, 0U, c2_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 1U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 2U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_u, 3U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_omega, 4U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_b_omegaDot, 5U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  c2_b_quad = c2_r1;
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 4);
  for (c2_i6 = 0; c2_i6 < 3; c2_i6++) {
    c2_b[c2_i6] = c2_b_omega[c2_i6];
  }

  c2_eml_scalar_eg(chartInstance);
  c2_eml_scalar_eg(chartInstance);
  for (c2_i7 = 0; c2_i7 < 3; c2_i7++) {
    c2_y[c2_i7] = 0.0;
  }

  for (c2_i8 = 0; c2_i8 < 9; c2_i8++) {
    c2_b_a[c2_i8] = c2_a[c2_i8];
  }

  for (c2_i9 = 0; c2_i9 < 3; c2_i9++) {
    c2_b_b[c2_i9] = c2_b[c2_i9];
  }

  c2_b_eml_xgemm(chartInstance, c2_b_a, c2_b_b, c2_y);
  for (c2_i10 = 0; c2_i10 < 3; c2_i10++) {
    c2_b[c2_i10] = c2_b_omega[c2_i10];
  }

  c2_c1 = c2_b[1] * c2_y[2] - c2_b[2] * c2_y[1];
  c2_c2 = c2_b[2] * c2_y[0] - c2_b[0] * c2_y[2];
  c2_c3 = c2_b[0] * c2_y[1] - c2_b[1] * c2_y[0];
  c2_b[0] = c2_c1;
  c2_b[1] = c2_c2;
  c2_b[2] = c2_c3;
  for (c2_i11 = 0; c2_i11 < 3; c2_i11++) {
    c2_b_omegaDot[c2_i11] = -c2_b[c2_i11] + c2_b_u[c2_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, -4);
  _SFD_SYMBOL_SCOPE_POP();
  for (c2_i12 = 0; c2_i12 < 3; c2_i12++) {
    (*chartInstance->c2_omegaDot)[c2_i12] = c2_b_omegaDot[c2_i12];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulation2MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c2_i13 = 0; c2_i13 < 3; c2_i13++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_omegaDot)[c2_i13], 3U);
  }
}

static void mdl_start_c2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void initSimStructsc2_simulation2(SFc2_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber)
{
  (void)c2_machineNumber;
  (void)c2_chartNumber;
  (void)c2_instanceNumber;
}

static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i14;
  real_T c2_b_inData[3];
  int32_T c2_i15;
  real_T c2_b_u[3];
  const mxArray *c2_y = NULL;
  SFc2_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc2_simulation2InstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i14 = 0; c2_i14 < 3; c2_i14++) {
    c2_b_inData[c2_i14] = (*(real_T (*)[3])c2_inData)[c2_i14];
  }

  for (c2_i15 = 0; c2_i15 < 3; c2_i15++) {
    c2_b_u[c2_i15] = c2_b_inData[c2_i15];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_b_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_omegaDot, const char_T *c2_identifier, real_T c2_y[3])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_omegaDot), &c2_thisId,
                        c2_y);
  sf_mex_destroy(&c2_b_omegaDot);
}

static void c2_b_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[3])
{
  real_T c2_dv1[3];
  int32_T c2_i16;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), c2_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c2_i16 = 0; c2_i16 < 3; c2_i16++) {
    c2_y[c2_i16] = c2_dv1[c2_i16];
  }

  sf_mex_destroy(&c2_b_u);
}

static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_omegaDot;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[3];
  int32_T c2_i17;
  SFc2_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc2_simulation2InstanceStruct *)chartInstanceVoid;
  c2_b_omegaDot = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_omegaDot), &c2_thisId,
                        c2_y);
  sf_mex_destroy(&c2_b_omegaDot);
  for (c2_i17 = 0; c2_i17 < 3; c2_i17++) {
    (*(real_T (*)[3])c2_outData)[c2_i17] = c2_y[c2_i17];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  real_T c2_b_u;
  const mxArray *c2_y = NULL;
  SFc2_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc2_simulation2InstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_b_u = *(real_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static real_T c2_c_emlrt_marshallIn(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), &c2_d0, 1, 0, 0U, 0, 0U, 0);
  c2_y = c2_d0;
  sf_mex_destroy(&c2_b_u);
  return c2_y;
}

static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_nargout;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y;
  SFc2_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc2_simulation2InstanceStruct *)chartInstanceVoid;
  c2_nargout = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_nargout), &c2_thisId);
  sf_mex_destroy(&c2_nargout);
  *(real_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  c2_sdz1UxT0ktEtEZMZuZvzEsH c2_b_u;
  const mxArray *c2_y = NULL;
  int32_T c2_i18;
  real_T c2_c_u[8];
  const mxArray *c2_b_y = NULL;
  int32_T c2_i19;
  real_T c2_d_u[4];
  const mxArray *c2_c_y = NULL;
  int32_T c2_i20;
  real_T c2_e_u[12];
  const mxArray *c2_d_y = NULL;
  int32_T c2_i21;
  real_T c2_f_u[12];
  const mxArray *c2_e_y = NULL;
  real_T c2_g_u;
  const mxArray *c2_f_y = NULL;
  real_T c2_h_u;
  const mxArray *c2_g_y = NULL;
  real_T c2_i_u;
  const mxArray *c2_h_y = NULL;
  real_T c2_j_u;
  const mxArray *c2_i_y = NULL;
  int32_T c2_i22;
  real_T c2_k_u[9];
  const mxArray *c2_j_y = NULL;
  real_T c2_l_u;
  const mxArray *c2_k_y = NULL;
  real_T c2_m_u;
  const mxArray *c2_l_y = NULL;
  int32_T c2_i23;
  real_T c2_n_u[16];
  const mxArray *c2_m_y = NULL;
  real_T c2_o_u;
  const mxArray *c2_n_y = NULL;
  int32_T c2_i24;
  real_T c2_p_u[2];
  const mxArray *c2_o_y = NULL;
  real_T c2_q_u;
  const mxArray *c2_p_y = NULL;
  int32_T c2_i25;
  real_T c2_r_u[3];
  const mxArray *c2_q_y = NULL;
  int32_T c2_i26;
  real_T c2_s_u[3];
  const mxArray *c2_r_y = NULL;
  real_T c2_t_u;
  const mxArray *c2_s_y = NULL;
  SFc2_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc2_simulation2InstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_b_u = *(c2_sdz1UxT0ktEtEZMZuZvzEsH *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  for (c2_i18 = 0; c2_i18 < 8; c2_i18++) {
    c2_c_u[c2_i18] = c2_b_u.thrustLimits[c2_i18];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 2, 2, 4),
                false);
  sf_mex_addfield(c2_y, c2_b_y, "thrustLimits", "thrustLimits", 0);
  for (c2_i19 = 0; c2_i19 < 4; c2_i19++) {
    c2_d_u[c2_i19] = c2_b_u.thrustRateLimits[c2_i19];
  }

  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", c2_d_u, 0, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_addfield(c2_y, c2_c_y, "thrustRateLimits", "thrustRateLimits", 0);
  for (c2_i20 = 0; c2_i20 < 12; c2_i20++) {
    c2_e_u[c2_i20] = c2_b_u.thrustDir[c2_i20];
  }

  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", c2_e_u, 0, 0U, 1U, 0U, 2, 3, 4),
                false);
  sf_mex_addfield(c2_y, c2_d_y, "thrustDir", "thrustDir", 0);
  for (c2_i21 = 0; c2_i21 < 12; c2_i21++) {
    c2_f_u[c2_i21] = c2_b_u.L[c2_i21];
  }

  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", c2_f_u, 0, 0U, 1U, 0U, 2, 3, 4),
                false);
  sf_mex_addfield(c2_y, c2_e_y, "L", "L", 0);
  c2_g_u = c2_b_u.nL;
  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", &c2_g_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_f_y, "nL", "nL", 0);
  c2_h_u = c2_b_u.rad;
  c2_g_y = NULL;
  sf_mex_assign(&c2_g_y, sf_mex_create("y", &c2_h_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_g_y, "rad", "rad", 0);
  c2_i_u = c2_b_u.bladeRad;
  c2_h_y = NULL;
  sf_mex_assign(&c2_h_y, sf_mex_create("y", &c2_i_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_h_y, "bladeRad", "bladeRad", 0);
  c2_j_u = c2_b_u.mass;
  c2_i_y = NULL;
  sf_mex_assign(&c2_i_y, sf_mex_create("y", &c2_j_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_i_y, "mass", "mass", 0);
  for (c2_i22 = 0; c2_i22 < 9; c2_i22++) {
    c2_k_u[c2_i22] = c2_b_u.I[c2_i22];
  }

  c2_j_y = NULL;
  sf_mex_assign(&c2_j_y, sf_mex_create("y", c2_k_u, 0, 0U, 1U, 0U, 2, 3, 3),
                false);
  sf_mex_addfield(c2_y, c2_j_y, "I", "I", 0);
  c2_l_u = c2_b_u.km;
  c2_k_y = NULL;
  sf_mex_assign(&c2_k_y, sf_mex_create("y", &c2_l_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_k_y, "km", "km", 0);
  c2_m_u = c2_b_u.kf;
  c2_l_y = NULL;
  sf_mex_assign(&c2_l_y, sf_mex_create("y", &c2_m_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_l_y, "kf", "kf", 0);
  for (c2_i23 = 0; c2_i23 < 16; c2_i23++) {
    c2_n_u[c2_i23] = c2_b_u.K[c2_i23];
  }

  c2_m_y = NULL;
  sf_mex_assign(&c2_m_y, sf_mex_create("y", c2_n_u, 0, 0U, 1U, 0U, 2, 4, 4),
                false);
  sf_mex_addfield(c2_y, c2_m_y, "K", "K", 0);
  c2_o_u = c2_b_u.velocityGain;
  c2_n_y = NULL;
  sf_mex_assign(&c2_n_y, sf_mex_create("y", &c2_o_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_n_y, "velocityGain", "velocityGain", 0);
  for (c2_i24 = 0; c2_i24 < 2; c2_i24++) {
    c2_p_u[c2_i24] = c2_b_u.vGain[c2_i24];
  }

  c2_o_y = NULL;
  sf_mex_assign(&c2_o_y, sf_mex_create("y", c2_p_u, 0, 0U, 1U, 0U, 1, 2), false);
  sf_mex_addfield(c2_y, c2_o_y, "vGain", "vGain", 0);
  c2_q_u = c2_b_u.zGain;
  c2_p_y = NULL;
  sf_mex_assign(&c2_p_y, sf_mex_create("y", &c2_q_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_p_y, "zGain", "zGain", 0);
  for (c2_i25 = 0; c2_i25 < 3; c2_i25++) {
    c2_r_u[c2_i25] = c2_b_u.thetaGain[c2_i25];
  }

  c2_q_y = NULL;
  sf_mex_assign(&c2_q_y, sf_mex_create("y", c2_r_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_addfield(c2_y, c2_q_y, "thetaGain", "thetaGain", 0);
  for (c2_i26 = 0; c2_i26 < 3; c2_i26++) {
    c2_s_u[c2_i26] = c2_b_u.omegaGain[c2_i26];
  }

  c2_r_y = NULL;
  sf_mex_assign(&c2_r_y, sf_mex_create("y", c2_s_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_addfield(c2_y, c2_r_y, "omegaGain", "omegaGain", 0);
  c2_t_u = c2_b_u.Ts2;
  c2_s_y = NULL;
  sf_mex_assign(&c2_s_y, sf_mex_create("y", &c2_t_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_s_y, "Ts2", "Ts2", 0);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

const mxArray *sf_c2_simulation2_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  sf_mex_assign(&c2_nameCaptureInfo, sf_mex_createstruct("structure", 2, 14, 1),
                false);
  c2_info_helper(&c2_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c2_nameCaptureInfo);
  return c2_nameCaptureInfo;
}

static void c2_info_helper(const mxArray **c2_info)
{
  const mxArray *c2_rhs0 = NULL;
  const mxArray *c2_lhs0 = NULL;
  const mxArray *c2_rhs1 = NULL;
  const mxArray *c2_lhs1 = NULL;
  const mxArray *c2_rhs2 = NULL;
  const mxArray *c2_lhs2 = NULL;
  const mxArray *c2_rhs3 = NULL;
  const mxArray *c2_lhs3 = NULL;
  const mxArray *c2_rhs4 = NULL;
  const mxArray *c2_lhs4 = NULL;
  const mxArray *c2_rhs5 = NULL;
  const mxArray *c2_lhs5 = NULL;
  const mxArray *c2_rhs6 = NULL;
  const mxArray *c2_lhs6 = NULL;
  const mxArray *c2_rhs7 = NULL;
  const mxArray *c2_lhs7 = NULL;
  const mxArray *c2_rhs8 = NULL;
  const mxArray *c2_lhs8 = NULL;
  const mxArray *c2_rhs9 = NULL;
  const mxArray *c2_lhs9 = NULL;
  const mxArray *c2_rhs10 = NULL;
  const mxArray *c2_lhs10 = NULL;
  const mxArray *c2_rhs11 = NULL;
  const mxArray *c2_lhs11 = NULL;
  const mxArray *c2_rhs12 = NULL;
  const mxArray *c2_lhs12 = NULL;
  const mxArray *c2_rhs13 = NULL;
  const mxArray *c2_lhs13 = NULL;
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c2_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c2_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c2_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c2_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c2_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_xgemm"), "name", "name", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375980690U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c2_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c2_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c2_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c2_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c2_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c2_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c2_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c2_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("cross"), "name", "name", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/cross.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286818842U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c2_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs13), "lhs", "lhs",
                  13);
  sf_mex_destroy(&c2_rhs0);
  sf_mex_destroy(&c2_lhs0);
  sf_mex_destroy(&c2_rhs1);
  sf_mex_destroy(&c2_lhs1);
  sf_mex_destroy(&c2_rhs2);
  sf_mex_destroy(&c2_lhs2);
  sf_mex_destroy(&c2_rhs3);
  sf_mex_destroy(&c2_lhs3);
  sf_mex_destroy(&c2_rhs4);
  sf_mex_destroy(&c2_lhs4);
  sf_mex_destroy(&c2_rhs5);
  sf_mex_destroy(&c2_lhs5);
  sf_mex_destroy(&c2_rhs6);
  sf_mex_destroy(&c2_lhs6);
  sf_mex_destroy(&c2_rhs7);
  sf_mex_destroy(&c2_lhs7);
  sf_mex_destroy(&c2_rhs8);
  sf_mex_destroy(&c2_lhs8);
  sf_mex_destroy(&c2_rhs9);
  sf_mex_destroy(&c2_lhs9);
  sf_mex_destroy(&c2_rhs10);
  sf_mex_destroy(&c2_lhs10);
  sf_mex_destroy(&c2_rhs11);
  sf_mex_destroy(&c2_lhs11);
  sf_mex_destroy(&c2_rhs12);
  sf_mex_destroy(&c2_lhs12);
  sf_mex_destroy(&c2_rhs13);
  sf_mex_destroy(&c2_lhs13);
}

static const mxArray *c2_emlrt_marshallOut(const char * c2_b_u)
{
  const mxArray *c2_y = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_b_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c2_b_u)), false);
  return c2_y;
}

static const mxArray *c2_b_emlrt_marshallOut(const uint32_T c2_b_u)
{
  const mxArray *c2_y = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_b_u, 7, 0U, 0U, 0U, 0), false);
  return c2_y;
}

static void c2_eml_scalar_eg(SFc2_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c2_eml_xgemm(SFc2_simulation2InstanceStruct *chartInstance, real_T
  c2_A[9], real_T c2_B[3], real_T c2_C[3], real_T c2_b_C[3])
{
  int32_T c2_i27;
  int32_T c2_i28;
  real_T c2_b_A[9];
  int32_T c2_i29;
  real_T c2_b_B[3];
  for (c2_i27 = 0; c2_i27 < 3; c2_i27++) {
    c2_b_C[c2_i27] = c2_C[c2_i27];
  }

  for (c2_i28 = 0; c2_i28 < 9; c2_i28++) {
    c2_b_A[c2_i28] = c2_A[c2_i28];
  }

  for (c2_i29 = 0; c2_i29 < 3; c2_i29++) {
    c2_b_B[c2_i29] = c2_B[c2_i29];
  }

  c2_b_eml_xgemm(chartInstance, c2_b_A, c2_b_B, c2_b_C);
}

static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_b_u;
  const mxArray *c2_y = NULL;
  SFc2_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc2_simulation2InstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_b_u = *(int32_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_b_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static int32_T c2_d_emlrt_marshallIn(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId)
{
  int32_T c2_y;
  int32_T c2_i30;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), &c2_i30, 1, 6, 0U, 0, 0U, 0);
  c2_y = c2_i30;
  sf_mex_destroy(&c2_b_u);
  return c2_y;
}

static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_sfEvent;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  int32_T c2_y;
  SFc2_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc2_simulation2InstanceStruct *)chartInstanceVoid;
  c2_b_sfEvent = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_sfEvent),
    &c2_thisId);
  sf_mex_destroy(&c2_b_sfEvent);
  *(int32_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static void c2_e_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId,
  c2_sdz1UxT0ktEtEZMZuZvzEsH *c2_y)
{
  emlrtMsgIdentifier c2_thisId;
  static const char * c2_fieldNames[18] = { "thrustLimits", "thrustRateLimits",
    "thrustDir", "L", "nL", "rad", "bladeRad", "mass", "I", "km", "kf", "K",
    "velocityGain", "vGain", "zGain", "thetaGain", "omegaGain", "Ts2" };

  c2_thisId.fParent = c2_parentId;
  sf_mex_check_struct(c2_parentId, c2_b_u, 18, c2_fieldNames, 0U, NULL);
  c2_thisId.fIdentifier = "thrustLimits";
  c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u,
    "thrustLimits", "thrustLimits", 0)), &c2_thisId, c2_y->thrustLimits);
  c2_thisId.fIdentifier = "thrustRateLimits";
  c2_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u,
    "thrustRateLimits", "thrustRateLimits", 0)), &c2_thisId,
                        c2_y->thrustRateLimits);
  c2_thisId.fIdentifier = "thrustDir";
  c2_h_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u,
    "thrustDir", "thrustDir", 0)), &c2_thisId, c2_y->thrustDir);
  c2_thisId.fIdentifier = "L";
  c2_h_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u, "L",
    "L", 0)), &c2_thisId, c2_y->L);
  c2_thisId.fIdentifier = "nL";
  c2_y->nL = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_b_u, "nL", "nL", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rad";
  c2_y->rad = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_b_u, "rad", "rad", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "bladeRad";
  c2_y->bladeRad = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c2_b_u, "bladeRad", "bladeRad", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "mass";
  c2_y->mass = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_b_u, "mass", "mass", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "I";
  c2_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u, "I",
    "I", 0)), &c2_thisId, c2_y->I);
  c2_thisId.fIdentifier = "km";
  c2_y->km = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_b_u, "km", "km", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "kf";
  c2_y->kf = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_b_u, "kf", "kf", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "K";
  c2_j_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u, "K",
    "K", 0)), &c2_thisId, c2_y->K);
  c2_thisId.fIdentifier = "velocityGain";
  c2_y->velocityGain = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c2_b_u, "velocityGain", "velocityGain", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "vGain";
  c2_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u,
    "vGain", "vGain", 0)), &c2_thisId, c2_y->vGain);
  c2_thisId.fIdentifier = "zGain";
  c2_y->zGain = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_b_u, "zGain", "zGain", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "thetaGain";
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u,
    "thetaGain", "thetaGain", 0)), &c2_thisId, c2_y->thetaGain);
  c2_thisId.fIdentifier = "omegaGain";
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_b_u,
    "omegaGain", "omegaGain", 0)), &c2_thisId, c2_y->omegaGain);
  c2_thisId.fIdentifier = "Ts2";
  c2_y->Ts2 = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_b_u, "Ts2", "Ts2", 0)), &c2_thisId);
  sf_mex_destroy(&c2_b_u);
}

static void c2_f_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[8])
{
  real_T c2_dv2[8];
  int32_T c2_i31;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), c2_dv2, 1, 0, 0U, 1, 0U, 2, 2,
                4);
  for (c2_i31 = 0; c2_i31 < 8; c2_i31++) {
    c2_y[c2_i31] = c2_dv2[c2_i31];
  }

  sf_mex_destroy(&c2_b_u);
}

static void c2_g_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[4])
{
  real_T c2_dv3[4];
  int32_T c2_i32;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), c2_dv3, 1, 0, 0U, 1, 0U, 2, 1,
                4);
  for (c2_i32 = 0; c2_i32 < 4; c2_i32++) {
    c2_y[c2_i32] = c2_dv3[c2_i32];
  }

  sf_mex_destroy(&c2_b_u);
}

static void c2_h_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[12])
{
  real_T c2_dv4[12];
  int32_T c2_i33;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), c2_dv4, 1, 0, 0U, 1, 0U, 2, 3,
                4);
  for (c2_i33 = 0; c2_i33 < 12; c2_i33++) {
    c2_y[c2_i33] = c2_dv4[c2_i33];
  }

  sf_mex_destroy(&c2_b_u);
}

static void c2_i_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[9])
{
  real_T c2_dv5[9];
  int32_T c2_i34;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), c2_dv5, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c2_i34 = 0; c2_i34 < 9; c2_i34++) {
    c2_y[c2_i34] = c2_dv5[c2_i34];
  }

  sf_mex_destroy(&c2_b_u);
}

static void c2_j_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[16])
{
  real_T c2_dv6[16];
  int32_T c2_i35;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), c2_dv6, 1, 0, 0U, 1, 0U, 2, 4,
                4);
  for (c2_i35 = 0; c2_i35 < 16; c2_i35++) {
    c2_y[c2_i35] = c2_dv6[c2_i35];
  }

  sf_mex_destroy(&c2_b_u);
}

static void c2_k_emlrt_marshallIn(SFc2_simulation2InstanceStruct *chartInstance,
  const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[2])
{
  real_T c2_dv7[2];
  int32_T c2_i36;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), c2_dv7, 1, 0, 0U, 1, 0U, 1, 2);
  for (c2_i36 = 0; c2_i36 < 2; c2_i36++) {
    c2_y[c2_i36] = c2_dv7[c2_i36];
  }

  sf_mex_destroy(&c2_b_u);
}

static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_quad;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  c2_sdz1UxT0ktEtEZMZuZvzEsH c2_y;
  SFc2_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc2_simulation2InstanceStruct *)chartInstanceVoid;
  c2_b_quad = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_quad), &c2_thisId, &c2_y);
  sf_mex_destroy(&c2_b_quad);
  *(c2_sdz1UxT0ktEtEZMZuZvzEsH *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static uint8_T c2_l_emlrt_marshallIn(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_simulation2, const char_T
  *c2_identifier)
{
  uint8_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_m_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c2_b_is_active_c2_simulation2), &c2_thisId);
  sf_mex_destroy(&c2_b_is_active_c2_simulation2);
  return c2_y;
}

static uint8_T c2_m_emlrt_marshallIn(SFc2_simulation2InstanceStruct
  *chartInstance, const mxArray *c2_b_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_y;
  uint8_T c2_u0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_b_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_y = c2_u0;
  sf_mex_destroy(&c2_b_u);
  return c2_y;
}

static void c2_b_eml_xgemm(SFc2_simulation2InstanceStruct *chartInstance, real_T
  c2_A[9], real_T c2_B[3], real_T c2_C[3])
{
  int32_T c2_i37;
  int32_T c2_i38;
  int32_T c2_i39;
  (void)chartInstance;
  for (c2_i37 = 0; c2_i37 < 3; c2_i37++) {
    c2_C[c2_i37] = 0.0;
    c2_i38 = 0;
    for (c2_i39 = 0; c2_i39 < 3; c2_i39++) {
      c2_C[c2_i37] += c2_A[c2_i38 + c2_i37] * c2_B[c2_i39];
      c2_i38 += 3;
    }
  }
}

static void init_dsm_address_info(SFc2_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc2_simulation2InstanceStruct
  *chartInstance)
{
  chartInstance->c2_u = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c2_omega = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c2_omegaDot = (real_T (*)[3])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c2_simulation2_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1088572093U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(200331150U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2301982599U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(976112478U);
}

mxArray* sf_c2_simulation2_get_post_codegen_info(void);
mxArray *sf_c2_simulation2_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("bbmbnBTwP4e4JQbxGg9zuD");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(13));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo = sf_c2_simulation2_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c2_simulation2_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c2_simulation2_jit_fallback_info(void)
{
  const char *infoFields[] = { "fallbackType", "fallbackReason",
    "incompatibleSymbol", };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 3, infoFields);
  mxArray *fallbackReason = mxCreateString("feature_off");
  mxArray *incompatibleSymbol = mxCreateString("");
  mxArray *fallbackType = mxCreateString("early");
  mxSetField(mxInfo, 0, infoFields[0], fallbackType);
  mxSetField(mxInfo, 0, infoFields[1], fallbackReason);
  mxSetField(mxInfo, 0, infoFields[2], incompatibleSymbol);
  return mxInfo;
}

mxArray *sf_c2_simulation2_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c2_simulation2_get_post_codegen_info(void)
{
  const char* fieldNames[] = { "exportedFunctionsUsedByThisChart",
    "exportedFunctionsChecksum" };

  mwSize dims[2] = { 1, 1 };

  mxArray* mxPostCodegenInfo = mxCreateStructArray(2, dims, sizeof(fieldNames)/
    sizeof(fieldNames[0]), fieldNames);

  {
    mxArray* mxExportedFunctionsChecksum = mxCreateString("");
    mwSize exp_dims[2] = { 0, 1 };

    mxArray* mxExportedFunctionsUsedByThisChart = mxCreateCellArray(2, exp_dims);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsUsedByThisChart",
               mxExportedFunctionsUsedByThisChart);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsChecksum",
               mxExportedFunctionsChecksum);
  }

  return mxPostCodegenInfo;
}

static const mxArray *sf_get_sim_state_info_c2_simulation2(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"omegaDot\",},{M[8],M[0],T\"is_active_c2_simulation2\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_simulation2_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_simulation2InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc2_simulation2InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulation2MachineNumber_,
           2,
           1,
           1,
           0,
           4,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation(_simulation2MachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_simulation2MachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _simulation2MachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"u");
          _SFD_SET_DATA_PROPS(1,1,1,0,"omega");
          _SFD_SET_DATA_PROPS(2,10,0,0,"quad");
          _SFD_SET_DATA_PROPS(3,2,0,1,"omegaDot");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,123);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_STRUCT,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)c2_d_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)
            c2_sf_marshallIn);
        }

        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c2_u);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c2_omega);
        _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c2_quad);
        _SFD_SET_DATA_VALUE_PTR(3U, *chartInstance->c2_omegaDot);
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _simulation2MachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "t60zV1JhPUSU8RgcLIXPBF";
}

static void sf_opaque_initialize_c2_simulation2(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc2_simulation2InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c2_simulation2((SFc2_simulation2InstanceStruct*)
    chartInstanceVar);
  initialize_c2_simulation2((SFc2_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c2_simulation2(void *chartInstanceVar)
{
  enable_c2_simulation2((SFc2_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c2_simulation2(void *chartInstanceVar)
{
  disable_c2_simulation2((SFc2_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c2_simulation2(void *chartInstanceVar)
{
  sf_gateway_c2_simulation2((SFc2_simulation2InstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c2_simulation2(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c2_simulation2((SFc2_simulation2InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c2_simulation2(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c2_simulation2((SFc2_simulation2InstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c2_simulation2(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_simulation2InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation2_optimization_info();
    }

    finalize_c2_simulation2((SFc2_simulation2InstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc2_simulation2((SFc2_simulation2InstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_simulation2(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c2_simulation2((SFc2_simulation2InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c2_simulation2(SimStruct *S)
{
  /* Actual parameters from chart:
     quad
   */
  const char_T *rtParamNames[] = { "quad" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation2_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,2,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,2);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,2,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,2,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1779554680U));
  ssSetChecksum1(S,(2923131833U));
  ssSetChecksum2(S,(1753268463U));
  ssSetChecksum3(S,(2041557501U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c2_simulation2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_simulation2(SimStruct *S)
{
  SFc2_simulation2InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc2_simulation2InstanceStruct *)utMalloc(sizeof
    (SFc2_simulation2InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc2_simulation2InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c2_simulation2;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c2_simulation2;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c2_simulation2;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c2_simulation2;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c2_simulation2;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c2_simulation2;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c2_simulation2;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c2_simulation2;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_simulation2;
  chartInstance->chartInfo.mdlStart = mdlStart_c2_simulation2;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c2_simulation2;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.callAtomicSubchartUserFcn = NULL;
  chartInstance->chartInfo.callAtomicSubchartAutoFcn = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->checksum = SF_RUNTIME_INFO_CHECKSUM;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  crtInfo->compiledInfo = NULL;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  init_simulink_io_address(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c2_simulation2_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_simulation2(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_simulation2(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_simulation2(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_simulation2_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
