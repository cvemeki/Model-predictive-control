/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation2_sfun.h"
#include "c5_simulation2.h"
#include "mwmathutil.h"
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
static const char * c5_debug_family_names[4] = { "nargin", "nargout", "dist",
  "d" };

/* Function Declarations */
static void initialize_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance);
static void initialize_params_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance);
static void enable_c5_simulation2(SFc5_simulation2InstanceStruct *chartInstance);
static void disable_c5_simulation2(SFc5_simulation2InstanceStruct *chartInstance);
static void c5_update_debugger_state_c5_simulation2
  (SFc5_simulation2InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c5_simulation2
  (SFc5_simulation2InstanceStruct *chartInstance);
static void set_sim_state_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_st);
static void finalize_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance);
static void sf_gateway_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance);
static void mdl_start_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance);
static void initSimStructsc5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber, uint32_T c5_instanceNumber);
static const mxArray *c5_sf_marshallOut(void *chartInstanceVoid, void *c5_inData);
static void c5_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_b_d, const char_T *c5_identifier, real_T c5_y[3]);
static void c5_b_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[3]);
static void c5_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static const mxArray *c5_b_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static real_T c5_c_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static void c5_info_helper(const mxArray **c5_info);
static const mxArray *c5_emlrt_marshallOut(const char * c5_u);
static const mxArray *c5_b_emlrt_marshallOut(const uint32_T c5_u);
static real_T c5_rand(SFc5_simulation2InstanceStruct *chartInstance);
static uint32_T c5_eml_rand_str2id(SFc5_simulation2InstanceStruct *chartInstance);
static uint32_T c5_b_eml_rand_str2id(SFc5_simulation2InstanceStruct
  *chartInstance);
static uint32_T c5_c_eml_rand_str2id(SFc5_simulation2InstanceStruct
  *chartInstance);
static void c5_eml_rand_mcg16807(SFc5_simulation2InstanceStruct *chartInstance,
  uint32_T c5_d_state, uint32_T *c5_e_state, real_T *c5_r);
static void c5_eml_rand_mt19937ar(SFc5_simulation2InstanceStruct *chartInstance,
  uint32_T c5_d_state[625]);
static void c5_twister_state_vector(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T c5_mt[625], real_T c5_seed, uint32_T c5_b_mt[625]);
static void c5_b_eml_rand_mt19937ar(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T c5_d_state[625], uint32_T c5_e_state[625], real_T
  *c5_r);
static void c5_eml_error(SFc5_simulation2InstanceStruct *chartInstance);
static void c5_b_rand(SFc5_simulation2InstanceStruct *chartInstance, real_T
                      c5_r[3]);
static const mxArray *c5_c_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static int32_T c5_d_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static uint32_T c5_e_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_b_method, const char_T *c5_identifier);
static uint32_T c5_f_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static uint32_T c5_g_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_d_state, const char_T *c5_identifier);
static uint32_T c5_h_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_i_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_d_state, const char_T *c5_identifier, uint32_T c5_y[625]);
static void c5_j_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, uint32_T c5_y[625]);
static void c5_k_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_d_state, const char_T *c5_identifier, uint32_T c5_y[2]);
static void c5_l_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, uint32_T c5_y[2]);
static uint8_T c5_m_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_b_is_active_c5_simulation2, const char_T
  *c5_identifier);
static uint8_T c5_n_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static real_T c5_b_eml_rand_mcg16807(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T *c5_d_state);
static void c5_b_twister_state_vector(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T c5_mt[625], real_T c5_seed);
static real_T c5_c_eml_rand_mt19937ar(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T c5_d_state[625]);
static void init_dsm_address_info(SFc5_simulation2InstanceStruct *chartInstance);
static void init_simulink_io_address(SFc5_simulation2InstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  chartInstance->c5_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c5_method_not_empty = false;
  chartInstance->c5_state_not_empty = false;
  chartInstance->c5_b_state_not_empty = false;
  chartInstance->c5_c_state_not_empty = false;
  chartInstance->c5_is_active_c5_simulation2 = 0U;
}

static void initialize_params_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c5_simulation2(SFc5_simulation2InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c5_simulation2(SFc5_simulation2InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c5_update_debugger_state_c5_simulation2
  (SFc5_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c5_simulation2
  (SFc5_simulation2InstanceStruct *chartInstance)
{
  const mxArray *c5_st;
  const mxArray *c5_y = NULL;
  int32_T c5_i0;
  real_T c5_u[3];
  const mxArray *c5_b_y = NULL;
  uint32_T c5_hoistedGlobal;
  uint32_T c5_b_u;
  const mxArray *c5_c_y = NULL;
  uint32_T c5_b_hoistedGlobal;
  uint32_T c5_c_u;
  const mxArray *c5_d_y = NULL;
  int32_T c5_i1;
  uint32_T c5_d_u[625];
  const mxArray *c5_e_y = NULL;
  int32_T c5_i2;
  uint32_T c5_e_u[2];
  const mxArray *c5_f_y = NULL;
  uint8_T c5_c_hoistedGlobal;
  uint8_T c5_f_u;
  const mxArray *c5_g_y = NULL;
  c5_st = NULL;
  c5_st = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_createcellmatrix(6, 1), false);
  for (c5_i0 = 0; c5_i0 < 3; c5_i0++) {
    c5_u[c5_i0] = (*chartInstance->c5_d)[c5_i0];
  }

  c5_b_y = NULL;
  sf_mex_assign(&c5_b_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c5_y, 0, c5_b_y);
  c5_hoistedGlobal = chartInstance->c5_method;
  c5_b_u = c5_hoistedGlobal;
  c5_c_y = NULL;
  if (!chartInstance->c5_method_not_empty) {
    sf_mex_assign(&c5_c_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c5_c_y, sf_mex_create("y", &c5_b_u, 7, 0U, 0U, 0U, 0), false);
  }

  sf_mex_setcell(c5_y, 1, c5_c_y);
  c5_b_hoistedGlobal = chartInstance->c5_state;
  c5_c_u = c5_b_hoistedGlobal;
  c5_d_y = NULL;
  if (!chartInstance->c5_state_not_empty) {
    sf_mex_assign(&c5_d_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c5_d_y, sf_mex_create("y", &c5_c_u, 7, 0U, 0U, 0U, 0), false);
  }

  sf_mex_setcell(c5_y, 2, c5_d_y);
  for (c5_i1 = 0; c5_i1 < 625; c5_i1++) {
    c5_d_u[c5_i1] = chartInstance->c5_c_state[c5_i1];
  }

  c5_e_y = NULL;
  if (!chartInstance->c5_c_state_not_empty) {
    sf_mex_assign(&c5_e_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c5_e_y, sf_mex_create("y", c5_d_u, 7, 0U, 1U, 0U, 1, 625),
                  false);
  }

  sf_mex_setcell(c5_y, 3, c5_e_y);
  for (c5_i2 = 0; c5_i2 < 2; c5_i2++) {
    c5_e_u[c5_i2] = chartInstance->c5_b_state[c5_i2];
  }

  c5_f_y = NULL;
  if (!chartInstance->c5_b_state_not_empty) {
    sf_mex_assign(&c5_f_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c5_f_y, sf_mex_create("y", c5_e_u, 7, 0U, 1U, 0U, 1, 2),
                  false);
  }

  sf_mex_setcell(c5_y, 4, c5_f_y);
  c5_c_hoistedGlobal = chartInstance->c5_is_active_c5_simulation2;
  c5_f_u = c5_c_hoistedGlobal;
  c5_g_y = NULL;
  sf_mex_assign(&c5_g_y, sf_mex_create("y", &c5_f_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c5_y, 5, c5_g_y);
  sf_mex_assign(&c5_st, c5_y, false);
  return c5_st;
}

static void set_sim_state_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_st)
{
  const mxArray *c5_u;
  real_T c5_dv0[3];
  int32_T c5_i3;
  uint32_T c5_uv0[625];
  int32_T c5_i4;
  uint32_T c5_uv1[2];
  int32_T c5_i5;
  chartInstance->c5_doneDoubleBufferReInit = true;
  c5_u = sf_mex_dup(c5_st);
  c5_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 0)), "d",
                      c5_dv0);
  for (c5_i3 = 0; c5_i3 < 3; c5_i3++) {
    (*chartInstance->c5_d)[c5_i3] = c5_dv0[c5_i3];
  }

  chartInstance->c5_method = c5_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c5_u, 1)), "method");
  chartInstance->c5_state = c5_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c5_u, 2)), "state");
  c5_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 3)),
                        "state", c5_uv0);
  for (c5_i4 = 0; c5_i4 < 625; c5_i4++) {
    chartInstance->c5_c_state[c5_i4] = c5_uv0[c5_i4];
  }

  c5_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 4)),
                        "state", c5_uv1);
  for (c5_i5 = 0; c5_i5 < 2; c5_i5++) {
    chartInstance->c5_b_state[c5_i5] = c5_uv1[c5_i5];
  }

  chartInstance->c5_is_active_c5_simulation2 = c5_m_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 5)),
     "is_active_c5_simulation2");
  sf_mex_destroy(&c5_u);
  c5_update_debugger_state_c5_simulation2(chartInstance);
  sf_mex_destroy(&c5_st);
}

static void finalize_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  real_T c5_hoistedGlobal;
  real_T c5_b_dist;
  uint32_T c5_debug_family_var_map[4];
  real_T c5_nargin = 1.0;
  real_T c5_nargout = 1.0;
  real_T c5_b_d[3];
  real_T c5_b;
  int32_T c5_i6;
  static real_T c5_a[3] = { 0.0, 0.0, 20.0 };

  real_T c5_y[3];
  real_T c5_b_b[3];
  int32_T c5_i7;
  int32_T c5_i8;
  int32_T c5_i9;
  int32_T c5_i10;
  int32_T c5_i11;
  int32_T c5_i12;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  chartInstance->c5_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  c5_hoistedGlobal = *chartInstance->c5_dist;
  c5_b_dist = c5_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c5_debug_family_names,
    c5_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_nargin, 0U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_nargout, 1U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_b_dist, 2U, c5_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_b_d, 3U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 3);
  if (CV_EML_IF(0, 1, 0, c5_b_dist != 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 4);
    c5_b = c5_rand(chartInstance) - 0.5;
    for (c5_i6 = 0; c5_i6 < 3; c5_i6++) {
      c5_y[c5_i6] = c5_a[c5_i6] * c5_b;
    }

    c5_b_rand(chartInstance, c5_b_b);
    for (c5_i7 = 0; c5_i7 < 3; c5_i7++) {
      c5_b_b[c5_i7] -= 0.5;
    }

    for (c5_i8 = 0; c5_i8 < 3; c5_i8++) {
      c5_b_b[c5_i8] *= 2.0;
    }

    for (c5_i9 = 0; c5_i9 < 3; c5_i9++) {
      c5_b_d[c5_i9] = (1.0 + c5_y[c5_i9]) + c5_b_b[c5_i9];
    }
  } else {
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 6);
    for (c5_i10 = 0; c5_i10 < 3; c5_i10++) {
      c5_b_d[c5_i10] = 0.0;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, -6);
  _SFD_SYMBOL_SCOPE_POP();
  for (c5_i11 = 0; c5_i11 < 3; c5_i11++) {
    (*chartInstance->c5_d)[c5_i11] = c5_b_d[c5_i11];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulation2MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c5_i12 = 0; c5_i12 < 3; c5_i12++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c5_d)[c5_i12], 0U);
  }

  _SFD_DATA_RANGE_CHECK(*chartInstance->c5_dist, 1U);
}

static void mdl_start_c5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void initSimStructsc5_simulation2(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber, uint32_T c5_instanceNumber)
{
  (void)c5_machineNumber;
  (void)c5_chartNumber;
  (void)c5_instanceNumber;
}

static const mxArray *c5_sf_marshallOut(void *chartInstanceVoid, void *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_i13;
  real_T c5_b_inData[3];
  int32_T c5_i14;
  real_T c5_u[3];
  const mxArray *c5_y = NULL;
  SFc5_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc5_simulation2InstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  for (c5_i13 = 0; c5_i13 < 3; c5_i13++) {
    c5_b_inData[c5_i13] = (*(real_T (*)[3])c5_inData)[c5_i13];
  }

  for (c5_i14 = 0; c5_i14 < 3; c5_i14++) {
    c5_u[c5_i14] = c5_b_inData[c5_i14];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static void c5_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_b_d, const char_T *c5_identifier, real_T c5_y[3])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_d), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_d);
}

static void c5_b_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[3])
{
  real_T c5_dv1[3];
  int32_T c5_i15;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c5_i15 = 0; c5_i15 < 3; c5_i15++) {
    c5_y[c5_i15] = c5_dv1[c5_i15];
  }

  sf_mex_destroy(&c5_u);
}

static void c5_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_b_d;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y[3];
  int32_T c5_i16;
  SFc5_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc5_simulation2InstanceStruct *)chartInstanceVoid;
  c5_b_d = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_d), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_d);
  for (c5_i16 = 0; c5_i16 < 3; c5_i16++) {
    (*(real_T (*)[3])c5_outData)[c5_i16] = c5_y[c5_i16];
  }

  sf_mex_destroy(&c5_mxArrayInData);
}

static const mxArray *c5_b_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  real_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc5_simulation2InstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(real_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static real_T c5_c_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d0;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d0, 1, 0, 0U, 0, 0U, 0);
  c5_y = c5_d0;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_nargout;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y;
  SFc5_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc5_simulation2InstanceStruct *)chartInstanceVoid;
  c5_nargout = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_nargout), &c5_thisId);
  sf_mex_destroy(&c5_nargout);
  *(real_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

const mxArray *sf_c5_simulation2_get_eml_resolved_functions_info(void)
{
  const mxArray *c5_nameCaptureInfo = NULL;
  c5_nameCaptureInfo = NULL;
  sf_mex_assign(&c5_nameCaptureInfo, sf_mex_createstruct("structure", 2, 27, 1),
                false);
  c5_info_helper(&c5_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c5_nameCaptureInfo);
  return c5_nameCaptureInfo;
}

static void c5_info_helper(const mxArray **c5_info)
{
  const mxArray *c5_rhs0 = NULL;
  const mxArray *c5_lhs0 = NULL;
  const mxArray *c5_rhs1 = NULL;
  const mxArray *c5_lhs1 = NULL;
  const mxArray *c5_rhs2 = NULL;
  const mxArray *c5_lhs2 = NULL;
  const mxArray *c5_rhs3 = NULL;
  const mxArray *c5_lhs3 = NULL;
  const mxArray *c5_rhs4 = NULL;
  const mxArray *c5_lhs4 = NULL;
  const mxArray *c5_rhs5 = NULL;
  const mxArray *c5_lhs5 = NULL;
  const mxArray *c5_rhs6 = NULL;
  const mxArray *c5_lhs6 = NULL;
  const mxArray *c5_rhs7 = NULL;
  const mxArray *c5_lhs7 = NULL;
  const mxArray *c5_rhs8 = NULL;
  const mxArray *c5_lhs8 = NULL;
  const mxArray *c5_rhs9 = NULL;
  const mxArray *c5_lhs9 = NULL;
  const mxArray *c5_rhs10 = NULL;
  const mxArray *c5_lhs10 = NULL;
  const mxArray *c5_rhs11 = NULL;
  const mxArray *c5_lhs11 = NULL;
  const mxArray *c5_rhs12 = NULL;
  const mxArray *c5_lhs12 = NULL;
  const mxArray *c5_rhs13 = NULL;
  const mxArray *c5_lhs13 = NULL;
  const mxArray *c5_rhs14 = NULL;
  const mxArray *c5_lhs14 = NULL;
  const mxArray *c5_rhs15 = NULL;
  const mxArray *c5_lhs15 = NULL;
  const mxArray *c5_rhs16 = NULL;
  const mxArray *c5_lhs16 = NULL;
  const mxArray *c5_rhs17 = NULL;
  const mxArray *c5_lhs17 = NULL;
  const mxArray *c5_rhs18 = NULL;
  const mxArray *c5_lhs18 = NULL;
  const mxArray *c5_rhs19 = NULL;
  const mxArray *c5_lhs19 = NULL;
  const mxArray *c5_rhs20 = NULL;
  const mxArray *c5_lhs20 = NULL;
  const mxArray *c5_rhs21 = NULL;
  const mxArray *c5_lhs21 = NULL;
  const mxArray *c5_rhs22 = NULL;
  const mxArray *c5_lhs22 = NULL;
  const mxArray *c5_rhs23 = NULL;
  const mxArray *c5_lhs23 = NULL;
  const mxArray *c5_rhs24 = NULL;
  const mxArray *c5_lhs24 = NULL;
  const mxArray *c5_rhs25 = NULL;
  const mxArray *c5_lhs25 = NULL;
  const mxArray *c5_rhs26 = NULL;
  const mxArray *c5_lhs26 = NULL;
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("rand"), "name", "name", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/rand.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383877290U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c5_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/rand.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_is_rand_extrinsic"),
                  "name", "name", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_is_rand_extrinsic.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1368183032U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c5_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/rand.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand"), "name", "name", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c5_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_str2id"), "name",
                  "name", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_str2id.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1313347822U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c5_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_str2id.m"),
                  "context", "context", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c5_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_mcg16807_stateful"),
                  "name", "name", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807_stateful.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1366162244U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c5_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807_stateful.m"),
                  "context", "context", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_mcg16807"), "name",
                  "name", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c5_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807_stateful.m"),
                  "context", "context", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_mcg16807"), "name",
                  "name", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("uint32"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c5_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_shr3cong_stateful"),
                  "name", "name", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong_stateful.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1366162244U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c5_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong_stateful.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_shr3cong"), "name",
                  "name", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c5_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong_stateful.m"),
                  "context", "context", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_shr3cong"), "name",
                  "name", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("uint32"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c5_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_mt19937ar_stateful"),
                  "name", "name", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar_stateful.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1366162244U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c5_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar_stateful.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_mt19937ar"), "name",
                  "name", 12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1406813148U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c5_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar_stateful.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_rand_mt19937ar"), "name",
                  "name", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("uint32"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1406813148U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c5_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!genrandu"),
                  "context", "context", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eps"), "name", "name", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c5_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_eps"), "name", "name", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c5_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c5_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!is_valid_state"),
                  "context", "context", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isequal"), "name", "name", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286818758U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c5_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 18);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286818786U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c5_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m!isequal_scalar"),
                  "context", "context", 19);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isnan"), "name", "name", 19);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363713858U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c5_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 20);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c5_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!is_valid_state"),
                  "context", "context", 21);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 21);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c5_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!is_valid_state"),
                  "context", "context", 22);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 22);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c5_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 23);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c5_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!genrandu"),
                  "context", "context", 24);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_error"), "name", "name",
                  24);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343830358U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c5_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 25);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 25);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c5_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 26);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 26);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c5_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs26), "lhs", "lhs",
                  26);
  sf_mex_destroy(&c5_rhs0);
  sf_mex_destroy(&c5_lhs0);
  sf_mex_destroy(&c5_rhs1);
  sf_mex_destroy(&c5_lhs1);
  sf_mex_destroy(&c5_rhs2);
  sf_mex_destroy(&c5_lhs2);
  sf_mex_destroy(&c5_rhs3);
  sf_mex_destroy(&c5_lhs3);
  sf_mex_destroy(&c5_rhs4);
  sf_mex_destroy(&c5_lhs4);
  sf_mex_destroy(&c5_rhs5);
  sf_mex_destroy(&c5_lhs5);
  sf_mex_destroy(&c5_rhs6);
  sf_mex_destroy(&c5_lhs6);
  sf_mex_destroy(&c5_rhs7);
  sf_mex_destroy(&c5_lhs7);
  sf_mex_destroy(&c5_rhs8);
  sf_mex_destroy(&c5_lhs8);
  sf_mex_destroy(&c5_rhs9);
  sf_mex_destroy(&c5_lhs9);
  sf_mex_destroy(&c5_rhs10);
  sf_mex_destroy(&c5_lhs10);
  sf_mex_destroy(&c5_rhs11);
  sf_mex_destroy(&c5_lhs11);
  sf_mex_destroy(&c5_rhs12);
  sf_mex_destroy(&c5_lhs12);
  sf_mex_destroy(&c5_rhs13);
  sf_mex_destroy(&c5_lhs13);
  sf_mex_destroy(&c5_rhs14);
  sf_mex_destroy(&c5_lhs14);
  sf_mex_destroy(&c5_rhs15);
  sf_mex_destroy(&c5_lhs15);
  sf_mex_destroy(&c5_rhs16);
  sf_mex_destroy(&c5_lhs16);
  sf_mex_destroy(&c5_rhs17);
  sf_mex_destroy(&c5_lhs17);
  sf_mex_destroy(&c5_rhs18);
  sf_mex_destroy(&c5_lhs18);
  sf_mex_destroy(&c5_rhs19);
  sf_mex_destroy(&c5_lhs19);
  sf_mex_destroy(&c5_rhs20);
  sf_mex_destroy(&c5_lhs20);
  sf_mex_destroy(&c5_rhs21);
  sf_mex_destroy(&c5_lhs21);
  sf_mex_destroy(&c5_rhs22);
  sf_mex_destroy(&c5_lhs22);
  sf_mex_destroy(&c5_rhs23);
  sf_mex_destroy(&c5_lhs23);
  sf_mex_destroy(&c5_rhs24);
  sf_mex_destroy(&c5_lhs24);
  sf_mex_destroy(&c5_rhs25);
  sf_mex_destroy(&c5_lhs25);
  sf_mex_destroy(&c5_rhs26);
  sf_mex_destroy(&c5_lhs26);
}

static const mxArray *c5_emlrt_marshallOut(const char * c5_u)
{
  const mxArray *c5_y = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c5_u)), false);
  return c5_y;
}

static const mxArray *c5_b_emlrt_marshallOut(const uint32_T c5_u)
{
  const mxArray *c5_y = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 7, 0U, 0U, 0U, 0), false);
  return c5_y;
}

static real_T c5_rand(SFc5_simulation2InstanceStruct *chartInstance)
{
  real_T c5_r;
  uint32_T c5_V4;
  uint32_T c5_V5NORMAL;
  uint32_T c5_TWISTER;
  uint32_T c5_d_state;
  real_T c5_d1;
  int32_T c5_i17;
  uint32_T c5_icng;
  uint32_T c5_jsr;
  uint32_T c5_u0;
  uint32_T c5_u1;
  uint32_T c5_ui;
  uint32_T c5_b_ui;
  real_T c5_b_r;
  real_T c5_d2;
  uint32_T c5_uv2[625];
  int32_T c5_i18;
  real_T c5_d3;
  c5_V4 = c5_eml_rand_str2id(chartInstance);
  c5_V5NORMAL = c5_b_eml_rand_str2id(chartInstance);
  c5_TWISTER = c5_c_eml_rand_str2id(chartInstance);
  if (!chartInstance->c5_method_not_empty) {
    chartInstance->c5_method = c5_TWISTER;
    chartInstance->c5_method_not_empty = true;
  }

  if (chartInstance->c5_method == c5_V4) {
    if (!chartInstance->c5_state_not_empty) {
      chartInstance->c5_state = 1144108930U;
      chartInstance->c5_state_not_empty = true;
    }

    c5_d_state = chartInstance->c5_state;
    c5_d1 = c5_b_eml_rand_mcg16807(chartInstance, &c5_d_state);
    chartInstance->c5_state = c5_d_state;
    c5_r = c5_d1;
  } else if (chartInstance->c5_method == c5_V5NORMAL) {
    if (!chartInstance->c5_b_state_not_empty) {
      for (c5_i17 = 0; c5_i17 < 2; c5_i17++) {
        chartInstance->c5_b_state[c5_i17] = 362436069U + 158852560U * (uint32_T)
          c5_i17;
      }

      chartInstance->c5_b_state_not_empty = true;
    }

    c5_icng = chartInstance->c5_b_state[0];
    c5_jsr = chartInstance->c5_b_state[1];
    c5_u0 = c5_jsr;
    c5_u1 = c5_icng;
    c5_u1 = 69069U * c5_u1 + 1234567U;
    c5_u0 ^= c5_u0 << 13;
    c5_u0 ^= c5_u0 >> 17;
    c5_u0 ^= c5_u0 << 5;
    c5_ui = c5_u1 + c5_u0;
    chartInstance->c5_b_state[0] = c5_u1;
    chartInstance->c5_b_state[1] = c5_u0;
    c5_b_ui = c5_ui;
    c5_b_r = (real_T)c5_b_ui * 2.328306436538696E-10;
    c5_d2 = c5_b_r;
    c5_r = c5_d2;
  } else {
    if (!chartInstance->c5_c_state_not_empty) {
      c5_eml_rand_mt19937ar(chartInstance, c5_uv2);
      for (c5_i18 = 0; c5_i18 < 625; c5_i18++) {
        chartInstance->c5_c_state[c5_i18] = c5_uv2[c5_i18];
      }

      chartInstance->c5_c_state_not_empty = true;
    }

    c5_d3 = c5_c_eml_rand_mt19937ar(chartInstance, chartInstance->c5_c_state);
    c5_r = c5_d3;
  }

  return c5_r;
}

static uint32_T c5_eml_rand_str2id(SFc5_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
  return 4U;
}

static uint32_T c5_b_eml_rand_str2id(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
  return 5U;
}

static uint32_T c5_c_eml_rand_str2id(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
  return 7U;
}

static void c5_eml_rand_mcg16807(SFc5_simulation2InstanceStruct *chartInstance,
  uint32_T c5_d_state, uint32_T *c5_e_state, real_T *c5_r)
{
  *c5_e_state = c5_d_state;
  *c5_r = c5_b_eml_rand_mcg16807(chartInstance, c5_e_state);
}

static void c5_eml_rand_mt19937ar(SFc5_simulation2InstanceStruct *chartInstance,
  uint32_T c5_d_state[625])
{
  int32_T c5_i19;
  for (c5_i19 = 0; c5_i19 < 625; c5_i19++) {
    c5_d_state[c5_i19] = 0U;
  }

  c5_b_twister_state_vector(chartInstance, c5_d_state, 5489.0);
}

static void c5_twister_state_vector(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T c5_mt[625], real_T c5_seed, uint32_T c5_b_mt[625])
{
  int32_T c5_i20;
  for (c5_i20 = 0; c5_i20 < 625; c5_i20++) {
    c5_b_mt[c5_i20] = c5_mt[c5_i20];
  }

  c5_b_twister_state_vector(chartInstance, c5_b_mt, c5_seed);
}

static void c5_b_eml_rand_mt19937ar(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T c5_d_state[625], uint32_T c5_e_state[625], real_T
  *c5_r)
{
  int32_T c5_i21;
  for (c5_i21 = 0; c5_i21 < 625; c5_i21++) {
    c5_e_state[c5_i21] = c5_d_state[c5_i21];
  }

  *c5_r = c5_c_eml_rand_mt19937ar(chartInstance, c5_e_state);
}

static void c5_eml_error(SFc5_simulation2InstanceStruct *chartInstance)
{
  int32_T c5_i22;
  static char_T c5_cv0[37] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'r', 'a', 'n', 'd', '_', 'i', 'n', 'v', 'a', 'l', 'i', 'd',
    'T', 'w', 'i', 's', 't', 'e', 'r', 'S', 't', 'a', 't', 'e' };

  char_T c5_u[37];
  const mxArray *c5_y = NULL;
  (void)chartInstance;
  for (c5_i22 = 0; c5_i22 < 37; c5_i22++) {
    c5_u[c5_i22] = c5_cv0[c5_i22];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 37), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c5_y));
}

static void c5_b_rand(SFc5_simulation2InstanceStruct *chartInstance, real_T
                      c5_r[3])
{
  uint32_T c5_V4;
  uint32_T c5_V5NORMAL;
  uint32_T c5_TWISTER;
  int32_T c5_k;
  real_T c5_b_k;
  uint32_T c5_d_state;
  real_T c5_d4;
  int32_T c5_i23;
  int32_T c5_c_k;
  real_T c5_d_k;
  uint32_T c5_icng;
  uint32_T c5_jsr;
  uint32_T c5_u2;
  uint32_T c5_u3;
  uint32_T c5_ui;
  uint32_T c5_b_ui;
  real_T c5_b_r;
  real_T c5_d5;
  uint32_T c5_uv3[625];
  int32_T c5_i24;
  int32_T c5_e_k;
  real_T c5_f_k;
  real_T c5_d6;
  c5_V4 = c5_eml_rand_str2id(chartInstance);
  c5_V5NORMAL = c5_b_eml_rand_str2id(chartInstance);
  c5_TWISTER = c5_c_eml_rand_str2id(chartInstance);
  if (!chartInstance->c5_method_not_empty) {
    chartInstance->c5_method = c5_TWISTER;
    chartInstance->c5_method_not_empty = true;
  }

  if (chartInstance->c5_method == c5_V4) {
    if (!chartInstance->c5_state_not_empty) {
      chartInstance->c5_state = 1144108930U;
      chartInstance->c5_state_not_empty = true;
    }

    for (c5_k = 0; c5_k < 3; c5_k++) {
      c5_b_k = 1.0 + (real_T)c5_k;
      c5_d_state = chartInstance->c5_state;
      c5_d4 = c5_b_eml_rand_mcg16807(chartInstance, &c5_d_state);
      chartInstance->c5_state = c5_d_state;
      c5_r[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_k), 1, 3, 1, 0) - 1] = c5_d4;
    }
  } else if (chartInstance->c5_method == c5_V5NORMAL) {
    if (!chartInstance->c5_b_state_not_empty) {
      for (c5_i23 = 0; c5_i23 < 2; c5_i23++) {
        chartInstance->c5_b_state[c5_i23] = 362436069U + 158852560U * (uint32_T)
          c5_i23;
      }

      chartInstance->c5_b_state_not_empty = true;
    }

    for (c5_c_k = 0; c5_c_k < 3; c5_c_k++) {
      c5_d_k = 1.0 + (real_T)c5_c_k;
      c5_icng = chartInstance->c5_b_state[0];
      c5_jsr = chartInstance->c5_b_state[1];
      c5_u2 = c5_jsr;
      c5_u3 = c5_icng;
      c5_u3 = 69069U * c5_u3 + 1234567U;
      c5_u2 ^= c5_u2 << 13;
      c5_u2 ^= c5_u2 >> 17;
      c5_u2 ^= c5_u2 << 5;
      c5_ui = c5_u3 + c5_u2;
      chartInstance->c5_b_state[0] = c5_u3;
      chartInstance->c5_b_state[1] = c5_u2;
      c5_b_ui = c5_ui;
      c5_b_r = (real_T)c5_b_ui * 2.328306436538696E-10;
      c5_d5 = c5_b_r;
      c5_r[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_d_k), 1, 3, 1, 0) - 1] = c5_d5;
    }
  } else {
    if (!chartInstance->c5_c_state_not_empty) {
      c5_eml_rand_mt19937ar(chartInstance, c5_uv3);
      for (c5_i24 = 0; c5_i24 < 625; c5_i24++) {
        chartInstance->c5_c_state[c5_i24] = c5_uv3[c5_i24];
      }

      chartInstance->c5_c_state_not_empty = true;
    }

    for (c5_e_k = 0; c5_e_k < 3; c5_e_k++) {
      c5_f_k = 1.0 + (real_T)c5_e_k;
      c5_d6 = c5_c_eml_rand_mt19937ar(chartInstance, chartInstance->c5_c_state);
      c5_r[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_f_k), 1, 3, 1, 0) - 1] = c5_d6;
    }
  }
}

static const mxArray *c5_c_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc5_simulation2InstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(int32_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static int32_T c5_d_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  int32_T c5_y;
  int32_T c5_i25;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_i25, 1, 6, 0U, 0, 0U, 0);
  c5_y = c5_i25;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_b_sfEvent;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  int32_T c5_y;
  SFc5_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc5_simulation2InstanceStruct *)chartInstanceVoid;
  c5_b_sfEvent = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_sfEvent),
    &c5_thisId);
  sf_mex_destroy(&c5_b_sfEvent);
  *(int32_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

static uint32_T c5_e_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_b_method, const char_T *c5_identifier)
{
  uint32_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_method),
    &c5_thisId);
  sf_mex_destroy(&c5_b_method);
  return c5_y;
}

static uint32_T c5_f_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  uint32_T c5_y;
  uint32_T c5_u4;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_method_not_empty = false;
  } else {
    chartInstance->c5_method_not_empty = true;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_u4, 1, 7, 0U, 0, 0U, 0);
    c5_y = c5_u4;
  }

  sf_mex_destroy(&c5_u);
  return c5_y;
}

static uint32_T c5_g_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_d_state, const char_T *c5_identifier)
{
  uint32_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_d_state), &c5_thisId);
  sf_mex_destroy(&c5_d_state);
  return c5_y;
}

static uint32_T c5_h_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  uint32_T c5_y;
  uint32_T c5_u5;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_state_not_empty = false;
  } else {
    chartInstance->c5_state_not_empty = true;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_u5, 1, 7, 0U, 0, 0U, 0);
    c5_y = c5_u5;
  }

  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_i_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_d_state, const char_T *c5_identifier, uint32_T c5_y[625])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_d_state), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_d_state);
}

static void c5_j_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, uint32_T c5_y[625])
{
  uint32_T c5_uv4[625];
  int32_T c5_i26;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_c_state_not_empty = false;
  } else {
    chartInstance->c5_c_state_not_empty = true;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_uv4, 1, 7, 0U, 1, 0U, 1, 625);
    for (c5_i26 = 0; c5_i26 < 625; c5_i26++) {
      c5_y[c5_i26] = c5_uv4[c5_i26];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_k_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_d_state, const char_T *c5_identifier, uint32_T c5_y[2])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_d_state), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_d_state);
}

static void c5_l_emlrt_marshallIn(SFc5_simulation2InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, uint32_T c5_y[2])
{
  uint32_T c5_uv5[2];
  int32_T c5_i27;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_b_state_not_empty = false;
  } else {
    chartInstance->c5_b_state_not_empty = true;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_uv5, 1, 7, 0U, 1, 0U, 1, 2);
    for (c5_i27 = 0; c5_i27 < 2; c5_i27++) {
      c5_y[c5_i27] = c5_uv5[c5_i27];
    }
  }

  sf_mex_destroy(&c5_u);
}

static uint8_T c5_m_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_b_is_active_c5_simulation2, const char_T
  *c5_identifier)
{
  uint8_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_n_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c5_b_is_active_c5_simulation2), &c5_thisId);
  sf_mex_destroy(&c5_b_is_active_c5_simulation2);
  return c5_y;
}

static uint8_T c5_n_emlrt_marshallIn(SFc5_simulation2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  uint8_T c5_y;
  uint8_T c5_u6;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_u6, 1, 3, 0U, 0, 0U, 0);
  c5_y = c5_u6;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static real_T c5_b_eml_rand_mcg16807(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T *c5_d_state)
{
  uint32_T c5_s;
  uint32_T c5_u7;
  uint32_T c5_hi;
  uint32_T c5_lo;
  uint32_T c5_test1;
  uint32_T c5_test2;
  uint32_T c5_e_state;
  real_T c5_b_r;
  (void)chartInstance;
  c5_s = *c5_d_state;
  c5_u7 = 127773U;
  if (c5_u7 == 0U) {
    c5_hi = MAX_uint32_T;
  } else {
    c5_hi = c5_s / c5_u7;
  }

  c5_lo = c5_s - c5_hi * 127773U;
  c5_test1 = 16807U * c5_lo;
  c5_test2 = 2836U * c5_hi;
  if (c5_test1 < c5_test2) {
    c5_e_state = (c5_test1 - c5_test2) + 2147483647U;
  } else {
    c5_e_state = c5_test1 - c5_test2;
  }

  c5_b_r = (real_T)c5_e_state * 4.6566128752457969E-10;
  *c5_d_state = c5_e_state;
  return c5_b_r;
}

static void c5_b_twister_state_vector(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T c5_mt[625], real_T c5_seed)
{
  real_T c5_d7;
  uint32_T c5_u8;
  uint32_T c5_r;
  int32_T c5_mti;
  real_T c5_b_mti;
  real_T c5_d8;
  uint32_T c5_u9;
  (void)chartInstance;
  c5_d7 = c5_seed;
  if (c5_d7 < 4.294967296E+9) {
    if (c5_d7 >= 0.0) {
      c5_u8 = (uint32_T)c5_d7;
    } else {
      c5_u8 = 0U;
    }
  } else if (c5_d7 >= 4.294967296E+9) {
    c5_u8 = MAX_uint32_T;
  } else {
    c5_u8 = 0U;
  }

  c5_r = c5_u8;
  c5_mt[0] = c5_r;
  for (c5_mti = 0; c5_mti < 623; c5_mti++) {
    c5_b_mti = 2.0 + (real_T)c5_mti;
    c5_d8 = muDoubleScalarRound(c5_b_mti - 1.0);
    if (c5_d8 < 4.294967296E+9) {
      if (c5_d8 >= 0.0) {
        c5_u9 = (uint32_T)c5_d8;
      } else {
        c5_u9 = 0U;
      }
    } else if (c5_d8 >= 4.294967296E+9) {
      c5_u9 = MAX_uint32_T;
    } else {
      c5_u9 = 0U;
    }

    c5_r = (c5_r ^ c5_r >> 30U) * 1812433253U + c5_u9;
    c5_mt[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c5_b_mti), 1, 625, 1, 0) - 1] = c5_r;
  }

  c5_mt[624] = 624U;
}

static real_T c5_c_eml_rand_mt19937ar(SFc5_simulation2InstanceStruct
  *chartInstance, uint32_T c5_d_state[625])
{
  int32_T c5_i28;
  uint32_T c5_u[2];
  int32_T c5_j;
  real_T c5_b_j;
  uint32_T c5_mti;
  int32_T c5_kk;
  real_T c5_b_kk;
  uint32_T c5_y;
  uint32_T c5_b_y;
  uint32_T c5_c_y;
  int32_T c5_c_kk;
  uint32_T c5_d_y;
  uint32_T c5_e_y;
  uint32_T c5_f_y;
  uint32_T c5_g_y;
  real_T c5_b_r;
  boolean_T c5_b0;
  boolean_T c5_isvalid;
  int32_T c5_k;
  int32_T c5_a;
  int32_T c5_b_a;
  boolean_T guard1 = false;
  int32_T exitg1;
  boolean_T exitg2;

  /* ========================= COPYRIGHT NOTICE ============================ */
  /*  This is a uniform (0,1) pseudorandom number generator based on:        */
  /*                                                                         */
  /*  A C-program for MT19937, with initialization improved 2002/1/26.       */
  /*  Coded by Takuji Nishimura and Makoto Matsumoto.                        */
  /*                                                                         */
  /*  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,      */
  /*  All rights reserved.                                                   */
  /*                                                                         */
  /*  Redistribution and use in source and binary forms, with or without     */
  /*  modification, are permitted provided that the following conditions     */
  /*  are met:                                                               */
  /*                                                                         */
  /*    1. Redistributions of source code must retain the above copyright    */
  /*       notice, this list of conditions and the following disclaimer.     */
  /*                                                                         */
  /*    2. Redistributions in binary form must reproduce the above copyright */
  /*       notice, this list of conditions and the following disclaimer      */
  /*       in the documentation and/or other materials provided with the     */
  /*       distribution.                                                     */
  /*                                                                         */
  /*    3. The names of its contributors may not be used to endorse or       */
  /*       promote products derived from this software without specific      */
  /*       prior written permission.                                         */
  /*                                                                         */
  /*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS    */
  /*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT      */
  /*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  */
  /*  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT  */
  /*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,  */
  /*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT       */
  /*  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,  */
  /*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY  */
  /*  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT    */
  /*  (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE */
  /*  OF THIS  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */
  /*                                                                         */
  /* =============================   END   ================================= */
  do {
    exitg1 = 0;
    for (c5_i28 = 0; c5_i28 < 2; c5_i28++) {
      c5_u[c5_i28] = 0U;
    }

    for (c5_j = 0; c5_j < 2; c5_j++) {
      c5_b_j = 1.0 + (real_T)c5_j;
      c5_mti = c5_d_state[624] + 1U;
      if ((real_T)c5_mti >= 625.0) {
        for (c5_kk = 0; c5_kk < 227; c5_kk++) {
          c5_b_kk = 1.0 + (real_T)c5_kk;
          c5_y = (c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                   _SFD_INTEGER_CHECK("", c5_b_kk), 1, 625, 1, 0) - 1] &
                  2147483648U) | (c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", c5_b_kk + 1.0), 1, 625, 1, 0) - 1] &
            2147483647U);
          c5_b_y = c5_y;
          c5_c_y = c5_b_y;
          if ((real_T)(c5_c_y & 1U) == 0.0) {
            c5_c_y >>= 1U;
          } else {
            c5_c_y = c5_c_y >> 1U ^ 2567483615U;
          }

          c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c5_b_kk), 1, 625, 1, 0) - 1] =
            c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_kk + 397.0), 1, 625, 1, 0) - 1] ^ c5_c_y;
        }

        for (c5_c_kk = 0; c5_c_kk < 396; c5_c_kk++) {
          c5_b_kk = 228.0 + (real_T)c5_c_kk;
          c5_y = (c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                   _SFD_INTEGER_CHECK("", c5_b_kk), 1, 625, 1, 0) - 1] &
                  2147483648U) | (c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", c5_b_kk + 1.0), 1, 625, 1, 0) - 1] &
            2147483647U);
          c5_d_y = c5_y;
          c5_e_y = c5_d_y;
          if ((real_T)(c5_e_y & 1U) == 0.0) {
            c5_e_y >>= 1U;
          } else {
            c5_e_y = c5_e_y >> 1U ^ 2567483615U;
          }

          c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c5_b_kk), 1, 625, 1, 0) - 1] =
            c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (c5_b_kk + 1.0) - 228.0), 1, 625, 1, 0) - 1] ^
            c5_e_y;
        }

        c5_y = (c5_d_state[623] & 2147483648U) | (c5_d_state[0] & 2147483647U);
        c5_f_y = c5_y;
        c5_g_y = c5_f_y;
        if ((real_T)(c5_g_y & 1U) == 0.0) {
          c5_g_y >>= 1U;
        } else {
          c5_g_y = c5_g_y >> 1U ^ 2567483615U;
        }

        c5_d_state[623] = c5_d_state[396] ^ c5_g_y;
        c5_mti = 1U;
      }

      c5_y = c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)(uint32_T)
        _SFD_INTEGER_CHECK("", (real_T)c5_mti), 1, 625, 1, 0) - 1];
      c5_d_state[624] = c5_mti;
      c5_y ^= c5_y >> 11U;
      c5_y ^= c5_y << 7U & 2636928640U;
      c5_y ^= c5_y << 15U & 4022730752U;
      c5_y ^= c5_y >> 18U;
      c5_u[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_j), 1, 2, 1, 0) - 1] = c5_y;
    }

    c5_u[0] >>= 5U;
    c5_u[1] >>= 6U;
    c5_b_r = 1.1102230246251565E-16 * ((real_T)c5_u[0] * 6.7108864E+7 + (real_T)
      c5_u[1]);
    if (c5_b_r == 0.0) {
      guard1 = false;
      if ((real_T)c5_d_state[624] >= 1.0) {
        if ((real_T)c5_d_state[624] < 625.0) {
          c5_b0 = true;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1 == true) {
        c5_b0 = false;
      }

      c5_isvalid = c5_b0;
      if (c5_isvalid) {
        c5_isvalid = false;
        c5_k = 1;
        exitg2 = false;
        while ((exitg2 == false) && (c5_k < 625)) {
          if ((real_T)c5_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
               _SFD_INTEGER_CHECK("", (real_T)c5_k), 1, 625, 1, 0) - 1] == 0.0)
          {
            c5_a = c5_k;
            c5_b_a = c5_a + 1;
            c5_k = c5_b_a;
          } else {
            c5_isvalid = true;
            exitg2 = true;
          }
        }
      }

      if (!c5_isvalid) {
        c5_eml_error(chartInstance);
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return c5_b_r;
}

static void init_dsm_address_info(SFc5_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc5_simulation2InstanceStruct
  *chartInstance)
{
  chartInstance->c5_d = (real_T (*)[3])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c5_dist = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
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

void sf_c5_simulation2_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(83396376U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2944532718U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3001722619U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2954002397U);
}

mxArray* sf_c5_simulation2_get_post_codegen_info(void);
mxArray *sf_c5_simulation2_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("HzAtdXdcyFe7hUNkuJjveE");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
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
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
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
    mxArray* mxPostCodegenInfo = sf_c5_simulation2_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c5_simulation2_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c5_simulation2_jit_fallback_info(void)
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

mxArray *sf_c5_simulation2_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c5_simulation2_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c5_simulation2(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x6'type','srcId','name','auxInfo'{{M[1],M[5],T\"d\",},{M[4],M[0],T\"method\",S'l','i','p'{{M1x2[512 518],M[1],T\"/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/randfun/eml_rand.m\"}}},{M[4],M[0],T\"state\",S'l','i','p'{{M1x2[165 170],M[1],T\"/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807_stateful.m\"}}},{M[4],M[0],T\"state\",S'l','i','p'{{M1x2[166 171],M[1],T\"/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar_stateful.m\"}}},{M[4],M[0],T\"state\",S'l','i','p'{{M1x2[165 170],M[1],T\"/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong_stateful.m\"}}},{M[8],M[0],T\"is_active_c5_simulation2\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 6, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c5_simulation2_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc5_simulation2InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc5_simulation2InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulation2MachineNumber_,
           5,
           1,
           1,
           0,
           2,
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
          _SFD_SET_DATA_PROPS(0,2,0,1,"d");
          _SFD_SET_DATA_PROPS(1,1,1,0,"dist");
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
        _SFD_CV_INIT_EML(0,1,1,1,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,151);
        _SFD_CV_INIT_EML_IF(0,1,0,40,48,125,150);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)
            c5_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c5_d);
        _SFD_SET_DATA_VALUE_PTR(1U, chartInstance->c5_dist);
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
  return "YbHS5tcIBhBIvr0F90Fpt";
}

static void sf_opaque_initialize_c5_simulation2(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc5_simulation2InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c5_simulation2((SFc5_simulation2InstanceStruct*)
    chartInstanceVar);
  initialize_c5_simulation2((SFc5_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c5_simulation2(void *chartInstanceVar)
{
  enable_c5_simulation2((SFc5_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c5_simulation2(void *chartInstanceVar)
{
  disable_c5_simulation2((SFc5_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c5_simulation2(void *chartInstanceVar)
{
  sf_gateway_c5_simulation2((SFc5_simulation2InstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c5_simulation2(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c5_simulation2((SFc5_simulation2InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c5_simulation2(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c5_simulation2((SFc5_simulation2InstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c5_simulation2(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc5_simulation2InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation2_optimization_info();
    }

    finalize_c5_simulation2((SFc5_simulation2InstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc5_simulation2((SFc5_simulation2InstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c5_simulation2(SimStruct *S)
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
    initialize_params_c5_simulation2((SFc5_simulation2InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c5_simulation2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation2_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,5);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,5,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,5,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,5);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,5,1);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,5,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 1; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,5);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(637601677U));
  ssSetChecksum1(S,(1886813661U));
  ssSetChecksum2(S,(3618364336U));
  ssSetChecksum3(S,(3611621209U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c5_simulation2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c5_simulation2(SimStruct *S)
{
  SFc5_simulation2InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc5_simulation2InstanceStruct *)utMalloc(sizeof
    (SFc5_simulation2InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc5_simulation2InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c5_simulation2;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c5_simulation2;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c5_simulation2;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c5_simulation2;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c5_simulation2;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c5_simulation2;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c5_simulation2;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c5_simulation2;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c5_simulation2;
  chartInstance->chartInfo.mdlStart = mdlStart_c5_simulation2;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c5_simulation2;
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

void c5_simulation2_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c5_simulation2(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c5_simulation2(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c5_simulation2(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c5_simulation2_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
