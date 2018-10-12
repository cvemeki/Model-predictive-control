/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation1_sfun.h"
#include "c3_simulation1.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simulation1_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c3_debug_family_names[13] = { "roll", "pitch", "yaw", "R",
  "g", "d", "quad", "nargin", "nargout", "uTot", "dist", "theta", "vDot" };

/* Function Declarations */
static void initialize_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance);
static void initialize_params_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance);
static void enable_c3_simulation1(SFc3_simulation1InstanceStruct *chartInstance);
static void disable_c3_simulation1(SFc3_simulation1InstanceStruct *chartInstance);
static void c3_update_debugger_state_c3_simulation1
  (SFc3_simulation1InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c3_simulation1
  (SFc3_simulation1InstanceStruct *chartInstance);
static void set_sim_state_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_st);
static void finalize_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance);
static void sf_gateway_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance);
static void mdl_start_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance);
static void c3_chartstep_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance);
static void initSimStructsc3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber);
static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData);
static void c3_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_b_vDot, const char_T *c3_identifier, real_T c3_y[3]);
static void c3_b_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[3]);
static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static real_T c3_c_emlrt_marshallIn(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_d_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[9]);
static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_info_helper(const mxArray **c3_info);
static const mxArray *c3_emlrt_marshallOut(const char * c3_u);
static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u);
static void c3_eml_scalar_eg(SFc3_simulation1InstanceStruct *chartInstance);
static void c3_threshold(SFc3_simulation1InstanceStruct *chartInstance);
static void c3_b_eml_scalar_eg(SFc3_simulation1InstanceStruct *chartInstance);
static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static int32_T c3_e_emlrt_marshallIn(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_f_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  c3_sdz1UxT0ktEtEZMZuZvzEsH *c3_y);
static void c3_g_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[8]);
static void c3_h_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[4]);
static void c3_i_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[12]);
static void c3_j_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[16]);
static void c3_k_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[2]);
static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static uint8_T c3_l_emlrt_marshallIn(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_b_is_active_c3_simulation1, const char_T
  *c3_identifier);
static uint8_T c3_m_emlrt_marshallIn(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void init_dsm_address_info(SFc3_simulation1InstanceStruct *chartInstance);
static void init_simulink_io_address(SFc3_simulation1InstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance)
{
  chartInstance->c3_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c3_is_active_c3_simulation1 = 0U;
}

static void initialize_params_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance)
{
  const mxArray *c3_m0 = NULL;
  const mxArray *c3_mxField;
  c3_sdz1UxT0ktEtEZMZuZvzEsH c3_r0;
  c3_m0 = sf_mex_get_sfun_param(chartInstance->S, 0, 1);
  c3_mxField = sf_mex_getfield(c3_m0, "thrustLimits", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.thrustLimits, 1, 0,
                      0U, 1, 0U, 2, 2, 4);
  c3_mxField = sf_mex_getfield(c3_m0, "thrustRateLimits", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.thrustRateLimits, 1,
                      0, 0U, 1, 0U, 2, 1, 4);
  c3_mxField = sf_mex_getfield(c3_m0, "thrustDir", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.thrustDir, 1, 0, 0U,
                      1, 0U, 2, 3, 4);
  c3_mxField = sf_mex_getfield(c3_m0, "L", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.L, 1, 0, 0U, 1, 0U,
                      2, 3, 4);
  c3_mxField = sf_mex_getfield(c3_m0, "nL", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.nL, 1, 0, 0U, 0, 0U,
                      0);
  c3_mxField = sf_mex_getfield(c3_m0, "rad", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.rad, 1, 0, 0U, 0,
                      0U, 0);
  c3_mxField = sf_mex_getfield(c3_m0, "bladeRad", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.bladeRad, 1, 0, 0U,
                      0, 0U, 0);
  c3_mxField = sf_mex_getfield(c3_m0, "mass", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.mass, 1, 0, 0U, 0,
                      0U, 0);
  c3_mxField = sf_mex_getfield(c3_m0, "I", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.I, 1, 0, 0U, 1, 0U,
                      2, 3, 3);
  c3_mxField = sf_mex_getfield(c3_m0, "km", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.km, 1, 0, 0U, 0, 0U,
                      0);
  c3_mxField = sf_mex_getfield(c3_m0, "kf", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.kf, 1, 0, 0U, 0, 0U,
                      0);
  c3_mxField = sf_mex_getfield(c3_m0, "K", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.K, 1, 0, 0U, 1, 0U,
                      2, 4, 4);
  c3_mxField = sf_mex_getfield(c3_m0, "velocityGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.velocityGain, 1, 0,
                      0U, 0, 0U, 0);
  c3_mxField = sf_mex_getfield(c3_m0, "vGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.vGain, 1, 0, 0U, 1,
                      0U, 1, 2);
  c3_mxField = sf_mex_getfield(c3_m0, "zGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.zGain, 1, 0, 0U, 0,
                      0U, 0);
  c3_mxField = sf_mex_getfield(c3_m0, "thetaGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.thetaGain, 1, 0, 0U,
                      1, 0U, 1, 3);
  c3_mxField = sf_mex_getfield(c3_m0, "omegaGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), c3_r0.omegaGain, 1, 0, 0U,
                      1, 0U, 1, 3);
  c3_mxField = sf_mex_getfield(c3_m0, "Ts2", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c3_mxField), &c3_r0.Ts2, 1, 0, 0U, 0,
                      0U, 0);
  sf_mex_destroy(&c3_m0);
  chartInstance->c3_quad = c3_r0;
}

static void enable_c3_simulation1(SFc3_simulation1InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c3_simulation1(SFc3_simulation1InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c3_update_debugger_state_c3_simulation1
  (SFc3_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c3_simulation1
  (SFc3_simulation1InstanceStruct *chartInstance)
{
  const mxArray *c3_st;
  const mxArray *c3_y = NULL;
  int32_T c3_i0;
  real_T c3_u[3];
  const mxArray *c3_b_y = NULL;
  uint8_T c3_hoistedGlobal;
  uint8_T c3_b_u;
  const mxArray *c3_c_y = NULL;
  c3_st = NULL;
  c3_st = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createcellmatrix(2, 1), false);
  for (c3_i0 = 0; c3_i0 < 3; c3_i0++) {
    c3_u[c3_i0] = (*chartInstance->c3_vDot)[c3_i0];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c3_y, 0, c3_b_y);
  c3_hoistedGlobal = chartInstance->c3_is_active_c3_simulation1;
  c3_b_u = c3_hoistedGlobal;
  c3_c_y = NULL;
  sf_mex_assign(&c3_c_y, sf_mex_create("y", &c3_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c3_y, 1, c3_c_y);
  sf_mex_assign(&c3_st, c3_y, false);
  return c3_st;
}

static void set_sim_state_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_st)
{
  const mxArray *c3_u;
  real_T c3_dv0[3];
  int32_T c3_i1;
  chartInstance->c3_doneDoubleBufferReInit = true;
  c3_u = sf_mex_dup(c3_st);
  c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 0)), "vDot",
                      c3_dv0);
  for (c3_i1 = 0; c3_i1 < 3; c3_i1++) {
    (*chartInstance->c3_vDot)[c3_i1] = c3_dv0[c3_i1];
  }

  chartInstance->c3_is_active_c3_simulation1 = c3_l_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 1)),
     "is_active_c3_simulation1");
  sf_mex_destroy(&c3_u);
  c3_update_debugger_state_c3_simulation1(chartInstance);
  sf_mex_destroy(&c3_st);
}

static void finalize_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance)
{
  int32_T c3_i2;
  int32_T c3_i3;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c3_uTot, 0U);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c3_dist, 1U);
  chartInstance->c3_sfEvent = CALL_EVENT;
  c3_chartstep_c3_simulation1(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulation1MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c3_i2 = 0; c3_i2 < 3; c3_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c3_vDot)[c3_i2], 2U);
  }

  for (c3_i3 = 0; c3_i3 < 3; c3_i3++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c3_theta)[c3_i3], 3U);
  }
}

static void mdl_start_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c3_chartstep_c3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance)
{
  real_T c3_hoistedGlobal;
  real_T c3_b_hoistedGlobal;
  real_T c3_b_uTot;
  real_T c3_b_dist;
  int32_T c3_i4;
  real_T c3_b_theta[3];
  uint32_T c3_debug_family_var_map[13];
  real_T c3_roll;
  real_T c3_pitch;
  real_T c3_yaw;
  real_T c3_R[9];
  real_T c3_g;
  real_T c3_d[3];
  c3_sdz1UxT0ktEtEZMZuZvzEsH c3_b_quad;
  real_T c3_nargin = 4.0;
  real_T c3_nargout = 1.0;
  real_T c3_b_vDot[3];
  static c3_sdz1UxT0ktEtEZMZuZvzEsH c3_r1 = { { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
      0.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 1.0 }, { 0.0, 0.2, 0.0, -0.2, 0.0, 0.0, 0.0, -0.2, 0.0,
      0.2, 0.0, 0.0 }, 0.2, 0.040000000000000008, 0.080000000000000016, 8.0, {
      10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 15.0 }, 11.0, 28.0, { 28.0, 0.0,
      -5.6000000000000005, 11.0, 28.0, 5.6000000000000005, 0.0, -11.0, 28.0, 0.0,
      5.6000000000000005, 11.0, 28.0, -5.6000000000000005, 0.0, -11.0 }, 10.0, {
      0.1, 0.1 }, 1.0, { 1.0, 1.0, 1.0 }, { 0.2, 0.2, 0.2 }, 0.3 };

  real_T c3_x;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_d_x;
  real_T c3_e_x;
  real_T c3_f_x;
  real_T c3_g_x;
  real_T c3_h_x;
  int32_T c3_i5;
  int32_T c3_i6;
  static real_T c3_dv1[3] = { 1.0, 0.0, 0.0 };

  real_T c3_i_x;
  real_T c3_j_x;
  real_T c3_k_x;
  real_T c3_l_x;
  real_T c3_m_x;
  real_T c3_n_x;
  real_T c3_o_x;
  real_T c3_p_x;
  int32_T c3_i7;
  real_T c3_a[9];
  real_T c3_b[9];
  int32_T c3_i8;
  int32_T c3_i9;
  static real_T c3_dv2[3] = { 0.0, 1.0, 0.0 };

  int32_T c3_i10;
  int32_T c3_i11;
  int32_T c3_i12;
  real_T c3_C[9];
  int32_T c3_i13;
  int32_T c3_i14;
  int32_T c3_i15;
  int32_T c3_i16;
  int32_T c3_i17;
  int32_T c3_i18;
  int32_T c3_i19;
  int32_T c3_i20;
  real_T c3_q_x;
  real_T c3_r_x;
  real_T c3_s_x;
  real_T c3_t_x;
  real_T c3_u_x;
  real_T c3_v_x;
  real_T c3_w_x;
  real_T c3_x_x;
  int32_T c3_i21;
  int32_T c3_i22;
  int32_T c3_i23;
  static real_T c3_dv3[3] = { 0.0, 0.0, 1.0 };

  int32_T c3_i24;
  int32_T c3_i25;
  int32_T c3_i26;
  int32_T c3_i27;
  int32_T c3_i28;
  int32_T c3_i29;
  int32_T c3_i30;
  int32_T c3_i31;
  int32_T c3_i32;
  int32_T c3_i33;
  int32_T c3_i34;
  int32_T c3_i35;
  static real_T c3_dv4[3] = { -0.0, -0.0, -8.0 };

  int32_T c3_i36;
  real_T c3_b_a;
  int32_T c3_i37;
  int32_T c3_i38;
  int32_T c3_i39;
  real_T c3_y[3];
  int32_T c3_i40;
  int32_T c3_i41;
  static real_T c3_b_b[3] = { 0.0, 0.0, 1.0 };

  int32_T c3_i42;
  static real_T c3_b_y[3] = { -0.0, -0.0, -78.48 };

  int32_T c3_i43;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
  c3_hoistedGlobal = *chartInstance->c3_uTot;
  c3_b_hoistedGlobal = *chartInstance->c3_dist;
  c3_b_uTot = c3_hoistedGlobal;
  c3_b_dist = c3_b_hoistedGlobal;
  for (c3_i4 = 0; c3_i4 < 3; c3_i4++) {
    c3_b_theta[c3_i4] = (*chartInstance->c3_theta)[c3_i4];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 13U, 13U, c3_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_roll, 0U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_pitch, 1U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_yaw, 2U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_R, 3U, c3_d_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_g, 4U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_d, 5U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_b_quad, 6U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 7U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 8U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_b_uTot, 9U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_b_dist, 10U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_b_theta, 11U, c3_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_b_vDot, 12U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  c3_b_quad = c3_r1;
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 5);
  c3_roll = c3_b_theta[0];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 5);
  c3_pitch = c3_b_theta[1];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 5);
  c3_yaw = c3_b_theta[2];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 6);
  c3_x = c3_roll;
  c3_b_x = c3_x;
  c3_b_x = muDoubleScalarCos(c3_b_x);
  c3_c_x = c3_roll;
  c3_d_x = c3_c_x;
  c3_d_x = muDoubleScalarSin(c3_d_x);
  c3_e_x = c3_roll;
  c3_f_x = c3_e_x;
  c3_f_x = muDoubleScalarSin(c3_f_x);
  c3_g_x = c3_roll;
  c3_h_x = c3_g_x;
  c3_h_x = muDoubleScalarCos(c3_h_x);
  c3_i5 = 0;
  for (c3_i6 = 0; c3_i6 < 3; c3_i6++) {
    c3_R[c3_i5] = c3_dv1[c3_i6];
    c3_i5 += 3;
  }

  c3_R[1] = 0.0;
  c3_R[4] = c3_b_x;
  c3_R[7] = -c3_d_x;
  c3_R[2] = 0.0;
  c3_R[5] = c3_f_x;
  c3_R[8] = c3_h_x;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 7);
  c3_i_x = c3_pitch;
  c3_j_x = c3_i_x;
  c3_j_x = muDoubleScalarCos(c3_j_x);
  c3_k_x = c3_pitch;
  c3_l_x = c3_k_x;
  c3_l_x = muDoubleScalarSin(c3_l_x);
  c3_m_x = c3_pitch;
  c3_n_x = c3_m_x;
  c3_n_x = muDoubleScalarSin(c3_n_x);
  c3_o_x = c3_pitch;
  c3_p_x = c3_o_x;
  c3_p_x = muDoubleScalarCos(c3_p_x);
  for (c3_i7 = 0; c3_i7 < 9; c3_i7++) {
    c3_a[c3_i7] = c3_R[c3_i7];
  }

  c3_b[0] = c3_j_x;
  c3_b[3] = 0.0;
  c3_b[6] = c3_l_x;
  c3_i8 = 0;
  for (c3_i9 = 0; c3_i9 < 3; c3_i9++) {
    c3_b[c3_i8 + 1] = c3_dv2[c3_i9];
    c3_i8 += 3;
  }

  c3_b[2] = -c3_n_x;
  c3_b[5] = 0.0;
  c3_b[8] = c3_p_x;
  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i10 = 0; c3_i10 < 9; c3_i10++) {
    c3_R[c3_i10] = 0.0;
  }

  for (c3_i11 = 0; c3_i11 < 9; c3_i11++) {
    c3_R[c3_i11] = 0.0;
  }

  for (c3_i12 = 0; c3_i12 < 9; c3_i12++) {
    c3_C[c3_i12] = c3_R[c3_i12];
  }

  for (c3_i13 = 0; c3_i13 < 9; c3_i13++) {
    c3_R[c3_i13] = c3_C[c3_i13];
  }

  c3_threshold(chartInstance);
  for (c3_i14 = 0; c3_i14 < 9; c3_i14++) {
    c3_C[c3_i14] = c3_R[c3_i14];
  }

  for (c3_i15 = 0; c3_i15 < 9; c3_i15++) {
    c3_R[c3_i15] = c3_C[c3_i15];
  }

  for (c3_i16 = 0; c3_i16 < 3; c3_i16++) {
    c3_i17 = 0;
    for (c3_i18 = 0; c3_i18 < 3; c3_i18++) {
      c3_R[c3_i17 + c3_i16] = 0.0;
      c3_i19 = 0;
      for (c3_i20 = 0; c3_i20 < 3; c3_i20++) {
        c3_R[c3_i17 + c3_i16] += c3_a[c3_i19 + c3_i16] * c3_b[c3_i20 + c3_i17];
        c3_i19 += 3;
      }

      c3_i17 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 8);
  c3_q_x = c3_yaw;
  c3_r_x = c3_q_x;
  c3_r_x = muDoubleScalarCos(c3_r_x);
  c3_s_x = c3_yaw;
  c3_t_x = c3_s_x;
  c3_t_x = muDoubleScalarSin(c3_t_x);
  c3_u_x = c3_yaw;
  c3_v_x = c3_u_x;
  c3_v_x = muDoubleScalarSin(c3_v_x);
  c3_w_x = c3_yaw;
  c3_x_x = c3_w_x;
  c3_x_x = muDoubleScalarCos(c3_x_x);
  for (c3_i21 = 0; c3_i21 < 9; c3_i21++) {
    c3_a[c3_i21] = c3_R[c3_i21];
  }

  c3_b[0] = c3_r_x;
  c3_b[3] = -c3_t_x;
  c3_b[6] = 0.0;
  c3_b[1] = c3_v_x;
  c3_b[4] = c3_x_x;
  c3_b[7] = 0.0;
  c3_i22 = 0;
  for (c3_i23 = 0; c3_i23 < 3; c3_i23++) {
    c3_b[c3_i22 + 2] = c3_dv3[c3_i23];
    c3_i22 += 3;
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i24 = 0; c3_i24 < 9; c3_i24++) {
    c3_R[c3_i24] = 0.0;
  }

  for (c3_i25 = 0; c3_i25 < 9; c3_i25++) {
    c3_R[c3_i25] = 0.0;
  }

  for (c3_i26 = 0; c3_i26 < 9; c3_i26++) {
    c3_C[c3_i26] = c3_R[c3_i26];
  }

  for (c3_i27 = 0; c3_i27 < 9; c3_i27++) {
    c3_R[c3_i27] = c3_C[c3_i27];
  }

  c3_threshold(chartInstance);
  for (c3_i28 = 0; c3_i28 < 9; c3_i28++) {
    c3_C[c3_i28] = c3_R[c3_i28];
  }

  for (c3_i29 = 0; c3_i29 < 9; c3_i29++) {
    c3_R[c3_i29] = c3_C[c3_i29];
  }

  for (c3_i30 = 0; c3_i30 < 3; c3_i30++) {
    c3_i31 = 0;
    for (c3_i32 = 0; c3_i32 < 3; c3_i32++) {
      c3_R[c3_i31 + c3_i30] = 0.0;
      c3_i33 = 0;
      for (c3_i34 = 0; c3_i34 < 3; c3_i34++) {
        c3_R[c3_i31 + c3_i30] += c3_a[c3_i33 + c3_i30] * c3_b[c3_i34 + c3_i31];
        c3_i33 += 3;
      }

      c3_i31 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 11);
  c3_g = 9.81;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 12);
  if (CV_EML_IF(0, 1, 0, c3_b_dist != 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 13);
    for (c3_i35 = 0; c3_i35 < 3; c3_i35++) {
      c3_d[c3_i35] = c3_dv4[c3_i35];
    }
  } else {
    _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 15);
    for (c3_i36 = 0; c3_i36 < 3; c3_i36++) {
      c3_d[c3_i36] = 0.0;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 17);
  c3_b_a = c3_b_uTot;
  for (c3_i37 = 0; c3_i37 < 9; c3_i37++) {
    c3_b[c3_i37] = c3_R[c3_i37];
  }

  for (c3_i38 = 0; c3_i38 < 9; c3_i38++) {
    c3_b[c3_i38] *= c3_b_a;
  }

  c3_b_eml_scalar_eg(chartInstance);
  c3_b_eml_scalar_eg(chartInstance);
  c3_threshold(chartInstance);
  for (c3_i39 = 0; c3_i39 < 3; c3_i39++) {
    c3_y[c3_i39] = 0.0;
    c3_i40 = 0;
    for (c3_i41 = 0; c3_i41 < 3; c3_i41++) {
      c3_y[c3_i39] += c3_b[c3_i40 + c3_i39] * c3_b_b[c3_i41];
      c3_i40 += 3;
    }
  }

  for (c3_i42 = 0; c3_i42 < 3; c3_i42++) {
    c3_b_vDot[c3_i42] = (c3_b_y[c3_i42] + c3_y[c3_i42]) + c3_d[c3_i42];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, -17);
  _SFD_SYMBOL_SCOPE_POP();
  for (c3_i43 = 0; c3_i43 < 3; c3_i43++) {
    (*chartInstance->c3_vDot)[c3_i43] = c3_b_vDot[c3_i43];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
}

static void initSimStructsc3_simulation1(SFc3_simulation1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber)
{
  (void)c3_machineNumber;
  (void)c3_chartNumber;
  (void)c3_instanceNumber;
}

static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i44;
  real_T c3_b_inData[3];
  int32_T c3_i45;
  real_T c3_u[3];
  const mxArray *c3_y = NULL;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i44 = 0; c3_i44 < 3; c3_i44++) {
    c3_b_inData[c3_i44] = (*(real_T (*)[3])c3_inData)[c3_i44];
  }

  for (c3_i45 = 0; c3_i45 < 3; c3_i45++) {
    c3_u[c3_i45] = c3_b_inData[c3_i45];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_b_vDot, const char_T *c3_identifier, real_T c3_y[3])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_vDot), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_b_vDot);
}

static void c3_b_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[3])
{
  real_T c3_dv5[3];
  int32_T c3_i46;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv5, 1, 0, 0U, 1, 0U, 1, 3);
  for (c3_i46 = 0; c3_i46 < 3; c3_i46++) {
    c3_y[c3_i46] = c3_dv5[c3_i46];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_vDot;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[3];
  int32_T c3_i47;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_b_vDot = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_vDot), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_b_vDot);
  for (c3_i47 = 0; c3_i47 < 3; c3_i47++) {
    (*(real_T (*)[3])c3_outData)[c3_i47] = c3_y[c3_i47];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  real_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(real_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static real_T c3_c_emlrt_marshallIn(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  real_T c3_y;
  real_T c3_d0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_d0, 1, 0, 0U, 0, 0U, 0);
  c3_y = c3_d0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_nargout;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_nargout = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_nargout), &c3_thisId);
  sf_mex_destroy(&c3_nargout);
  *(real_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  c3_sdz1UxT0ktEtEZMZuZvzEsH c3_u;
  const mxArray *c3_y = NULL;
  int32_T c3_i48;
  real_T c3_b_u[8];
  const mxArray *c3_b_y = NULL;
  int32_T c3_i49;
  real_T c3_c_u[4];
  const mxArray *c3_c_y = NULL;
  int32_T c3_i50;
  real_T c3_d_u[12];
  const mxArray *c3_d_y = NULL;
  int32_T c3_i51;
  real_T c3_e_u[12];
  const mxArray *c3_e_y = NULL;
  real_T c3_f_u;
  const mxArray *c3_f_y = NULL;
  real_T c3_g_u;
  const mxArray *c3_g_y = NULL;
  real_T c3_h_u;
  const mxArray *c3_h_y = NULL;
  real_T c3_i_u;
  const mxArray *c3_i_y = NULL;
  int32_T c3_i52;
  real_T c3_j_u[9];
  const mxArray *c3_j_y = NULL;
  real_T c3_k_u;
  const mxArray *c3_k_y = NULL;
  real_T c3_l_u;
  const mxArray *c3_l_y = NULL;
  int32_T c3_i53;
  real_T c3_m_u[16];
  const mxArray *c3_m_y = NULL;
  real_T c3_n_u;
  const mxArray *c3_n_y = NULL;
  int32_T c3_i54;
  real_T c3_o_u[2];
  const mxArray *c3_o_y = NULL;
  real_T c3_p_u;
  const mxArray *c3_p_y = NULL;
  int32_T c3_i55;
  real_T c3_q_u[3];
  const mxArray *c3_q_y = NULL;
  int32_T c3_i56;
  real_T c3_r_u[3];
  const mxArray *c3_r_y = NULL;
  real_T c3_s_u;
  const mxArray *c3_s_y = NULL;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(c3_sdz1UxT0ktEtEZMZuZvzEsH *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  for (c3_i48 = 0; c3_i48 < 8; c3_i48++) {
    c3_b_u[c3_i48] = c3_u.thrustLimits[c3_i48];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_b_u, 0, 0U, 1U, 0U, 2, 2, 4),
                false);
  sf_mex_addfield(c3_y, c3_b_y, "thrustLimits", "thrustLimits", 0);
  for (c3_i49 = 0; c3_i49 < 4; c3_i49++) {
    c3_c_u[c3_i49] = c3_u.thrustRateLimits[c3_i49];
  }

  c3_c_y = NULL;
  sf_mex_assign(&c3_c_y, sf_mex_create("y", c3_c_u, 0, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_addfield(c3_y, c3_c_y, "thrustRateLimits", "thrustRateLimits", 0);
  for (c3_i50 = 0; c3_i50 < 12; c3_i50++) {
    c3_d_u[c3_i50] = c3_u.thrustDir[c3_i50];
  }

  c3_d_y = NULL;
  sf_mex_assign(&c3_d_y, sf_mex_create("y", c3_d_u, 0, 0U, 1U, 0U, 2, 3, 4),
                false);
  sf_mex_addfield(c3_y, c3_d_y, "thrustDir", "thrustDir", 0);
  for (c3_i51 = 0; c3_i51 < 12; c3_i51++) {
    c3_e_u[c3_i51] = c3_u.L[c3_i51];
  }

  c3_e_y = NULL;
  sf_mex_assign(&c3_e_y, sf_mex_create("y", c3_e_u, 0, 0U, 1U, 0U, 2, 3, 4),
                false);
  sf_mex_addfield(c3_y, c3_e_y, "L", "L", 0);
  c3_f_u = c3_u.nL;
  c3_f_y = NULL;
  sf_mex_assign(&c3_f_y, sf_mex_create("y", &c3_f_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_f_y, "nL", "nL", 0);
  c3_g_u = c3_u.rad;
  c3_g_y = NULL;
  sf_mex_assign(&c3_g_y, sf_mex_create("y", &c3_g_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_g_y, "rad", "rad", 0);
  c3_h_u = c3_u.bladeRad;
  c3_h_y = NULL;
  sf_mex_assign(&c3_h_y, sf_mex_create("y", &c3_h_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_h_y, "bladeRad", "bladeRad", 0);
  c3_i_u = c3_u.mass;
  c3_i_y = NULL;
  sf_mex_assign(&c3_i_y, sf_mex_create("y", &c3_i_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_i_y, "mass", "mass", 0);
  for (c3_i52 = 0; c3_i52 < 9; c3_i52++) {
    c3_j_u[c3_i52] = c3_u.I[c3_i52];
  }

  c3_j_y = NULL;
  sf_mex_assign(&c3_j_y, sf_mex_create("y", c3_j_u, 0, 0U, 1U, 0U, 2, 3, 3),
                false);
  sf_mex_addfield(c3_y, c3_j_y, "I", "I", 0);
  c3_k_u = c3_u.km;
  c3_k_y = NULL;
  sf_mex_assign(&c3_k_y, sf_mex_create("y", &c3_k_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_k_y, "km", "km", 0);
  c3_l_u = c3_u.kf;
  c3_l_y = NULL;
  sf_mex_assign(&c3_l_y, sf_mex_create("y", &c3_l_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_l_y, "kf", "kf", 0);
  for (c3_i53 = 0; c3_i53 < 16; c3_i53++) {
    c3_m_u[c3_i53] = c3_u.K[c3_i53];
  }

  c3_m_y = NULL;
  sf_mex_assign(&c3_m_y, sf_mex_create("y", c3_m_u, 0, 0U, 1U, 0U, 2, 4, 4),
                false);
  sf_mex_addfield(c3_y, c3_m_y, "K", "K", 0);
  c3_n_u = c3_u.velocityGain;
  c3_n_y = NULL;
  sf_mex_assign(&c3_n_y, sf_mex_create("y", &c3_n_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_n_y, "velocityGain", "velocityGain", 0);
  for (c3_i54 = 0; c3_i54 < 2; c3_i54++) {
    c3_o_u[c3_i54] = c3_u.vGain[c3_i54];
  }

  c3_o_y = NULL;
  sf_mex_assign(&c3_o_y, sf_mex_create("y", c3_o_u, 0, 0U, 1U, 0U, 1, 2), false);
  sf_mex_addfield(c3_y, c3_o_y, "vGain", "vGain", 0);
  c3_p_u = c3_u.zGain;
  c3_p_y = NULL;
  sf_mex_assign(&c3_p_y, sf_mex_create("y", &c3_p_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_p_y, "zGain", "zGain", 0);
  for (c3_i55 = 0; c3_i55 < 3; c3_i55++) {
    c3_q_u[c3_i55] = c3_u.thetaGain[c3_i55];
  }

  c3_q_y = NULL;
  sf_mex_assign(&c3_q_y, sf_mex_create("y", c3_q_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_addfield(c3_y, c3_q_y, "thetaGain", "thetaGain", 0);
  for (c3_i56 = 0; c3_i56 < 3; c3_i56++) {
    c3_r_u[c3_i56] = c3_u.omegaGain[c3_i56];
  }

  c3_r_y = NULL;
  sf_mex_assign(&c3_r_y, sf_mex_create("y", c3_r_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_addfield(c3_y, c3_r_y, "omegaGain", "omegaGain", 0);
  c3_s_u = c3_u.Ts2;
  c3_s_y = NULL;
  sf_mex_assign(&c3_s_y, sf_mex_create("y", &c3_s_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c3_y, c3_s_y, "Ts2", "Ts2", 0);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i57;
  int32_T c3_i58;
  int32_T c3_i59;
  real_T c3_b_inData[9];
  int32_T c3_i60;
  int32_T c3_i61;
  int32_T c3_i62;
  real_T c3_u[9];
  const mxArray *c3_y = NULL;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i57 = 0;
  for (c3_i58 = 0; c3_i58 < 3; c3_i58++) {
    for (c3_i59 = 0; c3_i59 < 3; c3_i59++) {
      c3_b_inData[c3_i59 + c3_i57] = (*(real_T (*)[9])c3_inData)[c3_i59 + c3_i57];
    }

    c3_i57 += 3;
  }

  c3_i60 = 0;
  for (c3_i61 = 0; c3_i61 < 3; c3_i61++) {
    for (c3_i62 = 0; c3_i62 < 3; c3_i62++) {
      c3_u[c3_i62 + c3_i60] = c3_b_inData[c3_i62 + c3_i60];
    }

    c3_i60 += 3;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_d_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[9])
{
  real_T c3_dv6[9];
  int32_T c3_i63;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv6, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c3_i63 = 0; c3_i63 < 9; c3_i63++) {
    c3_y[c3_i63] = c3_dv6[c3_i63];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_R;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[9];
  int32_T c3_i64;
  int32_T c3_i65;
  int32_T c3_i66;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_R = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_R), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_R);
  c3_i64 = 0;
  for (c3_i65 = 0; c3_i65 < 3; c3_i65++) {
    for (c3_i66 = 0; c3_i66 < 3; c3_i66++) {
      (*(real_T (*)[9])c3_outData)[c3_i66 + c3_i64] = c3_y[c3_i66 + c3_i64];
    }

    c3_i64 += 3;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

const mxArray *sf_c3_simulation1_get_eml_resolved_functions_info(void)
{
  const mxArray *c3_nameCaptureInfo = NULL;
  c3_nameCaptureInfo = NULL;
  sf_mex_assign(&c3_nameCaptureInfo, sf_mex_createstruct("structure", 2, 17, 1),
                false);
  c3_info_helper(&c3_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c3_nameCaptureInfo);
  return c3_nameCaptureInfo;
}

static void c3_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs0 = NULL;
  const mxArray *c3_lhs0 = NULL;
  const mxArray *c3_rhs1 = NULL;
  const mxArray *c3_lhs1 = NULL;
  const mxArray *c3_rhs2 = NULL;
  const mxArray *c3_lhs2 = NULL;
  const mxArray *c3_rhs3 = NULL;
  const mxArray *c3_lhs3 = NULL;
  const mxArray *c3_rhs4 = NULL;
  const mxArray *c3_lhs4 = NULL;
  const mxArray *c3_rhs5 = NULL;
  const mxArray *c3_lhs5 = NULL;
  const mxArray *c3_rhs6 = NULL;
  const mxArray *c3_lhs6 = NULL;
  const mxArray *c3_rhs7 = NULL;
  const mxArray *c3_lhs7 = NULL;
  const mxArray *c3_rhs8 = NULL;
  const mxArray *c3_lhs8 = NULL;
  const mxArray *c3_rhs9 = NULL;
  const mxArray *c3_lhs9 = NULL;
  const mxArray *c3_rhs10 = NULL;
  const mxArray *c3_lhs10 = NULL;
  const mxArray *c3_rhs11 = NULL;
  const mxArray *c3_lhs11 = NULL;
  const mxArray *c3_rhs12 = NULL;
  const mxArray *c3_lhs12 = NULL;
  const mxArray *c3_rhs13 = NULL;
  const mxArray *c3_lhs13 = NULL;
  const mxArray *c3_rhs14 = NULL;
  const mxArray *c3_lhs14 = NULL;
  const mxArray *c3_rhs15 = NULL;
  const mxArray *c3_lhs15 = NULL;
  const mxArray *c3_rhs16 = NULL;
  const mxArray *c3_lhs16 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c3_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286818722U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c3_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sin"), "name", "name", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395328504U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c3_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286818736U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c3_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c3_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c3_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c3_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c3_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c3_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgemm"), "name", "name", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375980690U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c3_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c3_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c3_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c3_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c3_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c3_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c3_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c3_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs16), "lhs", "lhs",
                  16);
  sf_mex_destroy(&c3_rhs0);
  sf_mex_destroy(&c3_lhs0);
  sf_mex_destroy(&c3_rhs1);
  sf_mex_destroy(&c3_lhs1);
  sf_mex_destroy(&c3_rhs2);
  sf_mex_destroy(&c3_lhs2);
  sf_mex_destroy(&c3_rhs3);
  sf_mex_destroy(&c3_lhs3);
  sf_mex_destroy(&c3_rhs4);
  sf_mex_destroy(&c3_lhs4);
  sf_mex_destroy(&c3_rhs5);
  sf_mex_destroy(&c3_lhs5);
  sf_mex_destroy(&c3_rhs6);
  sf_mex_destroy(&c3_lhs6);
  sf_mex_destroy(&c3_rhs7);
  sf_mex_destroy(&c3_lhs7);
  sf_mex_destroy(&c3_rhs8);
  sf_mex_destroy(&c3_lhs8);
  sf_mex_destroy(&c3_rhs9);
  sf_mex_destroy(&c3_lhs9);
  sf_mex_destroy(&c3_rhs10);
  sf_mex_destroy(&c3_lhs10);
  sf_mex_destroy(&c3_rhs11);
  sf_mex_destroy(&c3_lhs11);
  sf_mex_destroy(&c3_rhs12);
  sf_mex_destroy(&c3_lhs12);
  sf_mex_destroy(&c3_rhs13);
  sf_mex_destroy(&c3_lhs13);
  sf_mex_destroy(&c3_rhs14);
  sf_mex_destroy(&c3_lhs14);
  sf_mex_destroy(&c3_rhs15);
  sf_mex_destroy(&c3_lhs15);
  sf_mex_destroy(&c3_rhs16);
  sf_mex_destroy(&c3_lhs16);
}

static const mxArray *c3_emlrt_marshallOut(const char * c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c3_u)), false);
  return c3_y;
}

static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 7, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static void c3_eml_scalar_eg(SFc3_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_threshold(SFc3_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_b_eml_scalar_eg(SFc3_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(int32_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static int32_T c3_e_emlrt_marshallIn(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  int32_T c3_y;
  int32_T c3_i67;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_i67, 1, 6, 0U, 0, 0U, 0);
  c3_y = c3_i67;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_sfEvent;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  int32_T c3_y;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_b_sfEvent = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_sfEvent),
    &c3_thisId);
  sf_mex_destroy(&c3_b_sfEvent);
  *(int32_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static void c3_f_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  c3_sdz1UxT0ktEtEZMZuZvzEsH *c3_y)
{
  emlrtMsgIdentifier c3_thisId;
  static const char * c3_fieldNames[18] = { "thrustLimits", "thrustRateLimits",
    "thrustDir", "L", "nL", "rad", "bladeRad", "mass", "I", "km", "kf", "K",
    "velocityGain", "vGain", "zGain", "thetaGain", "omegaGain", "Ts2" };

  c3_thisId.fParent = c3_parentId;
  sf_mex_check_struct(c3_parentId, c3_u, 18, c3_fieldNames, 0U, NULL);
  c3_thisId.fIdentifier = "thrustLimits";
  c3_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u,
    "thrustLimits", "thrustLimits", 0)), &c3_thisId, c3_y->thrustLimits);
  c3_thisId.fIdentifier = "thrustRateLimits";
  c3_h_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u,
    "thrustRateLimits", "thrustRateLimits", 0)), &c3_thisId,
                        c3_y->thrustRateLimits);
  c3_thisId.fIdentifier = "thrustDir";
  c3_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u,
    "thrustDir", "thrustDir", 0)), &c3_thisId, c3_y->thrustDir);
  c3_thisId.fIdentifier = "L";
  c3_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u, "L", "L",
    0)), &c3_thisId, c3_y->L);
  c3_thisId.fIdentifier = "nL";
  c3_y->nL = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c3_u, "nL", "nL", 0)), &c3_thisId);
  c3_thisId.fIdentifier = "rad";
  c3_y->rad = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c3_u, "rad", "rad", 0)), &c3_thisId);
  c3_thisId.fIdentifier = "bladeRad";
  c3_y->bladeRad = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c3_u, "bladeRad", "bladeRad", 0)), &c3_thisId);
  c3_thisId.fIdentifier = "mass";
  c3_y->mass = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c3_u, "mass", "mass", 0)), &c3_thisId);
  c3_thisId.fIdentifier = "I";
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u, "I", "I",
    0)), &c3_thisId, c3_y->I);
  c3_thisId.fIdentifier = "km";
  c3_y->km = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c3_u, "km", "km", 0)), &c3_thisId);
  c3_thisId.fIdentifier = "kf";
  c3_y->kf = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c3_u, "kf", "kf", 0)), &c3_thisId);
  c3_thisId.fIdentifier = "K";
  c3_j_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u, "K", "K",
    0)), &c3_thisId, c3_y->K);
  c3_thisId.fIdentifier = "velocityGain";
  c3_y->velocityGain = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c3_u, "velocityGain", "velocityGain", 0)), &c3_thisId);
  c3_thisId.fIdentifier = "vGain";
  c3_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u, "vGain",
    "vGain", 0)), &c3_thisId, c3_y->vGain);
  c3_thisId.fIdentifier = "zGain";
  c3_y->zGain = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c3_u, "zGain", "zGain", 0)), &c3_thisId);
  c3_thisId.fIdentifier = "thetaGain";
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u,
    "thetaGain", "thetaGain", 0)), &c3_thisId, c3_y->thetaGain);
  c3_thisId.fIdentifier = "omegaGain";
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c3_u,
    "omegaGain", "omegaGain", 0)), &c3_thisId, c3_y->omegaGain);
  c3_thisId.fIdentifier = "Ts2";
  c3_y->Ts2 = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c3_u, "Ts2", "Ts2", 0)), &c3_thisId);
  sf_mex_destroy(&c3_u);
}

static void c3_g_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[8])
{
  real_T c3_dv7[8];
  int32_T c3_i68;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv7, 1, 0, 0U, 1, 0U, 2, 2, 4);
  for (c3_i68 = 0; c3_i68 < 8; c3_i68++) {
    c3_y[c3_i68] = c3_dv7[c3_i68];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_h_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[4])
{
  real_T c3_dv8[4];
  int32_T c3_i69;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv8, 1, 0, 0U, 1, 0U, 2, 1, 4);
  for (c3_i69 = 0; c3_i69 < 4; c3_i69++) {
    c3_y[c3_i69] = c3_dv8[c3_i69];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_i_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[12])
{
  real_T c3_dv9[12];
  int32_T c3_i70;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv9, 1, 0, 0U, 1, 0U, 2, 3, 4);
  for (c3_i70 = 0; c3_i70 < 12; c3_i70++) {
    c3_y[c3_i70] = c3_dv9[c3_i70];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_j_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[16])
{
  real_T c3_dv10[16];
  int32_T c3_i71;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv10, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c3_i71 = 0; c3_i71 < 16; c3_i71++) {
    c3_y[c3_i71] = c3_dv10[c3_i71];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_k_emlrt_marshallIn(SFc3_simulation1InstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[2])
{
  real_T c3_dv11[2];
  int32_T c3_i72;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv11, 1, 0, 0U, 1, 0U, 1, 2);
  for (c3_i72 = 0; c3_i72 < 2; c3_i72++) {
    c3_y[c3_i72] = c3_dv11[c3_i72];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_quad;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  c3_sdz1UxT0ktEtEZMZuZvzEsH c3_y;
  SFc3_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc3_simulation1InstanceStruct *)chartInstanceVoid;
  c3_b_quad = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_quad), &c3_thisId, &c3_y);
  sf_mex_destroy(&c3_b_quad);
  *(c3_sdz1UxT0ktEtEZMZuZvzEsH *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static uint8_T c3_l_emlrt_marshallIn(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_b_is_active_c3_simulation1, const char_T
  *c3_identifier)
{
  uint8_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_m_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_b_is_active_c3_simulation1), &c3_thisId);
  sf_mex_destroy(&c3_b_is_active_c3_simulation1);
  return c3_y;
}

static uint8_T c3_m_emlrt_marshallIn(SFc3_simulation1InstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  uint8_T c3_y;
  uint8_T c3_u0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_u0, 1, 3, 0U, 0, 0U, 0);
  c3_y = c3_u0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void init_dsm_address_info(SFc3_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc3_simulation1InstanceStruct
  *chartInstance)
{
  chartInstance->c3_uTot = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c3_dist = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c3_vDot = (real_T (*)[3])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c3_theta = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 2);
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

void sf_c3_simulation1_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1138095773U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2941551444U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2193808057U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(170421775U);
}

mxArray* sf_c3_simulation1_get_post_codegen_info(void);
mxArray *sf_c3_simulation1_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("D6Zn6x7S68uIwTCNeYfvoH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
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
    mxArray* mxPostCodegenInfo = sf_c3_simulation1_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c3_simulation1_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c3_simulation1_jit_fallback_info(void)
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

mxArray *sf_c3_simulation1_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c3_simulation1_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c3_simulation1(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"vDot\",},{M[8],M[0],T\"is_active_c3_simulation1\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c3_simulation1_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc3_simulation1InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc3_simulation1InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulation1MachineNumber_,
           3,
           1,
           1,
           0,
           5,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation(_simulation1MachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_simulation1MachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _simulation1MachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"uTot");
          _SFD_SET_DATA_PROPS(1,1,1,0,"dist");
          _SFD_SET_DATA_PROPS(2,2,0,1,"vDot");
          _SFD_SET_DATA_PROPS(3,1,1,0,"theta");
          _SFD_SET_DATA_PROPS(4,10,0,0,"quad");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,486);
        _SFD_CV_INIT_EML_IF(0,1,0,351,359,407,432);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_STRUCT,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_c_sf_marshallOut,(MexInFcnForType)c3_e_sf_marshallIn);
        _SFD_SET_DATA_VALUE_PTR(0U, chartInstance->c3_uTot);
        _SFD_SET_DATA_VALUE_PTR(1U, chartInstance->c3_dist);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c3_vDot);
        _SFD_SET_DATA_VALUE_PTR(3U, *chartInstance->c3_theta);
        _SFD_SET_DATA_VALUE_PTR(4U, &chartInstance->c3_quad);
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _simulation1MachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "w8yNkoQGLbWkf38zLGJD6C";
}

static void sf_opaque_initialize_c3_simulation1(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc3_simulation1InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c3_simulation1((SFc3_simulation1InstanceStruct*)
    chartInstanceVar);
  initialize_c3_simulation1((SFc3_simulation1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c3_simulation1(void *chartInstanceVar)
{
  enable_c3_simulation1((SFc3_simulation1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c3_simulation1(void *chartInstanceVar)
{
  disable_c3_simulation1((SFc3_simulation1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c3_simulation1(void *chartInstanceVar)
{
  sf_gateway_c3_simulation1((SFc3_simulation1InstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c3_simulation1(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c3_simulation1((SFc3_simulation1InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c3_simulation1(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c3_simulation1((SFc3_simulation1InstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c3_simulation1(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc3_simulation1InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation1_optimization_info();
    }

    finalize_c3_simulation1((SFc3_simulation1InstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc3_simulation1((SFc3_simulation1InstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c3_simulation1(SimStruct *S)
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
    initialize_params_c3_simulation1((SFc3_simulation1InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c3_simulation1(SimStruct *S)
{
  /* Actual parameters from chart:
     quad
   */
  const char_T *rtParamNames[] = { "quad" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation1_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,3);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,3,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,3,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,3);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,3,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,3,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,3);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(4193405058U));
  ssSetChecksum1(S,(1042939635U));
  ssSetChecksum2(S,(802959368U));
  ssSetChecksum3(S,(120983841U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c3_simulation1(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c3_simulation1(SimStruct *S)
{
  SFc3_simulation1InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc3_simulation1InstanceStruct *)utMalloc(sizeof
    (SFc3_simulation1InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc3_simulation1InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c3_simulation1;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c3_simulation1;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c3_simulation1;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c3_simulation1;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c3_simulation1;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c3_simulation1;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c3_simulation1;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c3_simulation1;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c3_simulation1;
  chartInstance->chartInfo.mdlStart = mdlStart_c3_simulation1;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c3_simulation1;
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

void c3_simulation1_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c3_simulation1(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c3_simulation1(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c3_simulation1(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c3_simulation1_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
