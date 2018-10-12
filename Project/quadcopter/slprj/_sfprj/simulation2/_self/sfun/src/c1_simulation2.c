/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation2_sfun.h"
#include "c1_simulation2.h"
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
static const char * c1_debug_family_names[11] = { "roll", "pitch", "yaw", "R",
  "g", "quad", "nargin", "nargout", "uTot", "theta", "vDot" };

/* Function Declarations */
static void initialize_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance);
static void initialize_params_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance);
static void enable_c1_simulation2(SFc1_simulation2InstanceStruct *chartInstance);
static void disable_c1_simulation2(SFc1_simulation2InstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_simulation2
  (SFc1_simulation2InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c1_simulation2
  (SFc1_simulation2InstanceStruct *chartInstance);
static void set_sim_state_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_st);
static void finalize_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance);
static void sf_gateway_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance);
static void mdl_start_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance);
static void c1_chartstep_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance);
static void initSimStructsc1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_b_vDot, const char_T *c1_identifier, real_T c1_y[3]);
static void c1_b_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3]);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static real_T c1_c_emlrt_marshallIn(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_d_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[9]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_info_helper(const mxArray **c1_info);
static const mxArray *c1_emlrt_marshallOut(const char * c1_u);
static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u);
static void c1_eml_scalar_eg(SFc1_simulation2InstanceStruct *chartInstance);
static void c1_threshold(SFc1_simulation2InstanceStruct *chartInstance);
static void c1_b_eml_scalar_eg(SFc1_simulation2InstanceStruct *chartInstance);
static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_e_emlrt_marshallIn(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_f_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  c1_sdz1UxT0ktEtEZMZuZvzEsH *c1_y);
static void c1_g_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[8]);
static void c1_h_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[4]);
static void c1_i_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[12]);
static void c1_j_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[16]);
static void c1_k_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[2]);
static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_l_emlrt_marshallIn(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_b_is_active_c1_simulation2, const char_T
  *c1_identifier);
static uint8_T c1_m_emlrt_marshallIn(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void init_dsm_address_info(SFc1_simulation2InstanceStruct *chartInstance);
static void init_simulink_io_address(SFc1_simulation2InstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance)
{
  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_simulation2 = 0U;
}

static void initialize_params_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance)
{
  const mxArray *c1_m0 = NULL;
  const mxArray *c1_mxField;
  c1_sdz1UxT0ktEtEZMZuZvzEsH c1_r0;
  c1_m0 = sf_mex_get_sfun_param(chartInstance->S, 0, 1);
  c1_mxField = sf_mex_getfield(c1_m0, "thrustLimits", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.thrustLimits, 1, 0,
                      0U, 1, 0U, 2, 2, 4);
  c1_mxField = sf_mex_getfield(c1_m0, "thrustRateLimits", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.thrustRateLimits, 1,
                      0, 0U, 1, 0U, 2, 1, 4);
  c1_mxField = sf_mex_getfield(c1_m0, "thrustDir", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.thrustDir, 1, 0, 0U,
                      1, 0U, 2, 3, 4);
  c1_mxField = sf_mex_getfield(c1_m0, "L", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.L, 1, 0, 0U, 1, 0U,
                      2, 3, 4);
  c1_mxField = sf_mex_getfield(c1_m0, "nL", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.nL, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "rad", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.rad, 1, 0, 0U, 0,
                      0U, 0);
  c1_mxField = sf_mex_getfield(c1_m0, "bladeRad", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.bladeRad, 1, 0, 0U,
                      0, 0U, 0);
  c1_mxField = sf_mex_getfield(c1_m0, "mass", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.mass, 1, 0, 0U, 0,
                      0U, 0);
  c1_mxField = sf_mex_getfield(c1_m0, "I", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.I, 1, 0, 0U, 1, 0U,
                      2, 3, 3);
  c1_mxField = sf_mex_getfield(c1_m0, "km", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.km, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "kf", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.kf, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "K", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.K, 1, 0, 0U, 1, 0U,
                      2, 4, 4);
  c1_mxField = sf_mex_getfield(c1_m0, "velocityGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.velocityGain, 1, 0,
                      0U, 0, 0U, 0);
  c1_mxField = sf_mex_getfield(c1_m0, "vGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.vGain, 1, 0, 0U, 1,
                      0U, 1, 2);
  c1_mxField = sf_mex_getfield(c1_m0, "zGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.zGain, 1, 0, 0U, 0,
                      0U, 0);
  c1_mxField = sf_mex_getfield(c1_m0, "thetaGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.thetaGain, 1, 0, 0U,
                      1, 0U, 1, 3);
  c1_mxField = sf_mex_getfield(c1_m0, "omegaGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), c1_r0.omegaGain, 1, 0, 0U,
                      1, 0U, 1, 3);
  c1_mxField = sf_mex_getfield(c1_m0, "Ts2", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c1_mxField), &c1_r0.Ts2, 1, 0, 0U, 0,
                      0U, 0);
  sf_mex_destroy(&c1_m0);
  chartInstance->c1_quad = c1_r0;
}

static void enable_c1_simulation2(SFc1_simulation2InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_simulation2(SFc1_simulation2InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_simulation2
  (SFc1_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_simulation2
  (SFc1_simulation2InstanceStruct *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  int32_T c1_i0;
  real_T c1_u[3];
  const mxArray *c1_b_y = NULL;
  uint8_T c1_hoistedGlobal;
  uint8_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellmatrix(2, 1), false);
  for (c1_i0 = 0; c1_i0 < 3; c1_i0++) {
    c1_u[c1_i0] = (*chartInstance->c1_vDot)[c1_i0];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  c1_hoistedGlobal = chartInstance->c1_is_active_c1_simulation2;
  c1_b_u = c1_hoistedGlobal;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  sf_mex_assign(&c1_st, c1_y, false);
  return c1_st;
}

static void set_sim_state_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_st)
{
  const mxArray *c1_u;
  real_T c1_dv0[3];
  int32_T c1_i1;
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_u = sf_mex_dup(c1_st);
  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 0)), "vDot",
                      c1_dv0);
  for (c1_i1 = 0; c1_i1 < 3; c1_i1++) {
    (*chartInstance->c1_vDot)[c1_i1] = c1_dv0[c1_i1];
  }

  chartInstance->c1_is_active_c1_simulation2 = c1_l_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 1)),
     "is_active_c1_simulation2");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_simulation2(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance)
{
  int32_T c1_i2;
  int32_T c1_i3;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c1_uTot, 0U);
  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_simulation2(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulation2MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c1_i2 = 0; c1_i2 < 3; c1_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_vDot)[c1_i2], 1U);
  }

  for (c1_i3 = 0; c1_i3 < 3; c1_i3++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_theta)[c1_i3], 2U);
  }
}

static void mdl_start_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_chartstep_c1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance)
{
  real_T c1_hoistedGlobal;
  real_T c1_b_uTot;
  int32_T c1_i4;
  real_T c1_b_theta[3];
  uint32_T c1_debug_family_var_map[11];
  real_T c1_roll;
  real_T c1_pitch;
  real_T c1_yaw;
  real_T c1_R[9];
  real_T c1_g;
  c1_sdz1UxT0ktEtEZMZuZvzEsH c1_b_quad;
  real_T c1_nargin = 3.0;
  real_T c1_nargout = 1.0;
  real_T c1_b_vDot[3];
  static c1_sdz1UxT0ktEtEZMZuZvzEsH c1_r1 = { { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
      0.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 1.0 }, { 0.0, 0.2, 0.0, -0.2, 0.0, 0.0, 0.0, -0.2, 0.0,
      0.2, 0.0, 0.0 }, 0.2, 0.040000000000000008, 0.080000000000000016, 8.0, {
      10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 15.0 }, 11.0, 28.0, { 28.0, 0.0,
      -5.6000000000000005, 11.0, 28.0, 5.6000000000000005, 0.0, -11.0, 28.0, 0.0,
      5.6000000000000005, 11.0, 28.0, -5.6000000000000005, 0.0, -11.0 }, 10.0, {
      0.1, 0.1 }, 1.0, { 1.0, 1.0, 1.0 }, { 0.2, 0.2, 0.2 }, 0.3 };

  real_T c1_x;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_h_x;
  int32_T c1_i5;
  int32_T c1_i6;
  static real_T c1_dv1[3] = { 1.0, 0.0, 0.0 };

  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_m_x;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_p_x;
  int32_T c1_i7;
  real_T c1_a[9];
  real_T c1_b[9];
  int32_T c1_i8;
  int32_T c1_i9;
  static real_T c1_dv2[3] = { 0.0, 1.0, 0.0 };

  int32_T c1_i10;
  int32_T c1_i11;
  int32_T c1_i12;
  real_T c1_C[9];
  int32_T c1_i13;
  int32_T c1_i14;
  int32_T c1_i15;
  int32_T c1_i16;
  int32_T c1_i17;
  int32_T c1_i18;
  int32_T c1_i19;
  int32_T c1_i20;
  real_T c1_q_x;
  real_T c1_r_x;
  real_T c1_s_x;
  real_T c1_t_x;
  real_T c1_u_x;
  real_T c1_v_x;
  real_T c1_w_x;
  real_T c1_x_x;
  int32_T c1_i21;
  int32_T c1_i22;
  int32_T c1_i23;
  static real_T c1_dv3[3] = { 0.0, 0.0, 1.0 };

  int32_T c1_i24;
  int32_T c1_i25;
  int32_T c1_i26;
  int32_T c1_i27;
  int32_T c1_i28;
  int32_T c1_i29;
  int32_T c1_i30;
  int32_T c1_i31;
  int32_T c1_i32;
  int32_T c1_i33;
  int32_T c1_i34;
  real_T c1_b_a;
  int32_T c1_i35;
  int32_T c1_i36;
  int32_T c1_i37;
  real_T c1_y[3];
  int32_T c1_i38;
  int32_T c1_i39;
  static real_T c1_b_b[3] = { 0.0, 0.0, 1.0 };

  int32_T c1_i40;
  static real_T c1_b_y[3] = { -0.0, -0.0, -78.48 };

  int32_T c1_i41;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  c1_hoistedGlobal = *chartInstance->c1_uTot;
  c1_b_uTot = c1_hoistedGlobal;
  for (c1_i4 = 0; c1_i4 < 3; c1_i4++) {
    c1_b_theta[c1_i4] = (*chartInstance->c1_theta)[c1_i4];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 11U, 11U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_roll, 0U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_pitch, 1U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_yaw, 2U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_R, 3U, c1_d_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_g, 4U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_b_quad, 5U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 6U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 7U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_b_uTot, 8U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_theta, 9U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_vDot, 10U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  c1_b_quad = c1_r1;
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 5);
  c1_roll = c1_b_theta[0];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 5);
  c1_pitch = c1_b_theta[1];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 5);
  c1_yaw = c1_b_theta[2];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 6);
  c1_x = c1_roll;
  c1_b_x = c1_x;
  c1_b_x = muDoubleScalarCos(c1_b_x);
  c1_c_x = c1_roll;
  c1_d_x = c1_c_x;
  c1_d_x = muDoubleScalarSin(c1_d_x);
  c1_e_x = c1_roll;
  c1_f_x = c1_e_x;
  c1_f_x = muDoubleScalarSin(c1_f_x);
  c1_g_x = c1_roll;
  c1_h_x = c1_g_x;
  c1_h_x = muDoubleScalarCos(c1_h_x);
  c1_i5 = 0;
  for (c1_i6 = 0; c1_i6 < 3; c1_i6++) {
    c1_R[c1_i5] = c1_dv1[c1_i6];
    c1_i5 += 3;
  }

  c1_R[1] = 0.0;
  c1_R[4] = c1_b_x;
  c1_R[7] = -c1_d_x;
  c1_R[2] = 0.0;
  c1_R[5] = c1_f_x;
  c1_R[8] = c1_h_x;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 7);
  c1_i_x = c1_pitch;
  c1_j_x = c1_i_x;
  c1_j_x = muDoubleScalarCos(c1_j_x);
  c1_k_x = c1_pitch;
  c1_l_x = c1_k_x;
  c1_l_x = muDoubleScalarSin(c1_l_x);
  c1_m_x = c1_pitch;
  c1_n_x = c1_m_x;
  c1_n_x = muDoubleScalarSin(c1_n_x);
  c1_o_x = c1_pitch;
  c1_p_x = c1_o_x;
  c1_p_x = muDoubleScalarCos(c1_p_x);
  for (c1_i7 = 0; c1_i7 < 9; c1_i7++) {
    c1_a[c1_i7] = c1_R[c1_i7];
  }

  c1_b[0] = c1_j_x;
  c1_b[3] = 0.0;
  c1_b[6] = c1_l_x;
  c1_i8 = 0;
  for (c1_i9 = 0; c1_i9 < 3; c1_i9++) {
    c1_b[c1_i8 + 1] = c1_dv2[c1_i9];
    c1_i8 += 3;
  }

  c1_b[2] = -c1_n_x;
  c1_b[5] = 0.0;
  c1_b[8] = c1_p_x;
  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  for (c1_i10 = 0; c1_i10 < 9; c1_i10++) {
    c1_R[c1_i10] = 0.0;
  }

  for (c1_i11 = 0; c1_i11 < 9; c1_i11++) {
    c1_R[c1_i11] = 0.0;
  }

  for (c1_i12 = 0; c1_i12 < 9; c1_i12++) {
    c1_C[c1_i12] = c1_R[c1_i12];
  }

  for (c1_i13 = 0; c1_i13 < 9; c1_i13++) {
    c1_R[c1_i13] = c1_C[c1_i13];
  }

  c1_threshold(chartInstance);
  for (c1_i14 = 0; c1_i14 < 9; c1_i14++) {
    c1_C[c1_i14] = c1_R[c1_i14];
  }

  for (c1_i15 = 0; c1_i15 < 9; c1_i15++) {
    c1_R[c1_i15] = c1_C[c1_i15];
  }

  for (c1_i16 = 0; c1_i16 < 3; c1_i16++) {
    c1_i17 = 0;
    for (c1_i18 = 0; c1_i18 < 3; c1_i18++) {
      c1_R[c1_i17 + c1_i16] = 0.0;
      c1_i19 = 0;
      for (c1_i20 = 0; c1_i20 < 3; c1_i20++) {
        c1_R[c1_i17 + c1_i16] += c1_a[c1_i19 + c1_i16] * c1_b[c1_i20 + c1_i17];
        c1_i19 += 3;
      }

      c1_i17 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 8);
  c1_q_x = c1_yaw;
  c1_r_x = c1_q_x;
  c1_r_x = muDoubleScalarCos(c1_r_x);
  c1_s_x = c1_yaw;
  c1_t_x = c1_s_x;
  c1_t_x = muDoubleScalarSin(c1_t_x);
  c1_u_x = c1_yaw;
  c1_v_x = c1_u_x;
  c1_v_x = muDoubleScalarSin(c1_v_x);
  c1_w_x = c1_yaw;
  c1_x_x = c1_w_x;
  c1_x_x = muDoubleScalarCos(c1_x_x);
  for (c1_i21 = 0; c1_i21 < 9; c1_i21++) {
    c1_a[c1_i21] = c1_R[c1_i21];
  }

  c1_b[0] = c1_r_x;
  c1_b[3] = -c1_t_x;
  c1_b[6] = 0.0;
  c1_b[1] = c1_v_x;
  c1_b[4] = c1_x_x;
  c1_b[7] = 0.0;
  c1_i22 = 0;
  for (c1_i23 = 0; c1_i23 < 3; c1_i23++) {
    c1_b[c1_i22 + 2] = c1_dv3[c1_i23];
    c1_i22 += 3;
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  for (c1_i24 = 0; c1_i24 < 9; c1_i24++) {
    c1_R[c1_i24] = 0.0;
  }

  for (c1_i25 = 0; c1_i25 < 9; c1_i25++) {
    c1_R[c1_i25] = 0.0;
  }

  for (c1_i26 = 0; c1_i26 < 9; c1_i26++) {
    c1_C[c1_i26] = c1_R[c1_i26];
  }

  for (c1_i27 = 0; c1_i27 < 9; c1_i27++) {
    c1_R[c1_i27] = c1_C[c1_i27];
  }

  c1_threshold(chartInstance);
  for (c1_i28 = 0; c1_i28 < 9; c1_i28++) {
    c1_C[c1_i28] = c1_R[c1_i28];
  }

  for (c1_i29 = 0; c1_i29 < 9; c1_i29++) {
    c1_R[c1_i29] = c1_C[c1_i29];
  }

  for (c1_i30 = 0; c1_i30 < 3; c1_i30++) {
    c1_i31 = 0;
    for (c1_i32 = 0; c1_i32 < 3; c1_i32++) {
      c1_R[c1_i31 + c1_i30] = 0.0;
      c1_i33 = 0;
      for (c1_i34 = 0; c1_i34 < 3; c1_i34++) {
        c1_R[c1_i31 + c1_i30] += c1_a[c1_i33 + c1_i30] * c1_b[c1_i34 + c1_i31];
        c1_i33 += 3;
      }

      c1_i31 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 11);
  c1_g = 9.81;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 12);
  c1_b_a = c1_b_uTot;
  for (c1_i35 = 0; c1_i35 < 9; c1_i35++) {
    c1_b[c1_i35] = c1_R[c1_i35];
  }

  for (c1_i36 = 0; c1_i36 < 9; c1_i36++) {
    c1_b[c1_i36] *= c1_b_a;
  }

  c1_b_eml_scalar_eg(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i37 = 0; c1_i37 < 3; c1_i37++) {
    c1_y[c1_i37] = 0.0;
    c1_i38 = 0;
    for (c1_i39 = 0; c1_i39 < 3; c1_i39++) {
      c1_y[c1_i37] += c1_b[c1_i38 + c1_i37] * c1_b_b[c1_i39];
      c1_i38 += 3;
    }
  }

  for (c1_i40 = 0; c1_i40 < 3; c1_i40++) {
    c1_b_vDot[c1_i40] = c1_b_y[c1_i40] + c1_y[c1_i40];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -12);
  _SFD_SYMBOL_SCOPE_POP();
  for (c1_i41 = 0; c1_i41 < 3; c1_i41++) {
    (*chartInstance->c1_vDot)[c1_i41] = c1_b_vDot[c1_i41];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_simulation2(SFc1_simulation2InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber)
{
  (void)c1_machineNumber;
  (void)c1_chartNumber;
  (void)c1_instanceNumber;
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i42;
  real_T c1_b_inData[3];
  int32_T c1_i43;
  real_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i42 = 0; c1_i42 < 3; c1_i42++) {
    c1_b_inData[c1_i42] = (*(real_T (*)[3])c1_inData)[c1_i42];
  }

  for (c1_i43 = 0; c1_i43 < 3; c1_i43++) {
    c1_u[c1_i43] = c1_b_inData[c1_i43];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_b_vDot, const char_T *c1_identifier, real_T c1_y[3])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_vDot), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_vDot);
}

static void c1_b_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3])
{
  real_T c1_dv4[3];
  int32_T c1_i44;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv4, 1, 0, 0U, 1, 0U, 1, 3);
  for (c1_i44 = 0; c1_i44 < 3; c1_i44++) {
    c1_y[c1_i44] = c1_dv4[c1_i44];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_vDot;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[3];
  int32_T c1_i45;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_b_vDot = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_vDot), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_vDot);
  for (c1_i45 = 0; c1_i45 < 3; c1_i45++) {
    (*(real_T (*)[3])c1_outData)[c1_i45] = c1_y[c1_i45];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(real_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static real_T c1_c_emlrt_marshallIn(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d0, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_nargout;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_nargout = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_nargout), &c1_thisId);
  sf_mex_destroy(&c1_nargout);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  c1_sdz1UxT0ktEtEZMZuZvzEsH c1_u;
  const mxArray *c1_y = NULL;
  int32_T c1_i46;
  real_T c1_b_u[8];
  const mxArray *c1_b_y = NULL;
  int32_T c1_i47;
  real_T c1_c_u[4];
  const mxArray *c1_c_y = NULL;
  int32_T c1_i48;
  real_T c1_d_u[12];
  const mxArray *c1_d_y = NULL;
  int32_T c1_i49;
  real_T c1_e_u[12];
  const mxArray *c1_e_y = NULL;
  real_T c1_f_u;
  const mxArray *c1_f_y = NULL;
  real_T c1_g_u;
  const mxArray *c1_g_y = NULL;
  real_T c1_h_u;
  const mxArray *c1_h_y = NULL;
  real_T c1_i_u;
  const mxArray *c1_i_y = NULL;
  int32_T c1_i50;
  real_T c1_j_u[9];
  const mxArray *c1_j_y = NULL;
  real_T c1_k_u;
  const mxArray *c1_k_y = NULL;
  real_T c1_l_u;
  const mxArray *c1_l_y = NULL;
  int32_T c1_i51;
  real_T c1_m_u[16];
  const mxArray *c1_m_y = NULL;
  real_T c1_n_u;
  const mxArray *c1_n_y = NULL;
  int32_T c1_i52;
  real_T c1_o_u[2];
  const mxArray *c1_o_y = NULL;
  real_T c1_p_u;
  const mxArray *c1_p_y = NULL;
  int32_T c1_i53;
  real_T c1_q_u[3];
  const mxArray *c1_q_y = NULL;
  int32_T c1_i54;
  real_T c1_r_u[3];
  const mxArray *c1_r_y = NULL;
  real_T c1_s_u;
  const mxArray *c1_s_y = NULL;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(c1_sdz1UxT0ktEtEZMZuZvzEsH *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  for (c1_i46 = 0; c1_i46 < 8; c1_i46++) {
    c1_b_u[c1_i46] = c1_u.thrustLimits[c1_i46];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 2, 4),
                false);
  sf_mex_addfield(c1_y, c1_b_y, "thrustLimits", "thrustLimits", 0);
  for (c1_i47 = 0; c1_i47 < 4; c1_i47++) {
    c1_c_u[c1_i47] = c1_u.thrustRateLimits[c1_i47];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_c_u, 0, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_addfield(c1_y, c1_c_y, "thrustRateLimits", "thrustRateLimits", 0);
  for (c1_i48 = 0; c1_i48 < 12; c1_i48++) {
    c1_d_u[c1_i48] = c1_u.thrustDir[c1_i48];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_d_u, 0, 0U, 1U, 0U, 2, 3, 4),
                false);
  sf_mex_addfield(c1_y, c1_d_y, "thrustDir", "thrustDir", 0);
  for (c1_i49 = 0; c1_i49 < 12; c1_i49++) {
    c1_e_u[c1_i49] = c1_u.L[c1_i49];
  }

  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", c1_e_u, 0, 0U, 1U, 0U, 2, 3, 4),
                false);
  sf_mex_addfield(c1_y, c1_e_y, "L", "L", 0);
  c1_f_u = c1_u.nL;
  c1_f_y = NULL;
  sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_f_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_f_y, "nL", "nL", 0);
  c1_g_u = c1_u.rad;
  c1_g_y = NULL;
  sf_mex_assign(&c1_g_y, sf_mex_create("y", &c1_g_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_g_y, "rad", "rad", 0);
  c1_h_u = c1_u.bladeRad;
  c1_h_y = NULL;
  sf_mex_assign(&c1_h_y, sf_mex_create("y", &c1_h_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_h_y, "bladeRad", "bladeRad", 0);
  c1_i_u = c1_u.mass;
  c1_i_y = NULL;
  sf_mex_assign(&c1_i_y, sf_mex_create("y", &c1_i_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_i_y, "mass", "mass", 0);
  for (c1_i50 = 0; c1_i50 < 9; c1_i50++) {
    c1_j_u[c1_i50] = c1_u.I[c1_i50];
  }

  c1_j_y = NULL;
  sf_mex_assign(&c1_j_y, sf_mex_create("y", c1_j_u, 0, 0U, 1U, 0U, 2, 3, 3),
                false);
  sf_mex_addfield(c1_y, c1_j_y, "I", "I", 0);
  c1_k_u = c1_u.km;
  c1_k_y = NULL;
  sf_mex_assign(&c1_k_y, sf_mex_create("y", &c1_k_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_k_y, "km", "km", 0);
  c1_l_u = c1_u.kf;
  c1_l_y = NULL;
  sf_mex_assign(&c1_l_y, sf_mex_create("y", &c1_l_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_l_y, "kf", "kf", 0);
  for (c1_i51 = 0; c1_i51 < 16; c1_i51++) {
    c1_m_u[c1_i51] = c1_u.K[c1_i51];
  }

  c1_m_y = NULL;
  sf_mex_assign(&c1_m_y, sf_mex_create("y", c1_m_u, 0, 0U, 1U, 0U, 2, 4, 4),
                false);
  sf_mex_addfield(c1_y, c1_m_y, "K", "K", 0);
  c1_n_u = c1_u.velocityGain;
  c1_n_y = NULL;
  sf_mex_assign(&c1_n_y, sf_mex_create("y", &c1_n_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_n_y, "velocityGain", "velocityGain", 0);
  for (c1_i52 = 0; c1_i52 < 2; c1_i52++) {
    c1_o_u[c1_i52] = c1_u.vGain[c1_i52];
  }

  c1_o_y = NULL;
  sf_mex_assign(&c1_o_y, sf_mex_create("y", c1_o_u, 0, 0U, 1U, 0U, 1, 2), false);
  sf_mex_addfield(c1_y, c1_o_y, "vGain", "vGain", 0);
  c1_p_u = c1_u.zGain;
  c1_p_y = NULL;
  sf_mex_assign(&c1_p_y, sf_mex_create("y", &c1_p_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_p_y, "zGain", "zGain", 0);
  for (c1_i53 = 0; c1_i53 < 3; c1_i53++) {
    c1_q_u[c1_i53] = c1_u.thetaGain[c1_i53];
  }

  c1_q_y = NULL;
  sf_mex_assign(&c1_q_y, sf_mex_create("y", c1_q_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_addfield(c1_y, c1_q_y, "thetaGain", "thetaGain", 0);
  for (c1_i54 = 0; c1_i54 < 3; c1_i54++) {
    c1_r_u[c1_i54] = c1_u.omegaGain[c1_i54];
  }

  c1_r_y = NULL;
  sf_mex_assign(&c1_r_y, sf_mex_create("y", c1_r_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_addfield(c1_y, c1_r_y, "omegaGain", "omegaGain", 0);
  c1_s_u = c1_u.Ts2;
  c1_s_y = NULL;
  sf_mex_assign(&c1_s_y, sf_mex_create("y", &c1_s_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_s_y, "Ts2", "Ts2", 0);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i55;
  int32_T c1_i56;
  int32_T c1_i57;
  real_T c1_b_inData[9];
  int32_T c1_i58;
  int32_T c1_i59;
  int32_T c1_i60;
  real_T c1_u[9];
  const mxArray *c1_y = NULL;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i55 = 0;
  for (c1_i56 = 0; c1_i56 < 3; c1_i56++) {
    for (c1_i57 = 0; c1_i57 < 3; c1_i57++) {
      c1_b_inData[c1_i57 + c1_i55] = (*(real_T (*)[9])c1_inData)[c1_i57 + c1_i55];
    }

    c1_i55 += 3;
  }

  c1_i58 = 0;
  for (c1_i59 = 0; c1_i59 < 3; c1_i59++) {
    for (c1_i60 = 0; c1_i60 < 3; c1_i60++) {
      c1_u[c1_i60 + c1_i58] = c1_b_inData[c1_i60 + c1_i58];
    }

    c1_i58 += 3;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_d_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[9])
{
  real_T c1_dv5[9];
  int32_T c1_i61;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv5, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c1_i61 = 0; c1_i61 < 9; c1_i61++) {
    c1_y[c1_i61] = c1_dv5[c1_i61];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_R;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[9];
  int32_T c1_i62;
  int32_T c1_i63;
  int32_T c1_i64;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_R = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_R), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_R);
  c1_i62 = 0;
  for (c1_i63 = 0; c1_i63 < 3; c1_i63++) {
    for (c1_i64 = 0; c1_i64 < 3; c1_i64++) {
      (*(real_T (*)[9])c1_outData)[c1_i64 + c1_i62] = c1_y[c1_i64 + c1_i62];
    }

    c1_i62 += 3;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_simulation2_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_createstruct("structure", 2, 17, 1),
                false);
  c1_info_helper(&c1_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c1_nameCaptureInfo);
  return c1_nameCaptureInfo;
}

static void c1_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs0 = NULL;
  const mxArray *c1_lhs0 = NULL;
  const mxArray *c1_rhs1 = NULL;
  const mxArray *c1_lhs1 = NULL;
  const mxArray *c1_rhs2 = NULL;
  const mxArray *c1_lhs2 = NULL;
  const mxArray *c1_rhs3 = NULL;
  const mxArray *c1_lhs3 = NULL;
  const mxArray *c1_rhs4 = NULL;
  const mxArray *c1_lhs4 = NULL;
  const mxArray *c1_rhs5 = NULL;
  const mxArray *c1_lhs5 = NULL;
  const mxArray *c1_rhs6 = NULL;
  const mxArray *c1_lhs6 = NULL;
  const mxArray *c1_rhs7 = NULL;
  const mxArray *c1_lhs7 = NULL;
  const mxArray *c1_rhs8 = NULL;
  const mxArray *c1_lhs8 = NULL;
  const mxArray *c1_rhs9 = NULL;
  const mxArray *c1_lhs9 = NULL;
  const mxArray *c1_rhs10 = NULL;
  const mxArray *c1_lhs10 = NULL;
  const mxArray *c1_rhs11 = NULL;
  const mxArray *c1_lhs11 = NULL;
  const mxArray *c1_rhs12 = NULL;
  const mxArray *c1_lhs12 = NULL;
  const mxArray *c1_rhs13 = NULL;
  const mxArray *c1_lhs13 = NULL;
  const mxArray *c1_rhs14 = NULL;
  const mxArray *c1_lhs14 = NULL;
  const mxArray *c1_rhs15 = NULL;
  const mxArray *c1_lhs15 = NULL;
  const mxArray *c1_rhs16 = NULL;
  const mxArray *c1_lhs16 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c1_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818722U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c1_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328504U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c1_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818736U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c1_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c1_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c1_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c1_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c1_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c1_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgemm"), "name", "name", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980690U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c1_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c1_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c1_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c1_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c1_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c1_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c1_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c1_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs16), "lhs", "lhs",
                  16);
  sf_mex_destroy(&c1_rhs0);
  sf_mex_destroy(&c1_lhs0);
  sf_mex_destroy(&c1_rhs1);
  sf_mex_destroy(&c1_lhs1);
  sf_mex_destroy(&c1_rhs2);
  sf_mex_destroy(&c1_lhs2);
  sf_mex_destroy(&c1_rhs3);
  sf_mex_destroy(&c1_lhs3);
  sf_mex_destroy(&c1_rhs4);
  sf_mex_destroy(&c1_lhs4);
  sf_mex_destroy(&c1_rhs5);
  sf_mex_destroy(&c1_lhs5);
  sf_mex_destroy(&c1_rhs6);
  sf_mex_destroy(&c1_lhs6);
  sf_mex_destroy(&c1_rhs7);
  sf_mex_destroy(&c1_lhs7);
  sf_mex_destroy(&c1_rhs8);
  sf_mex_destroy(&c1_lhs8);
  sf_mex_destroy(&c1_rhs9);
  sf_mex_destroy(&c1_lhs9);
  sf_mex_destroy(&c1_rhs10);
  sf_mex_destroy(&c1_lhs10);
  sf_mex_destroy(&c1_rhs11);
  sf_mex_destroy(&c1_lhs11);
  sf_mex_destroy(&c1_rhs12);
  sf_mex_destroy(&c1_lhs12);
  sf_mex_destroy(&c1_rhs13);
  sf_mex_destroy(&c1_lhs13);
  sf_mex_destroy(&c1_rhs14);
  sf_mex_destroy(&c1_lhs14);
  sf_mex_destroy(&c1_rhs15);
  sf_mex_destroy(&c1_lhs15);
  sf_mex_destroy(&c1_rhs16);
  sf_mex_destroy(&c1_lhs16);
}

static const mxArray *c1_emlrt_marshallOut(const char * c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c1_u)), false);
  return c1_y;
}

static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 7, 0U, 0U, 0U, 0), false);
  return c1_y;
}

static void c1_eml_scalar_eg(SFc1_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_threshold(SFc1_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_b_eml_scalar_eg(SFc1_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_e_emlrt_marshallIn(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i65;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i65, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i65;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static void c1_f_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  c1_sdz1UxT0ktEtEZMZuZvzEsH *c1_y)
{
  emlrtMsgIdentifier c1_thisId;
  static const char * c1_fieldNames[18] = { "thrustLimits", "thrustRateLimits",
    "thrustDir", "L", "nL", "rad", "bladeRad", "mass", "I", "km", "kf", "K",
    "velocityGain", "vGain", "zGain", "thetaGain", "omegaGain", "Ts2" };

  c1_thisId.fParent = c1_parentId;
  sf_mex_check_struct(c1_parentId, c1_u, 18, c1_fieldNames, 0U, NULL);
  c1_thisId.fIdentifier = "thrustLimits";
  c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u,
    "thrustLimits", "thrustLimits", 0)), &c1_thisId, c1_y->thrustLimits);
  c1_thisId.fIdentifier = "thrustRateLimits";
  c1_h_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u,
    "thrustRateLimits", "thrustRateLimits", 0)), &c1_thisId,
                        c1_y->thrustRateLimits);
  c1_thisId.fIdentifier = "thrustDir";
  c1_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u,
    "thrustDir", "thrustDir", 0)), &c1_thisId, c1_y->thrustDir);
  c1_thisId.fIdentifier = "L";
  c1_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u, "L", "L",
    0)), &c1_thisId, c1_y->L);
  c1_thisId.fIdentifier = "nL";
  c1_y->nL = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "nL", "nL", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rad";
  c1_y->rad = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rad", "rad", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "bladeRad";
  c1_y->bladeRad = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c1_u, "bladeRad", "bladeRad", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "mass";
  c1_y->mass = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "mass", "mass", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "I";
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u, "I", "I",
    0)), &c1_thisId, c1_y->I);
  c1_thisId.fIdentifier = "km";
  c1_y->km = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "km", "km", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "kf";
  c1_y->kf = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "kf", "kf", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "K";
  c1_j_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u, "K", "K",
    0)), &c1_thisId, c1_y->K);
  c1_thisId.fIdentifier = "velocityGain";
  c1_y->velocityGain = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c1_u, "velocityGain", "velocityGain", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "vGain";
  c1_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u, "vGain",
    "vGain", 0)), &c1_thisId, c1_y->vGain);
  c1_thisId.fIdentifier = "zGain";
  c1_y->zGain = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "zGain", "zGain", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "thetaGain";
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u,
    "thetaGain", "thetaGain", 0)), &c1_thisId, c1_y->thetaGain);
  c1_thisId.fIdentifier = "omegaGain";
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c1_u,
    "omegaGain", "omegaGain", 0)), &c1_thisId, c1_y->omegaGain);
  c1_thisId.fIdentifier = "Ts2";
  c1_y->Ts2 = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Ts2", "Ts2", 0)), &c1_thisId);
  sf_mex_destroy(&c1_u);
}

static void c1_g_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[8])
{
  real_T c1_dv6[8];
  int32_T c1_i66;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv6, 1, 0, 0U, 1, 0U, 2, 2, 4);
  for (c1_i66 = 0; c1_i66 < 8; c1_i66++) {
    c1_y[c1_i66] = c1_dv6[c1_i66];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_h_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[4])
{
  real_T c1_dv7[4];
  int32_T c1_i67;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv7, 1, 0, 0U, 1, 0U, 2, 1, 4);
  for (c1_i67 = 0; c1_i67 < 4; c1_i67++) {
    c1_y[c1_i67] = c1_dv7[c1_i67];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_i_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[12])
{
  real_T c1_dv8[12];
  int32_T c1_i68;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv8, 1, 0, 0U, 1, 0U, 2, 3, 4);
  for (c1_i68 = 0; c1_i68 < 12; c1_i68++) {
    c1_y[c1_i68] = c1_dv8[c1_i68];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_j_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[16])
{
  real_T c1_dv9[16];
  int32_T c1_i69;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv9, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c1_i69 = 0; c1_i69 < 16; c1_i69++) {
    c1_y[c1_i69] = c1_dv9[c1_i69];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_k_emlrt_marshallIn(SFc1_simulation2InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[2])
{
  real_T c1_dv10[2];
  int32_T c1_i70;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv10, 1, 0, 0U, 1, 0U, 1, 2);
  for (c1_i70 = 0; c1_i70 < 2; c1_i70++) {
    c1_y[c1_i70] = c1_dv10[c1_i70];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_quad;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  c1_sdz1UxT0ktEtEZMZuZvzEsH c1_y;
  SFc1_simulation2InstanceStruct *chartInstance;
  chartInstance = (SFc1_simulation2InstanceStruct *)chartInstanceVoid;
  c1_b_quad = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_quad), &c1_thisId, &c1_y);
  sf_mex_destroy(&c1_b_quad);
  *(c1_sdz1UxT0ktEtEZMZuZvzEsH *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_l_emlrt_marshallIn(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_b_is_active_c1_simulation2, const char_T
  *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_simulation2), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_simulation2);
  return c1_y;
}

static uint8_T c1_m_emlrt_marshallIn(SFc1_simulation2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void init_dsm_address_info(SFc1_simulation2InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc1_simulation2InstanceStruct
  *chartInstance)
{
  chartInstance->c1_uTot = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c1_vDot = (real_T (*)[3])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c1_theta = (real_T (*)[3])ssGetInputPortSignal_wrapper
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

void sf_c1_simulation2_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(894396090U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2717444065U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(711428929U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3332146119U);
}

mxArray* sf_c1_simulation2_get_post_codegen_info(void);
mxArray *sf_c1_simulation2_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("BqBfnT0gSqI7kInyWs9cFC");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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
    mxArray* mxPostCodegenInfo = sf_c1_simulation2_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c1_simulation2_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_simulation2_jit_fallback_info(void)
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

mxArray *sf_c1_simulation2_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c1_simulation2_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c1_simulation2(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"vDot\",},{M[8],M[0],T\"is_active_c1_simulation2\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_simulation2_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_simulation2InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc1_simulation2InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulation2MachineNumber_,
           1,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"uTot");
          _SFD_SET_DATA_PROPS(1,2,0,1,"vDot");
          _SFD_SET_DATA_PROPS(2,1,1,0,"theta");
          _SFD_SET_DATA_PROPS(3,10,0,0,"quad");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,395);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)
            c1_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(3,SF_STRUCT,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)c1_e_sf_marshallIn);
        _SFD_SET_DATA_VALUE_PTR(0U, chartInstance->c1_uTot);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c1_vDot);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c1_theta);
        _SFD_SET_DATA_VALUE_PTR(3U, &chartInstance->c1_quad);
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
  return "TMICya6UtSmg35nhpj4SDE";
}

static void sf_opaque_initialize_c1_simulation2(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_simulation2InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c1_simulation2((SFc1_simulation2InstanceStruct*)
    chartInstanceVar);
  initialize_c1_simulation2((SFc1_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_simulation2(void *chartInstanceVar)
{
  enable_c1_simulation2((SFc1_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_simulation2(void *chartInstanceVar)
{
  disable_c1_simulation2((SFc1_simulation2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_simulation2(void *chartInstanceVar)
{
  sf_gateway_c1_simulation2((SFc1_simulation2InstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c1_simulation2(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c1_simulation2((SFc1_simulation2InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c1_simulation2(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c1_simulation2((SFc1_simulation2InstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c1_simulation2(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_simulation2InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation2_optimization_info();
    }

    finalize_c1_simulation2((SFc1_simulation2InstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_simulation2((SFc1_simulation2InstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_simulation2(SimStruct *S)
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
    initialize_params_c1_simulation2((SFc1_simulation2InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_simulation2(SimStruct *S)
{
  /* Actual parameters from chart:
     quad
   */
  const char_T *rtParamNames[] = { "quad" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation2_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,1,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,1);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1489849645U));
  ssSetChecksum1(S,(1562766917U));
  ssSetChecksum2(S,(3465469377U));
  ssSetChecksum3(S,(1110494450U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_simulation2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_simulation2(SimStruct *S)
{
  SFc1_simulation2InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc1_simulation2InstanceStruct *)utMalloc(sizeof
    (SFc1_simulation2InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_simulation2InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c1_simulation2;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c1_simulation2;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c1_simulation2;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_simulation2;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_simulation2;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c1_simulation2;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c1_simulation2;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c1_simulation2;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_simulation2;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_simulation2;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_simulation2;
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

void c1_simulation2_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_simulation2(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_simulation2(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_simulation2(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_simulation2_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
