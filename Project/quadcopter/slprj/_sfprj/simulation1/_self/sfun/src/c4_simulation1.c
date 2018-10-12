/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation1_sfun.h"
#include "c4_simulation1.h"
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
static const char * c4_debug_family_names[8] = { "d", "quad", "nargin",
  "nargout", "u", "dist", "omega", "omegaDot" };

/* Function Declarations */
static void initialize_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance);
static void initialize_params_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance);
static void enable_c4_simulation1(SFc4_simulation1InstanceStruct *chartInstance);
static void disable_c4_simulation1(SFc4_simulation1InstanceStruct *chartInstance);
static void c4_update_debugger_state_c4_simulation1
  (SFc4_simulation1InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c4_simulation1
  (SFc4_simulation1InstanceStruct *chartInstance);
static void set_sim_state_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_st);
static void finalize_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance);
static void sf_gateway_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance);
static void mdl_start_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance);
static void initSimStructsc4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber, uint32_T c4_instanceNumber);
static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData);
static void c4_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_omegaDot, const char_T *c4_identifier, real_T c4_y[3]);
static void c4_b_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[3]);
static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static real_T c4_c_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_info_helper(const mxArray **c4_info);
static const mxArray *c4_emlrt_marshallOut(const char * c4_b_u);
static const mxArray *c4_b_emlrt_marshallOut(const uint32_T c4_b_u);
static real_T c4_rand(SFc4_simulation1InstanceStruct *chartInstance);
static uint32_T c4_eml_rand_str2id(SFc4_simulation1InstanceStruct *chartInstance);
static uint32_T c4_b_eml_rand_str2id(SFc4_simulation1InstanceStruct
  *chartInstance);
static uint32_T c4_c_eml_rand_str2id(SFc4_simulation1InstanceStruct
  *chartInstance);
static void c4_eml_rand_mcg16807(SFc4_simulation1InstanceStruct *chartInstance,
  uint32_T c4_d_state, uint32_T *c4_e_state, real_T *c4_r);
static void c4_eml_rand_mt19937ar(SFc4_simulation1InstanceStruct *chartInstance,
  uint32_T c4_d_state[625]);
static void c4_twister_state_vector(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T c4_mt[625], real_T c4_seed, uint32_T c4_b_mt[625]);
static void c4_b_eml_rand_mt19937ar(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T c4_d_state[625], uint32_T c4_e_state[625], real_T
  *c4_r);
static void c4_eml_error(SFc4_simulation1InstanceStruct *chartInstance);
static void c4_b_rand(SFc4_simulation1InstanceStruct *chartInstance, real_T
                      c4_r[3]);
static void c4_eml_scalar_eg(SFc4_simulation1InstanceStruct *chartInstance);
static void c4_eml_xgemm(SFc4_simulation1InstanceStruct *chartInstance, real_T
  c4_A[9], real_T c4_B[3], real_T c4_C[3], real_T c4_b_C[3]);
static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static int32_T c4_d_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static void c4_e_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId,
  c4_sdz1UxT0ktEtEZMZuZvzEsH *c4_y);
static void c4_f_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[8]);
static void c4_g_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[4]);
static void c4_h_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[12]);
static void c4_i_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[9]);
static void c4_j_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[16]);
static void c4_k_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[2]);
static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static uint32_T c4_l_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_method, const char_T *c4_identifier);
static uint32_T c4_m_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId);
static uint32_T c4_n_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_d_state, const char_T *c4_identifier);
static uint32_T c4_o_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_p_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_d_state, const char_T *c4_identifier, uint32_T c4_y[625]);
static void c4_q_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, uint32_T c4_y
  [625]);
static void c4_r_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_d_state, const char_T *c4_identifier, uint32_T c4_y[2]);
static void c4_s_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, uint32_T c4_y[2]);
static uint8_T c4_t_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_is_active_c4_simulation1, const char_T
  *c4_identifier);
static uint8_T c4_u_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId);
static real_T c4_b_eml_rand_mcg16807(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T *c4_d_state);
static void c4_b_twister_state_vector(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T c4_mt[625], real_T c4_seed);
static real_T c4_c_eml_rand_mt19937ar(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T c4_d_state[625]);
static void c4_b_eml_xgemm(SFc4_simulation1InstanceStruct *chartInstance, real_T
  c4_A[9], real_T c4_B[3], real_T c4_C[3]);
static void init_dsm_address_info(SFc4_simulation1InstanceStruct *chartInstance);
static void init_simulink_io_address(SFc4_simulation1InstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  chartInstance->c4_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c4_method_not_empty = false;
  chartInstance->c4_state_not_empty = false;
  chartInstance->c4_b_state_not_empty = false;
  chartInstance->c4_c_state_not_empty = false;
  chartInstance->c4_is_active_c4_simulation1 = 0U;
}

static void initialize_params_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  const mxArray *c4_m0 = NULL;
  const mxArray *c4_mxField;
  c4_sdz1UxT0ktEtEZMZuZvzEsH c4_r0;
  c4_m0 = sf_mex_get_sfun_param(chartInstance->S, 0, 1);
  c4_mxField = sf_mex_getfield(c4_m0, "thrustLimits", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.thrustLimits, 1, 0,
                      0U, 1, 0U, 2, 2, 4);
  c4_mxField = sf_mex_getfield(c4_m0, "thrustRateLimits", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.thrustRateLimits, 1,
                      0, 0U, 1, 0U, 2, 1, 4);
  c4_mxField = sf_mex_getfield(c4_m0, "thrustDir", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.thrustDir, 1, 0, 0U,
                      1, 0U, 2, 3, 4);
  c4_mxField = sf_mex_getfield(c4_m0, "L", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.L, 1, 0, 0U, 1, 0U,
                      2, 3, 4);
  c4_mxField = sf_mex_getfield(c4_m0, "nL", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.nL, 1, 0, 0U, 0, 0U,
                      0);
  c4_mxField = sf_mex_getfield(c4_m0, "rad", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.rad, 1, 0, 0U, 0,
                      0U, 0);
  c4_mxField = sf_mex_getfield(c4_m0, "bladeRad", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.bladeRad, 1, 0, 0U,
                      0, 0U, 0);
  c4_mxField = sf_mex_getfield(c4_m0, "mass", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.mass, 1, 0, 0U, 0,
                      0U, 0);
  c4_mxField = sf_mex_getfield(c4_m0, "I", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.I, 1, 0, 0U, 1, 0U,
                      2, 3, 3);
  c4_mxField = sf_mex_getfield(c4_m0, "km", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.km, 1, 0, 0U, 0, 0U,
                      0);
  c4_mxField = sf_mex_getfield(c4_m0, "kf", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.kf, 1, 0, 0U, 0, 0U,
                      0);
  c4_mxField = sf_mex_getfield(c4_m0, "K", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.K, 1, 0, 0U, 1, 0U,
                      2, 4, 4);
  c4_mxField = sf_mex_getfield(c4_m0, "velocityGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.velocityGain, 1, 0,
                      0U, 0, 0U, 0);
  c4_mxField = sf_mex_getfield(c4_m0, "vGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.vGain, 1, 0, 0U, 1,
                      0U, 1, 2);
  c4_mxField = sf_mex_getfield(c4_m0, "zGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.zGain, 1, 0, 0U, 0,
                      0U, 0);
  c4_mxField = sf_mex_getfield(c4_m0, "thetaGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.thetaGain, 1, 0, 0U,
                      1, 0U, 1, 3);
  c4_mxField = sf_mex_getfield(c4_m0, "omegaGain", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), c4_r0.omegaGain, 1, 0, 0U,
                      1, 0U, 1, 3);
  c4_mxField = sf_mex_getfield(c4_m0, "Ts2", "quad", 0);
  sf_mex_import_named("quad", sf_mex_dup(c4_mxField), &c4_r0.Ts2, 1, 0, 0U, 0,
                      0U, 0);
  sf_mex_destroy(&c4_m0);
  chartInstance->c4_quad = c4_r0;
}

static void enable_c4_simulation1(SFc4_simulation1InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c4_simulation1(SFc4_simulation1InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c4_update_debugger_state_c4_simulation1
  (SFc4_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c4_simulation1
  (SFc4_simulation1InstanceStruct *chartInstance)
{
  const mxArray *c4_st;
  const mxArray *c4_y = NULL;
  int32_T c4_i0;
  real_T c4_b_u[3];
  const mxArray *c4_b_y = NULL;
  uint32_T c4_hoistedGlobal;
  uint32_T c4_c_u;
  const mxArray *c4_c_y = NULL;
  uint32_T c4_b_hoistedGlobal;
  uint32_T c4_d_u;
  const mxArray *c4_d_y = NULL;
  int32_T c4_i1;
  uint32_T c4_e_u[625];
  const mxArray *c4_e_y = NULL;
  int32_T c4_i2;
  uint32_T c4_f_u[2];
  const mxArray *c4_f_y = NULL;
  uint8_T c4_c_hoistedGlobal;
  uint8_T c4_g_u;
  const mxArray *c4_g_y = NULL;
  c4_st = NULL;
  c4_st = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_createcellmatrix(6, 1), false);
  for (c4_i0 = 0; c4_i0 < 3; c4_i0++) {
    c4_b_u[c4_i0] = (*chartInstance->c4_omegaDot)[c4_i0];
  }

  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_b_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c4_y, 0, c4_b_y);
  c4_hoistedGlobal = chartInstance->c4_method;
  c4_c_u = c4_hoistedGlobal;
  c4_c_y = NULL;
  if (!chartInstance->c4_method_not_empty) {
    sf_mex_assign(&c4_c_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c4_c_y, sf_mex_create("y", &c4_c_u, 7, 0U, 0U, 0U, 0), false);
  }

  sf_mex_setcell(c4_y, 1, c4_c_y);
  c4_b_hoistedGlobal = chartInstance->c4_state;
  c4_d_u = c4_b_hoistedGlobal;
  c4_d_y = NULL;
  if (!chartInstance->c4_state_not_empty) {
    sf_mex_assign(&c4_d_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c4_d_y, sf_mex_create("y", &c4_d_u, 7, 0U, 0U, 0U, 0), false);
  }

  sf_mex_setcell(c4_y, 2, c4_d_y);
  for (c4_i1 = 0; c4_i1 < 625; c4_i1++) {
    c4_e_u[c4_i1] = chartInstance->c4_c_state[c4_i1];
  }

  c4_e_y = NULL;
  if (!chartInstance->c4_c_state_not_empty) {
    sf_mex_assign(&c4_e_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c4_e_y, sf_mex_create("y", c4_e_u, 7, 0U, 1U, 0U, 1, 625),
                  false);
  }

  sf_mex_setcell(c4_y, 3, c4_e_y);
  for (c4_i2 = 0; c4_i2 < 2; c4_i2++) {
    c4_f_u[c4_i2] = chartInstance->c4_b_state[c4_i2];
  }

  c4_f_y = NULL;
  if (!chartInstance->c4_b_state_not_empty) {
    sf_mex_assign(&c4_f_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c4_f_y, sf_mex_create("y", c4_f_u, 7, 0U, 1U, 0U, 1, 2),
                  false);
  }

  sf_mex_setcell(c4_y, 4, c4_f_y);
  c4_c_hoistedGlobal = chartInstance->c4_is_active_c4_simulation1;
  c4_g_u = c4_c_hoistedGlobal;
  c4_g_y = NULL;
  sf_mex_assign(&c4_g_y, sf_mex_create("y", &c4_g_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c4_y, 5, c4_g_y);
  sf_mex_assign(&c4_st, c4_y, false);
  return c4_st;
}

static void set_sim_state_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_st)
{
  const mxArray *c4_b_u;
  real_T c4_dv0[3];
  int32_T c4_i3;
  uint32_T c4_uv0[625];
  int32_T c4_i4;
  uint32_T c4_uv1[2];
  int32_T c4_i5;
  chartInstance->c4_doneDoubleBufferReInit = true;
  c4_b_u = sf_mex_dup(c4_st);
  c4_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_b_u, 0)),
                      "omegaDot", c4_dv0);
  for (c4_i3 = 0; c4_i3 < 3; c4_i3++) {
    (*chartInstance->c4_omegaDot)[c4_i3] = c4_dv0[c4_i3];
  }

  chartInstance->c4_method = c4_l_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c4_b_u, 1)), "method");
  chartInstance->c4_state = c4_n_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c4_b_u, 2)), "state");
  c4_p_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_b_u, 3)),
                        "state", c4_uv0);
  for (c4_i4 = 0; c4_i4 < 625; c4_i4++) {
    chartInstance->c4_c_state[c4_i4] = c4_uv0[c4_i4];
  }

  c4_r_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_b_u, 4)),
                        "state", c4_uv1);
  for (c4_i5 = 0; c4_i5 < 2; c4_i5++) {
    chartInstance->c4_b_state[c4_i5] = c4_uv1[c4_i5];
  }

  chartInstance->c4_is_active_c4_simulation1 = c4_t_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c4_b_u, 5)),
     "is_active_c4_simulation1");
  sf_mex_destroy(&c4_b_u);
  c4_update_debugger_state_c4_simulation1(chartInstance);
  sf_mex_destroy(&c4_st);
}

static void finalize_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  int32_T c4_i6;
  int32_T c4_i7;
  real_T c4_hoistedGlobal;
  int32_T c4_i8;
  real_T c4_b_u[3];
  real_T c4_b_dist;
  int32_T c4_i9;
  real_T c4_b_omega[3];
  uint32_T c4_debug_family_var_map[8];
  real_T c4_d[3];
  c4_sdz1UxT0ktEtEZMZuZvzEsH c4_b_quad;
  real_T c4_nargin = 4.0;
  real_T c4_nargout = 1.0;
  real_T c4_b_omegaDot[3];
  static c4_sdz1UxT0ktEtEZMZuZvzEsH c4_r1 = { { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
      0.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 1.0 }, { 0.0, 0.2, 0.0, -0.2, 0.0, 0.0, 0.0, -0.2, 0.0,
      0.2, 0.0, 0.0 }, 0.2, 0.040000000000000008, 0.080000000000000016, 8.0, {
      10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 15.0 }, 11.0, 28.0, { 28.0, 0.0,
      -5.6000000000000005, 11.0, 28.0, 5.6000000000000005, 0.0, -11.0, 28.0, 0.0,
      5.6000000000000005, 11.0, 28.0, -5.6000000000000005, 0.0, -11.0 }, 10.0, {
      0.1, 0.1 }, 1.0, { 1.0, 1.0, 1.0 }, { 0.2, 0.2, 0.2 }, 0.3 };

  real_T c4_b;
  int32_T c4_i10;
  static real_T c4_a[3] = { 0.0, 0.0, 40.0 };

  real_T c4_y[3];
  real_T c4_b_b[3];
  int32_T c4_i11;
  int32_T c4_i12;
  int32_T c4_i13;
  int32_T c4_i14;
  int32_T c4_i15;
  int32_T c4_i16;
  int32_T c4_i17;
  static real_T c4_b_a[9] = { 10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 15.0 };

  real_T c4_c_a[9];
  int32_T c4_i18;
  real_T c4_c_b[3];
  int32_T c4_i19;
  real_T c4_c1;
  real_T c4_c2;
  real_T c4_c3;
  int32_T c4_i20;
  int32_T c4_i21;
  int32_T c4_i22;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c4_sfEvent);
  for (c4_i6 = 0; c4_i6 < 3; c4_i6++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c4_u)[c4_i6], 0U);
  }

  _SFD_DATA_RANGE_CHECK(*chartInstance->c4_dist, 1U);
  for (c4_i7 = 0; c4_i7 < 3; c4_i7++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c4_omega)[c4_i7], 2U);
  }

  chartInstance->c4_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c4_sfEvent);
  c4_hoistedGlobal = *chartInstance->c4_dist;
  for (c4_i8 = 0; c4_i8 < 3; c4_i8++) {
    c4_b_u[c4_i8] = (*chartInstance->c4_u)[c4_i8];
  }

  c4_b_dist = c4_hoistedGlobal;
  for (c4_i9 = 0; c4_i9 < 3; c4_i9++) {
    c4_b_omega[c4_i9] = (*chartInstance->c4_omega)[c4_i9];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c4_debug_family_names,
    c4_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_d, 0U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_b_quad, 1U, c4_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargin, 2U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargout, 3U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_b_u, 4U, c4_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_b_dist, 5U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_b_omega, 6U, c4_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_b_omegaDot, 7U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  c4_b_quad = c4_r1;
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 3);
  if (CV_EML_IF(0, 1, 0, c4_b_dist != 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 4);
    c4_b = c4_rand(chartInstance) - 0.5;
    for (c4_i10 = 0; c4_i10 < 3; c4_i10++) {
      c4_y[c4_i10] = c4_a[c4_i10] * c4_b;
    }

    c4_b_rand(chartInstance, c4_b_b);
    for (c4_i11 = 0; c4_i11 < 3; c4_i11++) {
      c4_b_b[c4_i11] -= 0.5;
    }

    for (c4_i12 = 0; c4_i12 < 3; c4_i12++) {
      c4_b_b[c4_i12] *= 4.0;
    }

    for (c4_i13 = 0; c4_i13 < 3; c4_i13++) {
      c4_d[c4_i13] = (1.0 + c4_y[c4_i13]) + c4_b_b[c4_i13];
    }
  } else {
    _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 6);
    for (c4_i14 = 0; c4_i14 < 3; c4_i14++) {
      c4_d[c4_i14] = 0.0;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 9);
  for (c4_i15 = 0; c4_i15 < 3; c4_i15++) {
    c4_b_b[c4_i15] = c4_b_omega[c4_i15];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i16 = 0; c4_i16 < 3; c4_i16++) {
    c4_y[c4_i16] = 0.0;
  }

  for (c4_i17 = 0; c4_i17 < 9; c4_i17++) {
    c4_c_a[c4_i17] = c4_b_a[c4_i17];
  }

  for (c4_i18 = 0; c4_i18 < 3; c4_i18++) {
    c4_c_b[c4_i18] = c4_b_b[c4_i18];
  }

  c4_b_eml_xgemm(chartInstance, c4_c_a, c4_c_b, c4_y);
  for (c4_i19 = 0; c4_i19 < 3; c4_i19++) {
    c4_b_b[c4_i19] = c4_b_omega[c4_i19];
  }

  c4_c1 = c4_b_b[1] * c4_y[2] - c4_b_b[2] * c4_y[1];
  c4_c2 = c4_b_b[2] * c4_y[0] - c4_b_b[0] * c4_y[2];
  c4_c3 = c4_b_b[0] * c4_y[1] - c4_b_b[1] * c4_y[0];
  c4_b_b[0] = c4_c1;
  c4_b_b[1] = c4_c2;
  c4_b_b[2] = c4_c3;
  for (c4_i20 = 0; c4_i20 < 3; c4_i20++) {
    c4_b_omegaDot[c4_i20] = (-c4_b_b[c4_i20] + c4_b_u[c4_i20]) + c4_d[c4_i20];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, -9);
  _SFD_SYMBOL_SCOPE_POP();
  for (c4_i21 = 0; c4_i21 < 3; c4_i21++) {
    (*chartInstance->c4_omegaDot)[c4_i21] = c4_b_omegaDot[c4_i21];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c4_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulation1MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c4_i22 = 0; c4_i22 < 3; c4_i22++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c4_omegaDot)[c4_i22], 4U);
  }
}

static void mdl_start_c4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void initSimStructsc4_simulation1(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber, uint32_T c4_instanceNumber)
{
  (void)c4_machineNumber;
  (void)c4_chartNumber;
  (void)c4_instanceNumber;
}

static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i23;
  real_T c4_b_inData[3];
  int32_T c4_i24;
  real_T c4_b_u[3];
  const mxArray *c4_y = NULL;
  SFc4_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc4_simulation1InstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i23 = 0; c4_i23 < 3; c4_i23++) {
    c4_b_inData[c4_i23] = (*(real_T (*)[3])c4_inData)[c4_i23];
  }

  for (c4_i24 = 0; c4_i24 < 3; c4_i24++) {
    c4_b_u[c4_i24] = c4_b_inData[c4_i24];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_b_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_omegaDot, const char_T *c4_identifier, real_T c4_y[3])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_omegaDot), &c4_thisId,
                        c4_y);
  sf_mex_destroy(&c4_b_omegaDot);
}

static void c4_b_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[3])
{
  real_T c4_dv1[3];
  int32_T c4_i25;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c4_i25 = 0; c4_i25 < 3; c4_i25++) {
    c4_y[c4_i25] = c4_dv1[c4_i25];
  }

  sf_mex_destroy(&c4_b_u);
}

static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_b_omegaDot;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[3];
  int32_T c4_i26;
  SFc4_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc4_simulation1InstanceStruct *)chartInstanceVoid;
  c4_b_omegaDot = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_omegaDot), &c4_thisId,
                        c4_y);
  sf_mex_destroy(&c4_b_omegaDot);
  for (c4_i26 = 0; c4_i26 < 3; c4_i26++) {
    (*(real_T (*)[3])c4_outData)[c4_i26] = c4_y[c4_i26];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  real_T c4_b_u;
  const mxArray *c4_y = NULL;
  SFc4_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc4_simulation1InstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_b_u = *(real_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static real_T c4_c_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId)
{
  real_T c4_y;
  real_T c4_d0;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), &c4_d0, 1, 0, 0U, 0, 0U, 0);
  c4_y = c4_d0;
  sf_mex_destroy(&c4_b_u);
  return c4_y;
}

static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_nargout;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y;
  SFc4_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc4_simulation1InstanceStruct *)chartInstanceVoid;
  c4_nargout = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_nargout), &c4_thisId);
  sf_mex_destroy(&c4_nargout);
  *(real_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  c4_sdz1UxT0ktEtEZMZuZvzEsH c4_b_u;
  const mxArray *c4_y = NULL;
  int32_T c4_i27;
  real_T c4_c_u[8];
  const mxArray *c4_b_y = NULL;
  int32_T c4_i28;
  real_T c4_d_u[4];
  const mxArray *c4_c_y = NULL;
  int32_T c4_i29;
  real_T c4_e_u[12];
  const mxArray *c4_d_y = NULL;
  int32_T c4_i30;
  real_T c4_f_u[12];
  const mxArray *c4_e_y = NULL;
  real_T c4_g_u;
  const mxArray *c4_f_y = NULL;
  real_T c4_h_u;
  const mxArray *c4_g_y = NULL;
  real_T c4_i_u;
  const mxArray *c4_h_y = NULL;
  real_T c4_j_u;
  const mxArray *c4_i_y = NULL;
  int32_T c4_i31;
  real_T c4_k_u[9];
  const mxArray *c4_j_y = NULL;
  real_T c4_l_u;
  const mxArray *c4_k_y = NULL;
  real_T c4_m_u;
  const mxArray *c4_l_y = NULL;
  int32_T c4_i32;
  real_T c4_n_u[16];
  const mxArray *c4_m_y = NULL;
  real_T c4_o_u;
  const mxArray *c4_n_y = NULL;
  int32_T c4_i33;
  real_T c4_p_u[2];
  const mxArray *c4_o_y = NULL;
  real_T c4_q_u;
  const mxArray *c4_p_y = NULL;
  int32_T c4_i34;
  real_T c4_r_u[3];
  const mxArray *c4_q_y = NULL;
  int32_T c4_i35;
  real_T c4_s_u[3];
  const mxArray *c4_r_y = NULL;
  real_T c4_t_u;
  const mxArray *c4_s_y = NULL;
  SFc4_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc4_simulation1InstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_b_u = *(c4_sdz1UxT0ktEtEZMZuZvzEsH *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  for (c4_i27 = 0; c4_i27 < 8; c4_i27++) {
    c4_c_u[c4_i27] = c4_b_u.thrustLimits[c4_i27];
  }

  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_c_u, 0, 0U, 1U, 0U, 2, 2, 4),
                false);
  sf_mex_addfield(c4_y, c4_b_y, "thrustLimits", "thrustLimits", 0);
  for (c4_i28 = 0; c4_i28 < 4; c4_i28++) {
    c4_d_u[c4_i28] = c4_b_u.thrustRateLimits[c4_i28];
  }

  c4_c_y = NULL;
  sf_mex_assign(&c4_c_y, sf_mex_create("y", c4_d_u, 0, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_addfield(c4_y, c4_c_y, "thrustRateLimits", "thrustRateLimits", 0);
  for (c4_i29 = 0; c4_i29 < 12; c4_i29++) {
    c4_e_u[c4_i29] = c4_b_u.thrustDir[c4_i29];
  }

  c4_d_y = NULL;
  sf_mex_assign(&c4_d_y, sf_mex_create("y", c4_e_u, 0, 0U, 1U, 0U, 2, 3, 4),
                false);
  sf_mex_addfield(c4_y, c4_d_y, "thrustDir", "thrustDir", 0);
  for (c4_i30 = 0; c4_i30 < 12; c4_i30++) {
    c4_f_u[c4_i30] = c4_b_u.L[c4_i30];
  }

  c4_e_y = NULL;
  sf_mex_assign(&c4_e_y, sf_mex_create("y", c4_f_u, 0, 0U, 1U, 0U, 2, 3, 4),
                false);
  sf_mex_addfield(c4_y, c4_e_y, "L", "L", 0);
  c4_g_u = c4_b_u.nL;
  c4_f_y = NULL;
  sf_mex_assign(&c4_f_y, sf_mex_create("y", &c4_g_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_f_y, "nL", "nL", 0);
  c4_h_u = c4_b_u.rad;
  c4_g_y = NULL;
  sf_mex_assign(&c4_g_y, sf_mex_create("y", &c4_h_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_g_y, "rad", "rad", 0);
  c4_i_u = c4_b_u.bladeRad;
  c4_h_y = NULL;
  sf_mex_assign(&c4_h_y, sf_mex_create("y", &c4_i_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_h_y, "bladeRad", "bladeRad", 0);
  c4_j_u = c4_b_u.mass;
  c4_i_y = NULL;
  sf_mex_assign(&c4_i_y, sf_mex_create("y", &c4_j_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_i_y, "mass", "mass", 0);
  for (c4_i31 = 0; c4_i31 < 9; c4_i31++) {
    c4_k_u[c4_i31] = c4_b_u.I[c4_i31];
  }

  c4_j_y = NULL;
  sf_mex_assign(&c4_j_y, sf_mex_create("y", c4_k_u, 0, 0U, 1U, 0U, 2, 3, 3),
                false);
  sf_mex_addfield(c4_y, c4_j_y, "I", "I", 0);
  c4_l_u = c4_b_u.km;
  c4_k_y = NULL;
  sf_mex_assign(&c4_k_y, sf_mex_create("y", &c4_l_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_k_y, "km", "km", 0);
  c4_m_u = c4_b_u.kf;
  c4_l_y = NULL;
  sf_mex_assign(&c4_l_y, sf_mex_create("y", &c4_m_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_l_y, "kf", "kf", 0);
  for (c4_i32 = 0; c4_i32 < 16; c4_i32++) {
    c4_n_u[c4_i32] = c4_b_u.K[c4_i32];
  }

  c4_m_y = NULL;
  sf_mex_assign(&c4_m_y, sf_mex_create("y", c4_n_u, 0, 0U, 1U, 0U, 2, 4, 4),
                false);
  sf_mex_addfield(c4_y, c4_m_y, "K", "K", 0);
  c4_o_u = c4_b_u.velocityGain;
  c4_n_y = NULL;
  sf_mex_assign(&c4_n_y, sf_mex_create("y", &c4_o_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_n_y, "velocityGain", "velocityGain", 0);
  for (c4_i33 = 0; c4_i33 < 2; c4_i33++) {
    c4_p_u[c4_i33] = c4_b_u.vGain[c4_i33];
  }

  c4_o_y = NULL;
  sf_mex_assign(&c4_o_y, sf_mex_create("y", c4_p_u, 0, 0U, 1U, 0U, 1, 2), false);
  sf_mex_addfield(c4_y, c4_o_y, "vGain", "vGain", 0);
  c4_q_u = c4_b_u.zGain;
  c4_p_y = NULL;
  sf_mex_assign(&c4_p_y, sf_mex_create("y", &c4_q_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_p_y, "zGain", "zGain", 0);
  for (c4_i34 = 0; c4_i34 < 3; c4_i34++) {
    c4_r_u[c4_i34] = c4_b_u.thetaGain[c4_i34];
  }

  c4_q_y = NULL;
  sf_mex_assign(&c4_q_y, sf_mex_create("y", c4_r_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_addfield(c4_y, c4_q_y, "thetaGain", "thetaGain", 0);
  for (c4_i35 = 0; c4_i35 < 3; c4_i35++) {
    c4_s_u[c4_i35] = c4_b_u.omegaGain[c4_i35];
  }

  c4_r_y = NULL;
  sf_mex_assign(&c4_r_y, sf_mex_create("y", c4_s_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_addfield(c4_y, c4_r_y, "omegaGain", "omegaGain", 0);
  c4_t_u = c4_b_u.Ts2;
  c4_s_y = NULL;
  sf_mex_assign(&c4_s_y, sf_mex_create("y", &c4_t_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c4_y, c4_s_y, "Ts2", "Ts2", 0);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

const mxArray *sf_c4_simulation1_get_eml_resolved_functions_info(void)
{
  const mxArray *c4_nameCaptureInfo = NULL;
  c4_nameCaptureInfo = NULL;
  sf_mex_assign(&c4_nameCaptureInfo, sf_mex_createstruct("structure", 2, 39, 1),
                false);
  c4_info_helper(&c4_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c4_nameCaptureInfo);
  return c4_nameCaptureInfo;
}

static void c4_info_helper(const mxArray **c4_info)
{
  const mxArray *c4_rhs0 = NULL;
  const mxArray *c4_lhs0 = NULL;
  const mxArray *c4_rhs1 = NULL;
  const mxArray *c4_lhs1 = NULL;
  const mxArray *c4_rhs2 = NULL;
  const mxArray *c4_lhs2 = NULL;
  const mxArray *c4_rhs3 = NULL;
  const mxArray *c4_lhs3 = NULL;
  const mxArray *c4_rhs4 = NULL;
  const mxArray *c4_lhs4 = NULL;
  const mxArray *c4_rhs5 = NULL;
  const mxArray *c4_lhs5 = NULL;
  const mxArray *c4_rhs6 = NULL;
  const mxArray *c4_lhs6 = NULL;
  const mxArray *c4_rhs7 = NULL;
  const mxArray *c4_lhs7 = NULL;
  const mxArray *c4_rhs8 = NULL;
  const mxArray *c4_lhs8 = NULL;
  const mxArray *c4_rhs9 = NULL;
  const mxArray *c4_lhs9 = NULL;
  const mxArray *c4_rhs10 = NULL;
  const mxArray *c4_lhs10 = NULL;
  const mxArray *c4_rhs11 = NULL;
  const mxArray *c4_lhs11 = NULL;
  const mxArray *c4_rhs12 = NULL;
  const mxArray *c4_lhs12 = NULL;
  const mxArray *c4_rhs13 = NULL;
  const mxArray *c4_lhs13 = NULL;
  const mxArray *c4_rhs14 = NULL;
  const mxArray *c4_lhs14 = NULL;
  const mxArray *c4_rhs15 = NULL;
  const mxArray *c4_lhs15 = NULL;
  const mxArray *c4_rhs16 = NULL;
  const mxArray *c4_lhs16 = NULL;
  const mxArray *c4_rhs17 = NULL;
  const mxArray *c4_lhs17 = NULL;
  const mxArray *c4_rhs18 = NULL;
  const mxArray *c4_lhs18 = NULL;
  const mxArray *c4_rhs19 = NULL;
  const mxArray *c4_lhs19 = NULL;
  const mxArray *c4_rhs20 = NULL;
  const mxArray *c4_lhs20 = NULL;
  const mxArray *c4_rhs21 = NULL;
  const mxArray *c4_lhs21 = NULL;
  const mxArray *c4_rhs22 = NULL;
  const mxArray *c4_lhs22 = NULL;
  const mxArray *c4_rhs23 = NULL;
  const mxArray *c4_lhs23 = NULL;
  const mxArray *c4_rhs24 = NULL;
  const mxArray *c4_lhs24 = NULL;
  const mxArray *c4_rhs25 = NULL;
  const mxArray *c4_lhs25 = NULL;
  const mxArray *c4_rhs26 = NULL;
  const mxArray *c4_lhs26 = NULL;
  const mxArray *c4_rhs27 = NULL;
  const mxArray *c4_lhs27 = NULL;
  const mxArray *c4_rhs28 = NULL;
  const mxArray *c4_lhs28 = NULL;
  const mxArray *c4_rhs29 = NULL;
  const mxArray *c4_lhs29 = NULL;
  const mxArray *c4_rhs30 = NULL;
  const mxArray *c4_lhs30 = NULL;
  const mxArray *c4_rhs31 = NULL;
  const mxArray *c4_lhs31 = NULL;
  const mxArray *c4_rhs32 = NULL;
  const mxArray *c4_lhs32 = NULL;
  const mxArray *c4_rhs33 = NULL;
  const mxArray *c4_lhs33 = NULL;
  const mxArray *c4_rhs34 = NULL;
  const mxArray *c4_lhs34 = NULL;
  const mxArray *c4_rhs35 = NULL;
  const mxArray *c4_lhs35 = NULL;
  const mxArray *c4_rhs36 = NULL;
  const mxArray *c4_lhs36 = NULL;
  const mxArray *c4_rhs37 = NULL;
  const mxArray *c4_lhs37 = NULL;
  const mxArray *c4_rhs38 = NULL;
  const mxArray *c4_lhs38 = NULL;
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("rand"), "name", "name", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/rand.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383877290U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c4_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/rand.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_is_rand_extrinsic"),
                  "name", "name", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_is_rand_extrinsic.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1368183032U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c4_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/rand.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand"), "name", "name", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c4_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_str2id"), "name",
                  "name", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_str2id.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1313347822U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c4_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_str2id.m"),
                  "context", "context", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c4_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_mcg16807_stateful"),
                  "name", "name", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807_stateful.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1366162244U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c4_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807_stateful.m"),
                  "context", "context", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_mcg16807"), "name",
                  "name", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c4_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807_stateful.m"),
                  "context", "context", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_mcg16807"), "name",
                  "name", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("uint32"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c4_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_shr3cong_stateful"),
                  "name", "name", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong_stateful.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1366162244U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c4_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong_stateful.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_shr3cong"), "name",
                  "name", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c4_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong_stateful.m"),
                  "context", "context", 10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_shr3cong"), "name",
                  "name", 10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("uint32"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1313347820U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c4_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_mt19937ar_stateful"),
                  "name", "name", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar_stateful.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1366162244U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c4_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar_stateful.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_mt19937ar"), "name",
                  "name", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1406813148U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c4_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar_stateful.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_rand_mt19937ar"), "name",
                  "name", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("uint32"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1406813148U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c4_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!genrandu"),
                  "context", "context", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eps"), "name", "name", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c4_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_eps"), "name", "name", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c4_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c4_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!is_valid_state"),
                  "context", "context", 17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isequal"), "name", "name", 17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286818758U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c4_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286818786U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c4_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m!isequal_scalar"),
                  "context", "context", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isnan"), "name", "name", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363713858U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c4_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c4_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!is_valid_state"),
                  "context", "context", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c4_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!is_valid_state"),
                  "context", "context", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c4_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c4_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar.m!genrandu"),
                  "context", "context", 24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_error"), "name", "name",
                  24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343830358U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c4_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c4_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c4_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c4_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c4_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c4_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375980690U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c4_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c4_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c4_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 33);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 33);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c4_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 34);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 34);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c4_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 35);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 35);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c4_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 36);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 36);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c4_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 37);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 37);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c4_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 38);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("cross"), "name", "name", 38);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/cross.m"), "resolved",
                  "resolved", 38);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286818842U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c4_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs38), "lhs", "lhs",
                  38);
  sf_mex_destroy(&c4_rhs0);
  sf_mex_destroy(&c4_lhs0);
  sf_mex_destroy(&c4_rhs1);
  sf_mex_destroy(&c4_lhs1);
  sf_mex_destroy(&c4_rhs2);
  sf_mex_destroy(&c4_lhs2);
  sf_mex_destroy(&c4_rhs3);
  sf_mex_destroy(&c4_lhs3);
  sf_mex_destroy(&c4_rhs4);
  sf_mex_destroy(&c4_lhs4);
  sf_mex_destroy(&c4_rhs5);
  sf_mex_destroy(&c4_lhs5);
  sf_mex_destroy(&c4_rhs6);
  sf_mex_destroy(&c4_lhs6);
  sf_mex_destroy(&c4_rhs7);
  sf_mex_destroy(&c4_lhs7);
  sf_mex_destroy(&c4_rhs8);
  sf_mex_destroy(&c4_lhs8);
  sf_mex_destroy(&c4_rhs9);
  sf_mex_destroy(&c4_lhs9);
  sf_mex_destroy(&c4_rhs10);
  sf_mex_destroy(&c4_lhs10);
  sf_mex_destroy(&c4_rhs11);
  sf_mex_destroy(&c4_lhs11);
  sf_mex_destroy(&c4_rhs12);
  sf_mex_destroy(&c4_lhs12);
  sf_mex_destroy(&c4_rhs13);
  sf_mex_destroy(&c4_lhs13);
  sf_mex_destroy(&c4_rhs14);
  sf_mex_destroy(&c4_lhs14);
  sf_mex_destroy(&c4_rhs15);
  sf_mex_destroy(&c4_lhs15);
  sf_mex_destroy(&c4_rhs16);
  sf_mex_destroy(&c4_lhs16);
  sf_mex_destroy(&c4_rhs17);
  sf_mex_destroy(&c4_lhs17);
  sf_mex_destroy(&c4_rhs18);
  sf_mex_destroy(&c4_lhs18);
  sf_mex_destroy(&c4_rhs19);
  sf_mex_destroy(&c4_lhs19);
  sf_mex_destroy(&c4_rhs20);
  sf_mex_destroy(&c4_lhs20);
  sf_mex_destroy(&c4_rhs21);
  sf_mex_destroy(&c4_lhs21);
  sf_mex_destroy(&c4_rhs22);
  sf_mex_destroy(&c4_lhs22);
  sf_mex_destroy(&c4_rhs23);
  sf_mex_destroy(&c4_lhs23);
  sf_mex_destroy(&c4_rhs24);
  sf_mex_destroy(&c4_lhs24);
  sf_mex_destroy(&c4_rhs25);
  sf_mex_destroy(&c4_lhs25);
  sf_mex_destroy(&c4_rhs26);
  sf_mex_destroy(&c4_lhs26);
  sf_mex_destroy(&c4_rhs27);
  sf_mex_destroy(&c4_lhs27);
  sf_mex_destroy(&c4_rhs28);
  sf_mex_destroy(&c4_lhs28);
  sf_mex_destroy(&c4_rhs29);
  sf_mex_destroy(&c4_lhs29);
  sf_mex_destroy(&c4_rhs30);
  sf_mex_destroy(&c4_lhs30);
  sf_mex_destroy(&c4_rhs31);
  sf_mex_destroy(&c4_lhs31);
  sf_mex_destroy(&c4_rhs32);
  sf_mex_destroy(&c4_lhs32);
  sf_mex_destroy(&c4_rhs33);
  sf_mex_destroy(&c4_lhs33);
  sf_mex_destroy(&c4_rhs34);
  sf_mex_destroy(&c4_lhs34);
  sf_mex_destroy(&c4_rhs35);
  sf_mex_destroy(&c4_lhs35);
  sf_mex_destroy(&c4_rhs36);
  sf_mex_destroy(&c4_lhs36);
  sf_mex_destroy(&c4_rhs37);
  sf_mex_destroy(&c4_lhs37);
  sf_mex_destroy(&c4_rhs38);
  sf_mex_destroy(&c4_lhs38);
}

static const mxArray *c4_emlrt_marshallOut(const char * c4_b_u)
{
  const mxArray *c4_y = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_b_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c4_b_u)), false);
  return c4_y;
}

static const mxArray *c4_b_emlrt_marshallOut(const uint32_T c4_b_u)
{
  const mxArray *c4_y = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_b_u, 7, 0U, 0U, 0U, 0), false);
  return c4_y;
}

static real_T c4_rand(SFc4_simulation1InstanceStruct *chartInstance)
{
  real_T c4_r;
  uint32_T c4_V4;
  uint32_T c4_V5NORMAL;
  uint32_T c4_TWISTER;
  uint32_T c4_d_state;
  real_T c4_d1;
  int32_T c4_i36;
  uint32_T c4_icng;
  uint32_T c4_jsr;
  uint32_T c4_u0;
  uint32_T c4_u1;
  uint32_T c4_ui;
  uint32_T c4_b_ui;
  real_T c4_b_r;
  real_T c4_d2;
  uint32_T c4_uv2[625];
  int32_T c4_i37;
  real_T c4_d3;
  c4_V4 = c4_eml_rand_str2id(chartInstance);
  c4_V5NORMAL = c4_b_eml_rand_str2id(chartInstance);
  c4_TWISTER = c4_c_eml_rand_str2id(chartInstance);
  if (!chartInstance->c4_method_not_empty) {
    chartInstance->c4_method = c4_TWISTER;
    chartInstance->c4_method_not_empty = true;
  }

  if (chartInstance->c4_method == c4_V4) {
    if (!chartInstance->c4_state_not_empty) {
      chartInstance->c4_state = 1144108930U;
      chartInstance->c4_state_not_empty = true;
    }

    c4_d_state = chartInstance->c4_state;
    c4_d1 = c4_b_eml_rand_mcg16807(chartInstance, &c4_d_state);
    chartInstance->c4_state = c4_d_state;
    c4_r = c4_d1;
  } else if (chartInstance->c4_method == c4_V5NORMAL) {
    if (!chartInstance->c4_b_state_not_empty) {
      for (c4_i36 = 0; c4_i36 < 2; c4_i36++) {
        chartInstance->c4_b_state[c4_i36] = 362436069U + 158852560U * (uint32_T)
          c4_i36;
      }

      chartInstance->c4_b_state_not_empty = true;
    }

    c4_icng = chartInstance->c4_b_state[0];
    c4_jsr = chartInstance->c4_b_state[1];
    c4_u0 = c4_jsr;
    c4_u1 = c4_icng;
    c4_u1 = 69069U * c4_u1 + 1234567U;
    c4_u0 ^= c4_u0 << 13;
    c4_u0 ^= c4_u0 >> 17;
    c4_u0 ^= c4_u0 << 5;
    c4_ui = c4_u1 + c4_u0;
    chartInstance->c4_b_state[0] = c4_u1;
    chartInstance->c4_b_state[1] = c4_u0;
    c4_b_ui = c4_ui;
    c4_b_r = (real_T)c4_b_ui * 2.328306436538696E-10;
    c4_d2 = c4_b_r;
    c4_r = c4_d2;
  } else {
    if (!chartInstance->c4_c_state_not_empty) {
      c4_eml_rand_mt19937ar(chartInstance, c4_uv2);
      for (c4_i37 = 0; c4_i37 < 625; c4_i37++) {
        chartInstance->c4_c_state[c4_i37] = c4_uv2[c4_i37];
      }

      chartInstance->c4_c_state_not_empty = true;
    }

    c4_d3 = c4_c_eml_rand_mt19937ar(chartInstance, chartInstance->c4_c_state);
    c4_r = c4_d3;
  }

  return c4_r;
}

static uint32_T c4_eml_rand_str2id(SFc4_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
  return 4U;
}

static uint32_T c4_b_eml_rand_str2id(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
  return 5U;
}

static uint32_T c4_c_eml_rand_str2id(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
  return 7U;
}

static void c4_eml_rand_mcg16807(SFc4_simulation1InstanceStruct *chartInstance,
  uint32_T c4_d_state, uint32_T *c4_e_state, real_T *c4_r)
{
  *c4_e_state = c4_d_state;
  *c4_r = c4_b_eml_rand_mcg16807(chartInstance, c4_e_state);
}

static void c4_eml_rand_mt19937ar(SFc4_simulation1InstanceStruct *chartInstance,
  uint32_T c4_d_state[625])
{
  int32_T c4_i38;
  for (c4_i38 = 0; c4_i38 < 625; c4_i38++) {
    c4_d_state[c4_i38] = 0U;
  }

  c4_b_twister_state_vector(chartInstance, c4_d_state, 5489.0);
}

static void c4_twister_state_vector(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T c4_mt[625], real_T c4_seed, uint32_T c4_b_mt[625])
{
  int32_T c4_i39;
  for (c4_i39 = 0; c4_i39 < 625; c4_i39++) {
    c4_b_mt[c4_i39] = c4_mt[c4_i39];
  }

  c4_b_twister_state_vector(chartInstance, c4_b_mt, c4_seed);
}

static void c4_b_eml_rand_mt19937ar(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T c4_d_state[625], uint32_T c4_e_state[625], real_T
  *c4_r)
{
  int32_T c4_i40;
  for (c4_i40 = 0; c4_i40 < 625; c4_i40++) {
    c4_e_state[c4_i40] = c4_d_state[c4_i40];
  }

  *c4_r = c4_c_eml_rand_mt19937ar(chartInstance, c4_e_state);
}

static void c4_eml_error(SFc4_simulation1InstanceStruct *chartInstance)
{
  int32_T c4_i41;
  static char_T c4_cv0[37] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'r', 'a', 'n', 'd', '_', 'i', 'n', 'v', 'a', 'l', 'i', 'd',
    'T', 'w', 'i', 's', 't', 'e', 'r', 'S', 't', 'a', 't', 'e' };

  char_T c4_b_u[37];
  const mxArray *c4_y = NULL;
  (void)chartInstance;
  for (c4_i41 = 0; c4_i41 < 37; c4_i41++) {
    c4_b_u[c4_i41] = c4_cv0[c4_i41];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_b_u, 10, 0U, 1U, 0U, 2, 1, 37),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c4_y));
}

static void c4_b_rand(SFc4_simulation1InstanceStruct *chartInstance, real_T
                      c4_r[3])
{
  uint32_T c4_V4;
  uint32_T c4_V5NORMAL;
  uint32_T c4_TWISTER;
  int32_T c4_k;
  real_T c4_b_k;
  uint32_T c4_d_state;
  real_T c4_d4;
  int32_T c4_i42;
  int32_T c4_c_k;
  real_T c4_d_k;
  uint32_T c4_icng;
  uint32_T c4_jsr;
  uint32_T c4_u2;
  uint32_T c4_u3;
  uint32_T c4_ui;
  uint32_T c4_b_ui;
  real_T c4_b_r;
  real_T c4_d5;
  uint32_T c4_uv3[625];
  int32_T c4_i43;
  int32_T c4_e_k;
  real_T c4_f_k;
  real_T c4_d6;
  c4_V4 = c4_eml_rand_str2id(chartInstance);
  c4_V5NORMAL = c4_b_eml_rand_str2id(chartInstance);
  c4_TWISTER = c4_c_eml_rand_str2id(chartInstance);
  if (!chartInstance->c4_method_not_empty) {
    chartInstance->c4_method = c4_TWISTER;
    chartInstance->c4_method_not_empty = true;
  }

  if (chartInstance->c4_method == c4_V4) {
    if (!chartInstance->c4_state_not_empty) {
      chartInstance->c4_state = 1144108930U;
      chartInstance->c4_state_not_empty = true;
    }

    for (c4_k = 0; c4_k < 3; c4_k++) {
      c4_b_k = 1.0 + (real_T)c4_k;
      c4_d_state = chartInstance->c4_state;
      c4_d4 = c4_b_eml_rand_mcg16807(chartInstance, &c4_d_state);
      chartInstance->c4_state = c4_d_state;
      c4_r[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_k), 1, 3, 1, 0) - 1] = c4_d4;
    }
  } else if (chartInstance->c4_method == c4_V5NORMAL) {
    if (!chartInstance->c4_b_state_not_empty) {
      for (c4_i42 = 0; c4_i42 < 2; c4_i42++) {
        chartInstance->c4_b_state[c4_i42] = 362436069U + 158852560U * (uint32_T)
          c4_i42;
      }

      chartInstance->c4_b_state_not_empty = true;
    }

    for (c4_c_k = 0; c4_c_k < 3; c4_c_k++) {
      c4_d_k = 1.0 + (real_T)c4_c_k;
      c4_icng = chartInstance->c4_b_state[0];
      c4_jsr = chartInstance->c4_b_state[1];
      c4_u2 = c4_jsr;
      c4_u3 = c4_icng;
      c4_u3 = 69069U * c4_u3 + 1234567U;
      c4_u2 ^= c4_u2 << 13;
      c4_u2 ^= c4_u2 >> 17;
      c4_u2 ^= c4_u2 << 5;
      c4_ui = c4_u3 + c4_u2;
      chartInstance->c4_b_state[0] = c4_u3;
      chartInstance->c4_b_state[1] = c4_u2;
      c4_b_ui = c4_ui;
      c4_b_r = (real_T)c4_b_ui * 2.328306436538696E-10;
      c4_d5 = c4_b_r;
      c4_r[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_d_k), 1, 3, 1, 0) - 1] = c4_d5;
    }
  } else {
    if (!chartInstance->c4_c_state_not_empty) {
      c4_eml_rand_mt19937ar(chartInstance, c4_uv3);
      for (c4_i43 = 0; c4_i43 < 625; c4_i43++) {
        chartInstance->c4_c_state[c4_i43] = c4_uv3[c4_i43];
      }

      chartInstance->c4_c_state_not_empty = true;
    }

    for (c4_e_k = 0; c4_e_k < 3; c4_e_k++) {
      c4_f_k = 1.0 + (real_T)c4_e_k;
      c4_d6 = c4_c_eml_rand_mt19937ar(chartInstance, chartInstance->c4_c_state);
      c4_r[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_f_k), 1, 3, 1, 0) - 1] = c4_d6;
    }
  }
}

static void c4_eml_scalar_eg(SFc4_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_eml_xgemm(SFc4_simulation1InstanceStruct *chartInstance, real_T
  c4_A[9], real_T c4_B[3], real_T c4_C[3], real_T c4_b_C[3])
{
  int32_T c4_i44;
  int32_T c4_i45;
  real_T c4_b_A[9];
  int32_T c4_i46;
  real_T c4_b_B[3];
  for (c4_i44 = 0; c4_i44 < 3; c4_i44++) {
    c4_b_C[c4_i44] = c4_C[c4_i44];
  }

  for (c4_i45 = 0; c4_i45 < 9; c4_i45++) {
    c4_b_A[c4_i45] = c4_A[c4_i45];
  }

  for (c4_i46 = 0; c4_i46 < 3; c4_i46++) {
    c4_b_B[c4_i46] = c4_B[c4_i46];
  }

  c4_b_eml_xgemm(chartInstance, c4_b_A, c4_b_B, c4_b_C);
}

static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_b_u;
  const mxArray *c4_y = NULL;
  SFc4_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc4_simulation1InstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_b_u = *(int32_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_b_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static int32_T c4_d_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId)
{
  int32_T c4_y;
  int32_T c4_i47;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), &c4_i47, 1, 6, 0U, 0, 0U, 0);
  c4_y = c4_i47;
  sf_mex_destroy(&c4_b_u);
  return c4_y;
}

static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_b_sfEvent;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  int32_T c4_y;
  SFc4_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc4_simulation1InstanceStruct *)chartInstanceVoid;
  c4_b_sfEvent = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_sfEvent),
    &c4_thisId);
  sf_mex_destroy(&c4_b_sfEvent);
  *(int32_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static void c4_e_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId,
  c4_sdz1UxT0ktEtEZMZuZvzEsH *c4_y)
{
  emlrtMsgIdentifier c4_thisId;
  static const char * c4_fieldNames[18] = { "thrustLimits", "thrustRateLimits",
    "thrustDir", "L", "nL", "rad", "bladeRad", "mass", "I", "km", "kf", "K",
    "velocityGain", "vGain", "zGain", "thetaGain", "omegaGain", "Ts2" };

  c4_thisId.fParent = c4_parentId;
  sf_mex_check_struct(c4_parentId, c4_b_u, 18, c4_fieldNames, 0U, NULL);
  c4_thisId.fIdentifier = "thrustLimits";
  c4_f_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u,
    "thrustLimits", "thrustLimits", 0)), &c4_thisId, c4_y->thrustLimits);
  c4_thisId.fIdentifier = "thrustRateLimits";
  c4_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u,
    "thrustRateLimits", "thrustRateLimits", 0)), &c4_thisId,
                        c4_y->thrustRateLimits);
  c4_thisId.fIdentifier = "thrustDir";
  c4_h_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u,
    "thrustDir", "thrustDir", 0)), &c4_thisId, c4_y->thrustDir);
  c4_thisId.fIdentifier = "L";
  c4_h_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u, "L",
    "L", 0)), &c4_thisId, c4_y->L);
  c4_thisId.fIdentifier = "nL";
  c4_y->nL = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c4_b_u, "nL", "nL", 0)), &c4_thisId);
  c4_thisId.fIdentifier = "rad";
  c4_y->rad = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c4_b_u, "rad", "rad", 0)), &c4_thisId);
  c4_thisId.fIdentifier = "bladeRad";
  c4_y->bladeRad = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c4_b_u, "bladeRad", "bladeRad", 0)), &c4_thisId);
  c4_thisId.fIdentifier = "mass";
  c4_y->mass = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c4_b_u, "mass", "mass", 0)), &c4_thisId);
  c4_thisId.fIdentifier = "I";
  c4_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u, "I",
    "I", 0)), &c4_thisId, c4_y->I);
  c4_thisId.fIdentifier = "km";
  c4_y->km = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c4_b_u, "km", "km", 0)), &c4_thisId);
  c4_thisId.fIdentifier = "kf";
  c4_y->kf = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c4_b_u, "kf", "kf", 0)), &c4_thisId);
  c4_thisId.fIdentifier = "K";
  c4_j_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u, "K",
    "K", 0)), &c4_thisId, c4_y->K);
  c4_thisId.fIdentifier = "velocityGain";
  c4_y->velocityGain = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c4_b_u, "velocityGain", "velocityGain", 0)), &c4_thisId);
  c4_thisId.fIdentifier = "vGain";
  c4_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u,
    "vGain", "vGain", 0)), &c4_thisId, c4_y->vGain);
  c4_thisId.fIdentifier = "zGain";
  c4_y->zGain = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c4_b_u, "zGain", "zGain", 0)), &c4_thisId);
  c4_thisId.fIdentifier = "thetaGain";
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u,
    "thetaGain", "thetaGain", 0)), &c4_thisId, c4_y->thetaGain);
  c4_thisId.fIdentifier = "omegaGain";
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c4_b_u,
    "omegaGain", "omegaGain", 0)), &c4_thisId, c4_y->omegaGain);
  c4_thisId.fIdentifier = "Ts2";
  c4_y->Ts2 = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c4_b_u, "Ts2", "Ts2", 0)), &c4_thisId);
  sf_mex_destroy(&c4_b_u);
}

static void c4_f_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[8])
{
  real_T c4_dv2[8];
  int32_T c4_i48;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_dv2, 1, 0, 0U, 1, 0U, 2, 2,
                4);
  for (c4_i48 = 0; c4_i48 < 8; c4_i48++) {
    c4_y[c4_i48] = c4_dv2[c4_i48];
  }

  sf_mex_destroy(&c4_b_u);
}

static void c4_g_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[4])
{
  real_T c4_dv3[4];
  int32_T c4_i49;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_dv3, 1, 0, 0U, 1, 0U, 2, 1,
                4);
  for (c4_i49 = 0; c4_i49 < 4; c4_i49++) {
    c4_y[c4_i49] = c4_dv3[c4_i49];
  }

  sf_mex_destroy(&c4_b_u);
}

static void c4_h_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[12])
{
  real_T c4_dv4[12];
  int32_T c4_i50;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_dv4, 1, 0, 0U, 1, 0U, 2, 3,
                4);
  for (c4_i50 = 0; c4_i50 < 12; c4_i50++) {
    c4_y[c4_i50] = c4_dv4[c4_i50];
  }

  sf_mex_destroy(&c4_b_u);
}

static void c4_i_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[9])
{
  real_T c4_dv5[9];
  int32_T c4_i51;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_dv5, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c4_i51 = 0; c4_i51 < 9; c4_i51++) {
    c4_y[c4_i51] = c4_dv5[c4_i51];
  }

  sf_mex_destroy(&c4_b_u);
}

static void c4_j_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[16])
{
  real_T c4_dv6[16];
  int32_T c4_i52;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_dv6, 1, 0, 0U, 1, 0U, 2, 4,
                4);
  for (c4_i52 = 0; c4_i52 < 16; c4_i52++) {
    c4_y[c4_i52] = c4_dv6[c4_i52];
  }

  sf_mex_destroy(&c4_b_u);
}

static void c4_k_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[2])
{
  real_T c4_dv7[2];
  int32_T c4_i53;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_dv7, 1, 0, 0U, 1, 0U, 1, 2);
  for (c4_i53 = 0; c4_i53 < 2; c4_i53++) {
    c4_y[c4_i53] = c4_dv7[c4_i53];
  }

  sf_mex_destroy(&c4_b_u);
}

static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_b_quad;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  c4_sdz1UxT0ktEtEZMZuZvzEsH c4_y;
  SFc4_simulation1InstanceStruct *chartInstance;
  chartInstance = (SFc4_simulation1InstanceStruct *)chartInstanceVoid;
  c4_b_quad = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_quad), &c4_thisId, &c4_y);
  sf_mex_destroy(&c4_b_quad);
  *(c4_sdz1UxT0ktEtEZMZuZvzEsH *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static uint32_T c4_l_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_method, const char_T *c4_identifier)
{
  uint32_T c4_y;
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_method),
    &c4_thisId);
  sf_mex_destroy(&c4_b_method);
  return c4_y;
}

static uint32_T c4_m_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId)
{
  uint32_T c4_y;
  uint32_T c4_u4;
  if (mxIsEmpty(c4_b_u)) {
    chartInstance->c4_method_not_empty = false;
  } else {
    chartInstance->c4_method_not_empty = true;
    sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), &c4_u4, 1, 7, 0U, 0, 0U, 0);
    c4_y = c4_u4;
  }

  sf_mex_destroy(&c4_b_u);
  return c4_y;
}

static uint32_T c4_n_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_d_state, const char_T *c4_identifier)
{
  uint32_T c4_y;
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_d_state), &c4_thisId);
  sf_mex_destroy(&c4_d_state);
  return c4_y;
}

static uint32_T c4_o_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId)
{
  uint32_T c4_y;
  uint32_T c4_u5;
  if (mxIsEmpty(c4_b_u)) {
    chartInstance->c4_state_not_empty = false;
  } else {
    chartInstance->c4_state_not_empty = true;
    sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), &c4_u5, 1, 7, 0U, 0, 0U, 0);
    c4_y = c4_u5;
  }

  sf_mex_destroy(&c4_b_u);
  return c4_y;
}

static void c4_p_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_d_state, const char_T *c4_identifier, uint32_T c4_y[625])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_q_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_d_state), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_d_state);
}

static void c4_q_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, uint32_T c4_y
  [625])
{
  uint32_T c4_uv4[625];
  int32_T c4_i54;
  if (mxIsEmpty(c4_b_u)) {
    chartInstance->c4_c_state_not_empty = false;
  } else {
    chartInstance->c4_c_state_not_empty = true;
    sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_uv4, 1, 7, 0U, 1, 0U, 1,
                  625);
    for (c4_i54 = 0; c4_i54 < 625; c4_i54++) {
      c4_y[c4_i54] = c4_uv4[c4_i54];
    }
  }

  sf_mex_destroy(&c4_b_u);
}

static void c4_r_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_d_state, const char_T *c4_identifier, uint32_T c4_y[2])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_s_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_d_state), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_d_state);
}

static void c4_s_emlrt_marshallIn(SFc4_simulation1InstanceStruct *chartInstance,
  const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId, uint32_T c4_y[2])
{
  uint32_T c4_uv5[2];
  int32_T c4_i55;
  if (mxIsEmpty(c4_b_u)) {
    chartInstance->c4_b_state_not_empty = false;
  } else {
    chartInstance->c4_b_state_not_empty = true;
    sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), c4_uv5, 1, 7, 0U, 1, 0U, 1, 2);
    for (c4_i55 = 0; c4_i55 < 2; c4_i55++) {
      c4_y[c4_i55] = c4_uv5[c4_i55];
    }
  }

  sf_mex_destroy(&c4_b_u);
}

static uint8_T c4_t_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_is_active_c4_simulation1, const char_T
  *c4_identifier)
{
  uint8_T c4_y;
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_u_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c4_b_is_active_c4_simulation1), &c4_thisId);
  sf_mex_destroy(&c4_b_is_active_c4_simulation1);
  return c4_y;
}

static uint8_T c4_u_emlrt_marshallIn(SFc4_simulation1InstanceStruct
  *chartInstance, const mxArray *c4_b_u, const emlrtMsgIdentifier *c4_parentId)
{
  uint8_T c4_y;
  uint8_T c4_u6;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_b_u), &c4_u6, 1, 3, 0U, 0, 0U, 0);
  c4_y = c4_u6;
  sf_mex_destroy(&c4_b_u);
  return c4_y;
}

static real_T c4_b_eml_rand_mcg16807(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T *c4_d_state)
{
  uint32_T c4_s;
  uint32_T c4_u7;
  uint32_T c4_hi;
  uint32_T c4_lo;
  uint32_T c4_test1;
  uint32_T c4_test2;
  uint32_T c4_e_state;
  real_T c4_b_r;
  (void)chartInstance;
  c4_s = *c4_d_state;
  c4_u7 = 127773U;
  if (c4_u7 == 0U) {
    c4_hi = MAX_uint32_T;
  } else {
    c4_hi = c4_s / c4_u7;
  }

  c4_lo = c4_s - c4_hi * 127773U;
  c4_test1 = 16807U * c4_lo;
  c4_test2 = 2836U * c4_hi;
  if (c4_test1 < c4_test2) {
    c4_e_state = (c4_test1 - c4_test2) + 2147483647U;
  } else {
    c4_e_state = c4_test1 - c4_test2;
  }

  c4_b_r = (real_T)c4_e_state * 4.6566128752457969E-10;
  *c4_d_state = c4_e_state;
  return c4_b_r;
}

static void c4_b_twister_state_vector(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T c4_mt[625], real_T c4_seed)
{
  real_T c4_d7;
  uint32_T c4_u8;
  uint32_T c4_r;
  int32_T c4_mti;
  real_T c4_b_mti;
  real_T c4_d8;
  uint32_T c4_u9;
  (void)chartInstance;
  c4_d7 = c4_seed;
  if (c4_d7 < 4.294967296E+9) {
    if (c4_d7 >= 0.0) {
      c4_u8 = (uint32_T)c4_d7;
    } else {
      c4_u8 = 0U;
    }
  } else if (c4_d7 >= 4.294967296E+9) {
    c4_u8 = MAX_uint32_T;
  } else {
    c4_u8 = 0U;
  }

  c4_r = c4_u8;
  c4_mt[0] = c4_r;
  for (c4_mti = 0; c4_mti < 623; c4_mti++) {
    c4_b_mti = 2.0 + (real_T)c4_mti;
    c4_d8 = muDoubleScalarRound(c4_b_mti - 1.0);
    if (c4_d8 < 4.294967296E+9) {
      if (c4_d8 >= 0.0) {
        c4_u9 = (uint32_T)c4_d8;
      } else {
        c4_u9 = 0U;
      }
    } else if (c4_d8 >= 4.294967296E+9) {
      c4_u9 = MAX_uint32_T;
    } else {
      c4_u9 = 0U;
    }

    c4_r = (c4_r ^ c4_r >> 30U) * 1812433253U + c4_u9;
    c4_mt[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c4_b_mti), 1, 625, 1, 0) - 1] = c4_r;
  }

  c4_mt[624] = 624U;
}

static real_T c4_c_eml_rand_mt19937ar(SFc4_simulation1InstanceStruct
  *chartInstance, uint32_T c4_d_state[625])
{
  int32_T c4_i56;
  uint32_T c4_b_u[2];
  int32_T c4_j;
  real_T c4_b_j;
  uint32_T c4_mti;
  int32_T c4_kk;
  real_T c4_b_kk;
  uint32_T c4_y;
  uint32_T c4_b_y;
  uint32_T c4_c_y;
  int32_T c4_c_kk;
  uint32_T c4_d_y;
  uint32_T c4_e_y;
  uint32_T c4_f_y;
  uint32_T c4_g_y;
  real_T c4_b_r;
  boolean_T c4_b0;
  boolean_T c4_isvalid;
  int32_T c4_k;
  int32_T c4_a;
  int32_T c4_b_a;
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
    for (c4_i56 = 0; c4_i56 < 2; c4_i56++) {
      c4_b_u[c4_i56] = 0U;
    }

    for (c4_j = 0; c4_j < 2; c4_j++) {
      c4_b_j = 1.0 + (real_T)c4_j;
      c4_mti = c4_d_state[624] + 1U;
      if ((real_T)c4_mti >= 625.0) {
        for (c4_kk = 0; c4_kk < 227; c4_kk++) {
          c4_b_kk = 1.0 + (real_T)c4_kk;
          c4_y = (c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                   _SFD_INTEGER_CHECK("", c4_b_kk), 1, 625, 1, 0) - 1] &
                  2147483648U) | (c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", c4_b_kk + 1.0), 1, 625, 1, 0) - 1] &
            2147483647U);
          c4_b_y = c4_y;
          c4_c_y = c4_b_y;
          if ((real_T)(c4_c_y & 1U) == 0.0) {
            c4_c_y >>= 1U;
          } else {
            c4_c_y = c4_c_y >> 1U ^ 2567483615U;
          }

          c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c4_b_kk), 1, 625, 1, 0) - 1] =
            c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_kk + 397.0), 1, 625, 1, 0) - 1] ^ c4_c_y;
        }

        for (c4_c_kk = 0; c4_c_kk < 396; c4_c_kk++) {
          c4_b_kk = 228.0 + (real_T)c4_c_kk;
          c4_y = (c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                   _SFD_INTEGER_CHECK("", c4_b_kk), 1, 625, 1, 0) - 1] &
                  2147483648U) | (c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", c4_b_kk + 1.0), 1, 625, 1, 0) - 1] &
            2147483647U);
          c4_d_y = c4_y;
          c4_e_y = c4_d_y;
          if ((real_T)(c4_e_y & 1U) == 0.0) {
            c4_e_y >>= 1U;
          } else {
            c4_e_y = c4_e_y >> 1U ^ 2567483615U;
          }

          c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c4_b_kk), 1, 625, 1, 0) - 1] =
            c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (c4_b_kk + 1.0) - 228.0), 1, 625, 1, 0) - 1] ^
            c4_e_y;
        }

        c4_y = (c4_d_state[623] & 2147483648U) | (c4_d_state[0] & 2147483647U);
        c4_f_y = c4_y;
        c4_g_y = c4_f_y;
        if ((real_T)(c4_g_y & 1U) == 0.0) {
          c4_g_y >>= 1U;
        } else {
          c4_g_y = c4_g_y >> 1U ^ 2567483615U;
        }

        c4_d_state[623] = c4_d_state[396] ^ c4_g_y;
        c4_mti = 1U;
      }

      c4_y = c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)(uint32_T)
        _SFD_INTEGER_CHECK("", (real_T)c4_mti), 1, 625, 1, 0) - 1];
      c4_d_state[624] = c4_mti;
      c4_y ^= c4_y >> 11U;
      c4_y ^= c4_y << 7U & 2636928640U;
      c4_y ^= c4_y << 15U & 4022730752U;
      c4_y ^= c4_y >> 18U;
      c4_b_u[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_j), 1, 2, 1, 0) - 1] = c4_y;
    }

    c4_b_u[0] >>= 5U;
    c4_b_u[1] >>= 6U;
    c4_b_r = 1.1102230246251565E-16 * ((real_T)c4_b_u[0] * 6.7108864E+7 +
      (real_T)c4_b_u[1]);
    if (c4_b_r == 0.0) {
      guard1 = false;
      if ((real_T)c4_d_state[624] >= 1.0) {
        if ((real_T)c4_d_state[624] < 625.0) {
          c4_b0 = true;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1 == true) {
        c4_b0 = false;
      }

      c4_isvalid = c4_b0;
      if (c4_isvalid) {
        c4_isvalid = false;
        c4_k = 1;
        exitg2 = false;
        while ((exitg2 == false) && (c4_k < 625)) {
          if ((real_T)c4_d_state[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
               _SFD_INTEGER_CHECK("", (real_T)c4_k), 1, 625, 1, 0) - 1] == 0.0)
          {
            c4_a = c4_k;
            c4_b_a = c4_a + 1;
            c4_k = c4_b_a;
          } else {
            c4_isvalid = true;
            exitg2 = true;
          }
        }
      }

      if (!c4_isvalid) {
        c4_eml_error(chartInstance);
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return c4_b_r;
}

static void c4_b_eml_xgemm(SFc4_simulation1InstanceStruct *chartInstance, real_T
  c4_A[9], real_T c4_B[3], real_T c4_C[3])
{
  int32_T c4_i57;
  int32_T c4_i58;
  int32_T c4_i59;
  (void)chartInstance;
  for (c4_i57 = 0; c4_i57 < 3; c4_i57++) {
    c4_C[c4_i57] = 0.0;
    c4_i58 = 0;
    for (c4_i59 = 0; c4_i59 < 3; c4_i59++) {
      c4_C[c4_i57] += c4_A[c4_i58 + c4_i57] * c4_B[c4_i59];
      c4_i58 += 3;
    }
  }
}

static void init_dsm_address_info(SFc4_simulation1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc4_simulation1InstanceStruct
  *chartInstance)
{
  chartInstance->c4_u = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c4_dist = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c4_omega = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c4_omegaDot = (real_T (*)[3])ssGetOutputPortSignal_wrapper
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

void sf_c4_simulation1_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4001295465U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2800362758U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(256938266U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3627084120U);
}

mxArray* sf_c4_simulation1_get_post_codegen_info(void);
mxArray *sf_c4_simulation1_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("Gfk94EFMIxP6qtZDRQUp0F");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

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
    mxArray* mxPostCodegenInfo = sf_c4_simulation1_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c4_simulation1_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c4_simulation1_jit_fallback_info(void)
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

mxArray *sf_c4_simulation1_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c4_simulation1_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c4_simulation1(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x6'type','srcId','name','auxInfo'{{M[1],M[5],T\"omegaDot\",},{M[4],M[0],T\"method\",S'l','i','p'{{M1x2[512 518],M[1],T\"/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/randfun/eml_rand.m\"}}},{M[4],M[0],T\"state\",S'l','i','p'{{M1x2[165 170],M[1],T\"/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/randfun/eml_rand_mcg16807_stateful.m\"}}},{M[4],M[0],T\"state\",S'l','i','p'{{M1x2[166 171],M[1],T\"/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/randfun/eml_rand_mt19937ar_stateful.m\"}}},{M[4],M[0],T\"state\",S'l','i','p'{{M1x2[165 170],M[1],T\"/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/randfun/eml_rand_shr3cong_stateful.m\"}}},{M[8],M[0],T\"is_active_c4_simulation1\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 6, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c4_simulation1_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc4_simulation1InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc4_simulation1InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulation1MachineNumber_,
           4,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"u");
          _SFD_SET_DATA_PROPS(1,1,1,0,"dist");
          _SFD_SET_DATA_PROPS(2,1,1,0,"omega");
          _SFD_SET_DATA_PROPS(3,10,0,0,"quad");
          _SFD_SET_DATA_PROPS(4,2,0,1,"omegaDot");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,245);
        _SFD_CV_INIT_EML_IF(0,1,0,60,68,146,171);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(3,SF_STRUCT,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_c_sf_marshallOut,(MexInFcnForType)c4_d_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)
            c4_sf_marshallIn);
        }

        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c4_u);
        _SFD_SET_DATA_VALUE_PTR(1U, chartInstance->c4_dist);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c4_omega);
        _SFD_SET_DATA_VALUE_PTR(3U, &chartInstance->c4_quad);
        _SFD_SET_DATA_VALUE_PTR(4U, *chartInstance->c4_omegaDot);
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
  return "g6YE2gjTYnsVRBxQV6byvF";
}

static void sf_opaque_initialize_c4_simulation1(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc4_simulation1InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c4_simulation1((SFc4_simulation1InstanceStruct*)
    chartInstanceVar);
  initialize_c4_simulation1((SFc4_simulation1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c4_simulation1(void *chartInstanceVar)
{
  enable_c4_simulation1((SFc4_simulation1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c4_simulation1(void *chartInstanceVar)
{
  disable_c4_simulation1((SFc4_simulation1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c4_simulation1(void *chartInstanceVar)
{
  sf_gateway_c4_simulation1((SFc4_simulation1InstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c4_simulation1(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c4_simulation1((SFc4_simulation1InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c4_simulation1(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c4_simulation1((SFc4_simulation1InstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c4_simulation1(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc4_simulation1InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation1_optimization_info();
    }

    finalize_c4_simulation1((SFc4_simulation1InstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc4_simulation1((SFc4_simulation1InstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c4_simulation1(SimStruct *S)
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
    initialize_params_c4_simulation1((SFc4_simulation1InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c4_simulation1(SimStruct *S)
{
  /* Actual parameters from chart:
     quad
   */
  const char_T *rtParamNames[] = { "quad" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation1_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,4);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,4,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,4,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,4);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,4,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,4,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,4);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1681390798U));
  ssSetChecksum1(S,(997000517U));
  ssSetChecksum2(S,(3551634258U));
  ssSetChecksum3(S,(3754630092U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c4_simulation1(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c4_simulation1(SimStruct *S)
{
  SFc4_simulation1InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc4_simulation1InstanceStruct *)utMalloc(sizeof
    (SFc4_simulation1InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc4_simulation1InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c4_simulation1;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c4_simulation1;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c4_simulation1;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c4_simulation1;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c4_simulation1;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c4_simulation1;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c4_simulation1;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c4_simulation1;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c4_simulation1;
  chartInstance->chartInfo.mdlStart = mdlStart_c4_simulation1;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c4_simulation1;
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

void c4_simulation1_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c4_simulation1(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c4_simulation1(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c4_simulation1(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c4_simulation1_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
