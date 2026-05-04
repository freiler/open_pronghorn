//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED

#include "OpenPronghornLinearFVTurbulencePhysics.h"

#include "WCNSFVFlowPhysicsBase.h"
#include "NS.h"

#include <iomanip>
#include <sstream>

registerWCNSFVTurbulenceBaseTasks("OpenPronghornApp", OpenPronghornLinearFVTurbulencePhysics);
registerMooseAction("OpenPronghornApp",
                    OpenPronghornLinearFVTurbulencePhysics,
                    "add_functor_material");

InputParameters
OpenPronghornLinearFVTurbulencePhysics::validParams()
{
  InputParameters params = WCNSFVTurbulencePhysicsBase::validParams();
  params.addClassDescription(
      "Define an OpenPronghorn k-epsilon turbulence model for an incompressible or "
      "weakly-compressible Navier-Stokes flow with a linearly-segregated finite-volume "
      "discretization. Supports Standard, LowRe, TwoLayer and Realizable variants and "
      "optional corrections (Yap, buoyancy, compressibility, curvature).");

  // Default to k-epsilon; users should not normally change this for this physics class
  params.set<MooseEnum>("turbulence_handling") = "k-epsilon";

  // Solver system names for the two additional linear systems
  params.addParam<std::vector<SolverSystemName>>(
      "system_names",
      {"TKE_system", "TKED_system"},
      "Names of the linear solver systems for TKE and epsilon. Adapt when using non-default "
      "system naming in the Problem block.");

  // mu_t must be an auxiliary variable for the linear (SIMPLE) solver
  params.set<bool>("mu_t_as_aux_variable") = true;
  params.suppressParameter<bool>("mu_t_as_aux_variable");
  params.suppressParameter<bool>("turbulent_viscosity_two_term_bc_expansion");
  params.suppressParameter<bool>("tke_two_term_bc_expansion");
  params.suppressParameter<bool>("tked_two_term_bc_expansion");
  params.suppressParameter<MooseEnum>("wall_treatment_T");
  params.suppressParameter<MooseEnum>("tke_face_interpolation");
  params.suppressParameter<MooseEnum>("tked_face_interpolation");

  // Set sensible defaults for base-class k-epsilon coefficients (no defaults in base)
  params.set<MooseFunctorName>("C1_eps") = "1.44";
  params.set<MooseFunctorName>("C2_eps") = "1.92";
  params.set<MooseFunctorName>("sigma_k") = "1.0";
  params.set<MooseFunctorName>("sigma_eps") = "1.3";

  // ---- LinearFV numerical scheme ----
  params.addParam<bool>(
      "use_nonorthogonal_correction",
      true,
      "Whether to apply a non-orthogonal correction to the diffusion stencil. "
      "Reduces dispersion on skewed meshes; can safely be disabled on orthogonal meshes.");
  params.addParamNamesToGroup("use_nonorthogonal_correction", "Numerical scheme");

  // ---- Wall distance ----
  params.addParam<VariableName>(
      "wall_distance_name",
      "wall_distance",
      "Name of the auxiliary variable that stores cell-center-to-wall distances. "
      "Required by TwoLayer and LowRe variants; also used for Yap correction.");
  params.addParamNamesToGroup("wall_distance_name", "Wall treatment");

  // ---- Inlet boundary conditions ----
  // The inlet boundaries come from the coupled flow physics; one type must be listed per boundary.
  MultiMooseEnum inlet_bc_types("intensity-length-scale fixed-value");
  params.addParam<MultiMooseEnum>(
      "turbulence_inlet_types",
      inlet_bc_types,
      "Type of TKE/epsilon boundary condition at each inlet boundary defined in the coupled flow "
      "physics (same order). 'intensity-length-scale': compute k = 1.5*(I*U_ref)^2 and "
      "eps = C_mu^(3/4)*k^(3/2)/L from turbulence_intensity, turbulence_reference_velocity, "
      "and turbulence_length_scale. 'fixed-value': prescribe explicit functor values via "
      "tke_inlet_functors and tked_inlet_functors.");
  params.addParam<std::vector<Real>>(
      "turbulence_intensity",
      "Turbulence intensity I (dimensionless, e.g. 0.05) for each 'intensity-length-scale' inlet "
      "boundary, listed in flow-physics inlet-boundary order.");
  params.addParam<std::vector<Real>>(
      "turbulence_length_scale",
      "Turbulence length scale L (same length units as the mesh) for each "
      "'intensity-length-scale' inlet boundary.");
  params.addParam<std::vector<Real>>(
      "turbulence_reference_velocity",
      "Reference velocity magnitude U_ref for each 'intensity-length-scale' inlet boundary.");
  params.addParam<std::vector<MooseFunctorName>>(
      "tke_inlet_functors",
      "TKE functor value at each 'fixed-value' inlet boundary, in flow-physics inlet order.");
  params.addParam<std::vector<MooseFunctorName>>(
      "tked_inlet_functors",
      "Epsilon functor value at each 'fixed-value' inlet boundary, in flow-physics inlet order.");
  params.addParamNamesToGroup(
      "turbulence_inlet_types turbulence_intensity turbulence_length_scale "
      "turbulence_reference_velocity tke_inlet_functors tked_inlet_functors",
      "Inlet boundary conditions");

  // ---- Outlet boundary conditions ----
  params.addParam<bool>(
      "add_turbulence_outflow_bcs",
      false,
      "If true, add zero-gradient (LinearFVAdvectionDiffusionOutflowBC) for TKE and epsilon at "
      "all outlet boundaries defined in the coupled flow physics.");
  params.addParamNamesToGroup("add_turbulence_outflow_bcs", "Outlet boundary conditions");

  // ---- k-epsilon model variant ----
  MooseEnum ke_variant(
      "Standard StandardLowRe StandardTwoLayer Realizable RealizableTwoLayer", "Standard");
  params.addParam<MooseEnum>(
      "k_epsilon_variant",
      ke_variant,
      "k-epsilon model variant used for turbulent viscosity and source/sink terms.");

  MooseEnum two_layer_flavor("Wolfstein NorrisReynolds Xu", "Wolfstein");
  params.addParam<MooseEnum>(
      "two_layer_flavor",
      two_layer_flavor,
      "Two-layer length-scale formulation for TwoLayer variants (Wolfstein, NorrisReynolds, Xu).");

  MooseEnum scale_limiter("none standard", "standard");
  params.addParam<MooseEnum>(
      "scale_limiter",
      scale_limiter,
      "Time-scale limiter for turbulent viscosity: 'none' or 'standard' "
      "(max(Te, Ct*sqrt(nu/epsilon))).");

  params.addParamNamesToGroup("k_epsilon_variant two_layer_flavor scale_limiter",
                              "k-epsilon variant");

  // ---- Optional correction switches ----
  params.addParam<bool>("use_yap", false, "Include the Yap correction term in the epsilon equation.");
  params.addParam<bool>(
      "use_low_re_Gprime", false, "Include low-Re extra production G' in the epsilon equation.");
  params.addParam<bool>("use_buoyancy", false, "Include buoyancy production Gb.");
  params.addParam<bool>("use_compressibility",
                        false,
                        "Include the compressibility correction gamma_M in the k equation.");
  params.addParam<bool>("use_curvature_correction",
                        false,
                        "Apply the Spalart-Shur curvature/rotation correction to shear production.");
  params.addParamNamesToGroup(
      "use_yap use_low_re_Gprime use_buoyancy use_compressibility use_curvature_correction",
      "k-epsilon corrections");

  // ---- Correction constants ----
  params.addParam<Real>("Ct",
                        6.0,
                        "Ct coefficient used in the Yap correction and ambient epsilon source "
                        "term in the kEpsilonTKEDSourceSink kernel.");
  params.addParam<Real>("C_lowRe", 1.0, "Low-Re f2 constant C for the standard LowRe variant.");
  params.addParam<Real>("D_lowRe", 1.0, "Low-Re extra-production coefficient D for G'.");
  params.addParam<Real>("E_lowRe", 0.00375, "Low-Re extra-production coefficient E for G'.");
  params.addParam<Real>("Cw", 0.83, "Yap correction constant Cw.");
  params.addParam<Real>("C3_eps", 1.0, "C3 epsilon buoyancy coefficient.");
  params.addParam<Real>("C_M", 1.0, "Coefficient for the compressibility correction gamma_M.");
  params.addParamNamesToGroup("Ct C_lowRe D_lowRe E_lowRe Cw C3_eps C_M",
                              "Correction constants");

  // ---- Low-Re viscosity damping (kEpsilonViscosity, Standard LowRe variant) ----
  params.addParam<Real>(
      "Cd0", 0.091, "Low-Re f_mu damping coefficient Cd0 (Standard LowRe variant).");
  params.addParam<Real>(
      "Cd1", 0.0042, "Low-Re f_mu damping coefficient Cd1 (Standard LowRe variant).");
  params.addParam<Real>(
      "Cd2", 0.00011, "Low-Re f_mu damping coefficient Cd2 (Standard LowRe variant).");
  params.addParamNamesToGroup("Cd0 Cd1 Cd2", "Low-Re viscosity coefficients");

  // ---- Realizable viscosity coefficients (kEpsilonViscosity) ----
  params.addParam<Real>("Ca0", 0.667, "Realizable C_mu coefficient Ca0.");
  params.addParam<Real>("Ca1", 1.25, "Realizable C_mu coefficient Ca1.");
  params.addParam<Real>("Ca2", 1.0, "Realizable C_mu coefficient Ca2.");
  params.addParam<Real>("Ca3", 0.9, "Realizable C_mu coefficient Ca3.");
  params.addParamNamesToGroup("Ca0 Ca1 Ca2 Ca3", "Realizable viscosity coefficients");

  // ---- Non-linear constitutive model ----
  MooseEnum nonlinear_model("none quadratic cubic", "none");
  params.addParam<MooseEnum>(
      "nonlinear_model",
      nonlinear_model,
      "Non-linear constitutive relation for shear production: none, quadratic, or cubic.");

  MooseEnum curvature_model("none standard", "none");
  params.addParam<MooseEnum>(
      "curvature_model",
      curvature_model,
      "Curvature/rotation correction model: 'none' or 'standard' (Spalart-Shur style).");
  params.addParamNamesToGroup("nonlinear_model curvature_model", "k-epsilon corrections");

  // ---- Buoyancy and compressibility functors ----
  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0.0, 0.0, 0.0), "Gravity vector for buoyancy production.");
  params.addParam<MooseFunctorName>(
      "temperature", "Temperature functor; required when use_buoyancy = true.");
  params.addParam<MooseFunctorName>(
      "beta", "Thermal expansion coefficient functor; required when use_buoyancy = true.");
  params.addParam<MooseFunctorName>(
      "speed_of_sound",
      "Speed-of-sound functor; required when use_compressibility = true.");
  params.addParamNamesToGroup(
      "gravity temperature beta speed_of_sound C_M", "Buoyancy and compressibility");

  return params;
}

OpenPronghornLinearFVTurbulencePhysics::OpenPronghornLinearFVTurbulencePhysics(
    const InputParameters & parameters)
  : WCNSFVTurbulencePhysicsBase(parameters),
    _wall_distance_name(getParam<VariableName>("wall_distance_name"))
{
  if (_turbulence_model != "k-epsilon" && _turbulence_model != "none")
    paramError("turbulence_handling",
               "OpenPronghornLinearFVTurbulencePhysics only supports the 'k-epsilon' (or 'none') "
               "turbulence model.");
}

void
OpenPronghornLinearFVTurbulencePhysics::initializePhysicsAdditional()
{
  if (_turbulence_model == "k-epsilon")
    getProblem().needSolutionState(1, Moose::SolutionIterationType::Nonlinear);
}

void
OpenPronghornLinearFVTurbulencePhysics::checkIntegrity() const
{
  WCNSFVTurbulencePhysicsBase::checkIntegrity();

  if (_flow_equations_physics &&
      !_flow_equations_physics->getParam<bool>("include_symmetrized_viscous_stress"))
    _flow_equations_physics->paramWarning(
        "include_symmetrized_viscous_stress",
        "This should be set to true when using a turbulence model.");
}

unsigned short
OpenPronghornLinearFVTurbulencePhysics::getNumberAlgebraicGhostingLayersNeeded() const
{
  return _flow_equations_physics->getNumberAlgebraicGhostingLayersNeeded();
}

// ---------------------------------------------------------------------------
// Variable creation
// ---------------------------------------------------------------------------

void
OpenPronghornLinearFVTurbulencePhysics::addSolverVariables()
{
  if (_turbulence_model != "k-epsilon")
    return;

  for (const auto & [var_name, param_name] :
       std::vector<std::pair<VariableName, std::string>>{{_tke_name, "tke_name"},
                                                         {_tked_name, "tked_name"}})
  {
    if (!shouldCreateVariable(var_name, _blocks, /*error_if_aux=*/true))
      reportPotentiallyMissedParameters({"system_names"}, "MooseLinearVariableFVReal");
    else if (_define_variables)
    {
      auto vparams = getFactory().getValidParams("MooseLinearVariableFVReal");
      assignBlocks(vparams, _blocks);
      vparams.set<SolverSystemName>("solver_sys") = getSolverSystem(var_name);
      getProblem().addVariable("MooseLinearVariableFVReal", var_name, vparams);
    }
    else
      paramError(param_name,
                 "Variable (" + var_name +
                     ") supplied to OpenPronghornLinearFVTurbulencePhysics does not exist.");
  }
}

void
OpenPronghornLinearFVTurbulencePhysics::addAuxiliaryVariables()
{
  if (_turbulence_model != "k-epsilon")
    return;

  // Turbulent viscosity: MooseLinearVariableFVReal so it is visible to linear BCs
  if (shouldCreateVariable(_turbulent_viscosity_name, _blocks, /*error_if_aux=*/false))
  {
    auto vparams = getFactory().getValidParams("MooseLinearVariableFVReal");
    assignBlocks(vparams, _blocks);
    getProblem().addAuxVariable("MooseLinearVariableFVReal", _turbulent_viscosity_name, vparams);
  }

  // Wall distance: MooseVariableFVReal required by WallDistanceAux
  if (shouldCreateVariable(_wall_distance_name, _blocks, /*error_if_aux=*/false))
  {
    auto vparams = getFactory().getValidParams("MooseVariableFVReal");
    assignBlocks(vparams, _blocks);
    getProblem().addAuxVariable("MooseVariableFVReal", _wall_distance_name, vparams);
  }
}

void
OpenPronghornLinearFVTurbulencePhysics::addAuxiliaryKernels()
{
  if (_turbulence_model != "k-epsilon")
    return;

  const std::string u_names[3] = {"u", "v", "w"};

  // WallDistanceAux: compute distance to turbulence walls at initial time and each nonlinear iter
  {
    auto params = getFactory().getValidParams("WallDistanceAux");
    assignBlocks(params, _blocks);
    params.set<AuxVariableName>("variable") = _wall_distance_name;
    params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_NONLINEAR};
    getProblem().addAuxKernel("WallDistanceAux", name() + "_wall_distance_aux", params);
  }

  // kEpsilonViscosity: compute mu_t at each nonlinear iteration
  {
    auto params = getFactory().getValidParams("kEpsilonViscosity");
    assignBlocks(params, _blocks);
    params.set<AuxVariableName>("variable") = _turbulent_viscosity_name;
    params.set<MooseFunctorName>(NS::density) = _flow_equations_physics->densityName();
    params.set<MooseFunctorName>(NS::mu) = _flow_equations_physics->dynamicViscosityName();
    params.set<MooseFunctorName>(NS::TKE) = _tke_name;
    params.set<MooseFunctorName>(NS::TKED) = _tked_name;
    params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
    params.set<bool>("bulk_wall_treatment") = getParam<bool>("bulk_wall_treatment");
    params.set<MooseEnum>("wall_treatment") = _wall_treatment_eps;
    params.set<Real>("C_mu") = getParam<Real>("C_mu");
    params.set<Real>("mu_t_ratio_max") = getParam<Real>("mu_t_ratio_max");
    const MooseEnum variant = getParam<MooseEnum>("k_epsilon_variant");
    params.set<MooseEnum>("k_epsilon_variant") = variant;
    params.set<MooseEnum>("scale_limiter") = getParam<MooseEnum>("scale_limiter");
    // Only set variant-specific parameters so kEpsilonViscosity applicability checks pass
    const bool is_two_layer = (variant == "StandardTwoLayer" || variant == "RealizableTwoLayer");
    const bool is_lowre = (variant == "StandardLowRe");
    const bool is_realizable = (variant == "Realizable" || variant == "RealizableTwoLayer");
    if (is_two_layer)
      params.set<MooseEnum>("two_layer_flavor") = getParam<MooseEnum>("two_layer_flavor");
    if (is_lowre)
    {
      params.set<Real>("Cd0") = getParam<Real>("Cd0");
      params.set<Real>("Cd1") = getParam<Real>("Cd1");
      params.set<Real>("Cd2") = getParam<Real>("Cd2");
    }
    if (is_realizable)
    {
      params.set<Real>("Ca0") = getParam<Real>("Ca0");
      params.set<Real>("Ca1") = getParam<Real>("Ca1");
      params.set<Real>("Ca2") = getParam<Real>("Ca2");
      params.set<Real>("Ca3") = getParam<Real>("Ca3");
    }
    params.set<MooseFunctorName>("wall_distance") = _wall_distance_name;
    for (const auto d : make_range(dimension()))
      params.set<MooseFunctorName>(u_names[d]) = _velocity_names[d];
    params.set<ExecFlagEnum>("execute_on") = {EXEC_NONLINEAR};
    getProblem().addAuxKernel("kEpsilonViscosity", name() + "_mu_t_aux", params);
  }
}

// ---------------------------------------------------------------------------
// FV kernels
// ---------------------------------------------------------------------------

void
OpenPronghornLinearFVTurbulencePhysics::addFVKernels()
{
  if (_turbulence_model != "k-epsilon")
    return;

  if (isTransient())
    addKEpsilonTimeDerivatives();
  addKEpsilonAdvection();
  addKEpsilonDiffusion();
  addKEpsilonSourceSink();
}

void
OpenPronghornLinearFVTurbulencePhysics::addKEpsilonTimeDerivatives()
{
  const std::string kernel_type = "LinearFVTimeDerivative";
  InputParameters params = getFactory().getValidParams(kernel_type);
  assignBlocks(params, _blocks);
  params.set<MooseFunctorName>("factor") = _flow_equations_physics->densityName();

  params.set<LinearVariableName>("variable") = _tke_name;
  if (shouldCreateTimeDerivative(_tke_name, _blocks, /*error_if_defined=*/false))
    getProblem().addLinearFVKernel(kernel_type, prefix() + "tke_time", params);

  params.set<LinearVariableName>("variable") = _tked_name;
  if (shouldCreateTimeDerivative(_tked_name, _blocks, /*error_if_defined=*/false))
    getProblem().addLinearFVKernel(kernel_type, prefix() + "tked_time", params);
}

void
OpenPronghornLinearFVTurbulencePhysics::addKEpsilonAdvection()
{
  const std::string kernel_type = "LinearFVTurbulentAdvection";
  InputParameters params = getFactory().getValidParams(kernel_type);
  assignBlocks(params, _blocks);
  params.set<UserObjectName>("rhie_chow_user_object") = _flow_equations_physics->rhieChowUOName();
  params.set<MooseEnum>("advected_interp_method") =
      getParam<MooseEnum>("tke_advection_interpolation");

  params.set<LinearVariableName>("variable") = _tke_name;
  getProblem().addLinearFVKernel(kernel_type, prefix() + "tke_advection", params);

  params.set<LinearVariableName>("variable") = _tked_name;
  params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
  params.set<MooseEnum>("advected_interp_method") =
      getParam<MooseEnum>("tked_advection_interpolation");
  getProblem().addLinearFVKernel(kernel_type, prefix() + "tked_advection", params);
}

void
OpenPronghornLinearFVTurbulencePhysics::addKEpsilonDiffusion()
{
  // Single diffusion kernel per equation avoids double application of flux BCs
  const std::string kernel_type = "LinearFVTurbulentDiffusion";
  InputParameters params = getFactory().getValidParams(kernel_type);
  assignBlocks(params, _blocks);
  params.set<bool>("use_nonorthogonal_correction") =
      getParam<bool>("use_nonorthogonal_correction");

  params.set<LinearVariableName>("variable") = _tke_name;
  params.set<MooseFunctorName>("diffusion_coeff") = "mu_eff_tke";
  getProblem().addLinearFVKernel(kernel_type, prefix() + "tke_diffusion", params);

  params.set<LinearVariableName>("variable") = _tked_name;
  params.set<MooseFunctorName>("diffusion_coeff") = "mu_eff_tked";
  params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
  getProblem().addLinearFVKernel(kernel_type, prefix() + "tked_diffusion", params);
}

void
OpenPronghornLinearFVTurbulencePhysics::addKEpsilonSourceSink()
{
  const std::string u_names[3] = {"u", "v", "w"};

  const bool use_buoyancy = getParam<bool>("use_buoyancy");
  const bool use_compressibility = getParam<bool>("use_compressibility");

  // Shared lambda to set velocity and optional correction functors
  auto setVelocityAndCorrections = [&](InputParameters & p)
  {
    for (const auto d : make_range(dimension()))
      p.set<MooseFunctorName>(u_names[d]) = _velocity_names[d];
    if (use_buoyancy)
    {
      p.set<RealVectorValue>("gravity") = getParam<RealVectorValue>("gravity");
      if (isParamValid("temperature"))
        p.set<MooseFunctorName>("temperature") = getParam<MooseFunctorName>("temperature");
      if (isParamValid("beta"))
        p.set<MooseFunctorName>("beta") = getParam<MooseFunctorName>("beta");
    }
    if (use_compressibility && isParamValid("speed_of_sound"))
      p.set<MooseFunctorName>("speed_of_sound") = getParam<MooseFunctorName>("speed_of_sound");
  };

  // TKE source/sink (OpenPronghorn custom kernel)
  {
    const std::string kernel_type = "kEpsilonTKESourceSink";
    InputParameters params = getFactory().getValidParams(kernel_type);
    assignBlocks(params, _blocks);
    params.set<LinearVariableName>("variable") = _tke_name;
    params.set<MooseFunctorName>(NS::TKED) = _tked_name;
    params.set<MooseFunctorName>(NS::density) = _flow_equations_physics->densityName();
    params.set<MooseFunctorName>(NS::mu) = _flow_equations_physics->dynamicViscosityName();
    params.set<MooseFunctorName>(NS::mu_t) = _turbulent_viscosity_name;
    params.set<Real>("C_mu") = getParam<Real>("C_mu");
    params.set<Real>("C_pl") = getParam<Real>("C_pl");
    params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
    params.set<MooseEnum>("wall_treatment") = _wall_treatment_eps;
    params.set<MooseEnum>("k_epsilon_variant") = getParam<MooseEnum>("k_epsilon_variant");
    params.set<MooseEnum>("nonlinear_model") = getParam<MooseEnum>("nonlinear_model");
    params.set<MooseEnum>("curvature_model") = getParam<MooseEnum>("curvature_model");
    params.set<bool>("use_buoyancy") = use_buoyancy;
    params.set<bool>("use_compressibility") = use_compressibility;
    params.set<bool>("use_curvature_correction") = getParam<bool>("use_curvature_correction");
    setVelocityAndCorrections(params);
    getProblem().addLinearFVKernel(kernel_type, prefix() + "tke_source_sink", params);
  }

  // TKED (epsilon) source/sink (OpenPronghorn custom kernel)
  // C1_eps and C2_eps are stored as MooseFunctorName in the base class but the kernel
  // takes Real, so we convert via std::stod (they are always constant strings e.g. "1.44").
  {
    const std::string kernel_type = "kEpsilonTKEDSourceSink";
    InputParameters params = getFactory().getValidParams(kernel_type);
    assignBlocks(params, _blocks);
    params.set<LinearVariableName>("variable") = _tked_name;
    params.set<MooseFunctorName>(NS::TKE) = _tke_name;
    params.set<MooseFunctorName>(NS::density) = _flow_equations_physics->densityName();
    params.set<MooseFunctorName>(NS::mu) = _flow_equations_physics->dynamicViscosityName();
    params.set<MooseFunctorName>(NS::mu_t) = _turbulent_viscosity_name;
    params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
    params.set<MooseEnum>("wall_treatment") = _wall_treatment_eps;
    params.set<Real>("C_mu") = getParam<Real>("C_mu");
    params.set<Real>("C_pl") = getParam<Real>("C_pl");
    params.set<Real>("C1_eps") =
        std::stod(std::string(getParam<MooseFunctorName>("C1_eps")));
    params.set<Real>("C2_eps") =
        std::stod(std::string(getParam<MooseFunctorName>("C2_eps")));
    params.set<Real>("C3_eps") = getParam<Real>("C3_eps");
    params.set<Real>("Ct") = getParam<Real>("Ct");
    params.set<Real>("C_lowRe") = getParam<Real>("C_lowRe");
    params.set<Real>("D_lowRe") = getParam<Real>("D_lowRe");
    params.set<Real>("E_lowRe") = getParam<Real>("E_lowRe");
    params.set<Real>("Cw") = getParam<Real>("Cw");
    params.set<MooseEnum>("k_epsilon_variant") = getParam<MooseEnum>("k_epsilon_variant");
    params.set<MooseEnum>("nonlinear_model") = getParam<MooseEnum>("nonlinear_model");
    params.set<MooseEnum>("curvature_model") = getParam<MooseEnum>("curvature_model");
    params.set<bool>("use_yap") = getParam<bool>("use_yap");
    params.set<bool>("use_low_re_Gprime") = getParam<bool>("use_low_re_Gprime");
    params.set<bool>("use_buoyancy") = use_buoyancy;
    params.set<bool>("use_compressibility") = use_compressibility;
    params.set<bool>("use_curvature_correction") = getParam<bool>("use_curvature_correction");
    params.set<MooseFunctorName>("wall_distance") = _wall_distance_name;
    setVelocityAndCorrections(params);
    getProblem().addLinearFVKernel(kernel_type, prefix() + "tked_source_sink", params);
  }
}

// ---------------------------------------------------------------------------
// Boundary conditions
// ---------------------------------------------------------------------------

void
OpenPronghornLinearFVTurbulencePhysics::addFVBCs()
{
  if (_turbulence_model != "k-epsilon" || !getParam<bool>("mu_t_as_aux_variable"))
    return;

  const std::string u_names[3] = {"u", "v", "w"};
  const std::string bc_type = "LinearFVTurbulentViscosityWallFunctionBC";

  InputParameters params = getFactory().getValidParams(bc_type);
  params.set<std::vector<BoundaryName>>("boundary") = _turbulence_walls;
  params.set<LinearVariableName>("variable") = _turbulent_viscosity_name;
  params.set<MooseFunctorName>(NS::density) = _flow_equations_physics->densityName();
  params.set<MooseFunctorName>(NS::mu) = _flow_equations_physics->dynamicViscosityName();
  params.set<MooseFunctorName>(NS::TKE) = _tke_name;
  params.set<Real>("C_mu") = getParam<Real>("C_mu");
  params.set<MooseEnum>("wall_treatment") = _wall_treatment_eps;
  for (const auto d : make_range(dimension()))
    params.set<MooseFunctorName>(u_names[d]) = _velocity_names[d];

  getProblem().addLinearFVBC(bc_type, prefix() + "turbulence_walls", params);

  if (getParam<MultiMooseEnum>("turbulence_inlet_types").isValid() && _flow_equations_physics)
    addTurbulenceInletBCs();

  if (getParam<bool>("add_turbulence_outflow_bcs") && _flow_equations_physics)
    addTurbulenceOutletBCs();
}

void
OpenPronghornLinearFVTurbulencePhysics::addTurbulenceInletBCs()
{
  const auto & inlet_boundaries = _flow_equations_physics->getInletBoundaries();
  const auto & inlet_types = getParam<MultiMooseEnum>("turbulence_inlet_types");

  if (inlet_types.size() != inlet_boundaries.size())
    paramError("turbulence_inlet_types",
               "Expected one entry per flow-physics inlet boundary (" +
                   std::to_string(inlet_boundaries.size()) + " found), got " +
                   std::to_string(inlet_types.size()) + ".");

  // Count each type to validate companion parameter sizes
  unsigned int ils_count = 0, fv_count = 0;
  for (const auto i : index_range(inlet_types))
  {
    if (inlet_types[i] == "intensity-length-scale")
      ++ils_count;
    else
      ++fv_count;
  }

  if (ils_count > 0)
  {
    if (!isParamValid("turbulence_intensity") || !isParamValid("turbulence_length_scale") ||
        !isParamValid("turbulence_reference_velocity"))
      paramError("turbulence_inlet_types",
                 "Parameters turbulence_intensity, turbulence_length_scale, and "
                 "turbulence_reference_velocity must all be provided for the "
                 "'intensity-length-scale' inlet type.");

    const auto & I_vals = getParam<std::vector<Real>>("turbulence_intensity");
    const auto & L_vals = getParam<std::vector<Real>>("turbulence_length_scale");
    const auto & U_vals = getParam<std::vector<Real>>("turbulence_reference_velocity");
    if (I_vals.size() != ils_count)
      paramError("turbulence_intensity",
                 "Expected " + std::to_string(ils_count) +
                     " value(s) (one per intensity-length-scale inlet), got " +
                     std::to_string(I_vals.size()) + ".");
    if (L_vals.size() != ils_count)
      paramError("turbulence_length_scale",
                 "Expected " + std::to_string(ils_count) +
                     " value(s) (one per intensity-length-scale inlet), got " +
                     std::to_string(L_vals.size()) + ".");
    if (U_vals.size() != ils_count)
      paramError("turbulence_reference_velocity",
                 "Expected " + std::to_string(ils_count) +
                     " value(s) (one per intensity-length-scale inlet), got " +
                     std::to_string(U_vals.size()) + ".");
  }

  if (fv_count > 0)
  {
    if (!isParamValid("tke_inlet_functors") || !isParamValid("tked_inlet_functors"))
      paramError("turbulence_inlet_types",
                 "Parameters tke_inlet_functors and tked_inlet_functors must be provided "
                 "for the 'fixed-value' inlet type.");
    const auto & tke_f = getParam<std::vector<MooseFunctorName>>("tke_inlet_functors");
    const auto & eps_f = getParam<std::vector<MooseFunctorName>>("tked_inlet_functors");
    if (tke_f.size() != fv_count)
      paramError("tke_inlet_functors",
                 "Expected " + std::to_string(fv_count) +
                     " value(s) (one per fixed-value inlet), got " +
                     std::to_string(tke_f.size()) + ".");
    if (eps_f.size() != fv_count)
      paramError("tked_inlet_functors",
                 "Expected " + std::to_string(fv_count) +
                     " value(s) (one per fixed-value inlet), got " +
                     std::to_string(eps_f.size()) + ".");
  }

  // Retrieve params once (guard avoids calling getParam on unset params)
  std::vector<Real> I_vals, L_vals, U_vals;
  std::vector<MooseFunctorName> tke_functors_fv, eps_functors_fv;
  if (ils_count > 0)
  {
    I_vals = getParam<std::vector<Real>>("turbulence_intensity");
    L_vals = getParam<std::vector<Real>>("turbulence_length_scale");
    U_vals = getParam<std::vector<Real>>("turbulence_reference_velocity");
  }
  if (fv_count > 0)
  {
    tke_functors_fv = getParam<std::vector<MooseFunctorName>>("tke_inlet_functors");
    eps_functors_fv = getParam<std::vector<MooseFunctorName>>("tked_inlet_functors");
  }

  const Real c_mu = getParam<Real>("C_mu");
  unsigned int ils_i = 0, fv_i = 0;

  for (const auto i : index_range(inlet_boundaries))
  {
    const auto & bdy = inlet_boundaries[i];
    MooseFunctorName tke_functor, eps_functor;

    if (inlet_types[i] == "intensity-length-scale")
    {
      const Real k = 1.5 * I_vals[ils_i] * I_vals[ils_i] * U_vals[ils_i] * U_vals[ils_i];
      const Real eps = std::pow(c_mu, 0.75) * std::pow(k, 1.5) / L_vals[ils_i];
      // Represent computed constants as high-precision strings for the MOOSE functor system
      std::ostringstream oss_k, oss_e;
      oss_k << std::setprecision(15) << k;
      oss_e << std::setprecision(15) << eps;
      tke_functor = oss_k.str();
      eps_functor = oss_e.str();
      ++ils_i;
    }
    else
    {
      tke_functor = tke_functors_fv[fv_i];
      eps_functor = eps_functors_fv[fv_i];
      ++fv_i;
    }

    const std::string bc_type = "LinearFVAdvectionDiffusionFunctorDirichletBC";
    {
      auto params = getFactory().getValidParams(bc_type);
      params.set<std::vector<BoundaryName>>("boundary") = {bdy};
      params.set<LinearVariableName>("variable") = _tke_name;
      params.set<MooseFunctorName>("functor") = tke_functor;
      getProblem().addLinearFVBC(bc_type, prefix() + "tke_inlet_" + bdy, params);
    }
    {
      auto params = getFactory().getValidParams(bc_type);
      params.set<std::vector<BoundaryName>>("boundary") = {bdy};
      params.set<LinearVariableName>("variable") = _tked_name;
      params.set<MooseFunctorName>("functor") = eps_functor;
      getProblem().addLinearFVBC(bc_type, prefix() + "tked_inlet_" + bdy, params);
    }
  }
}

void
OpenPronghornLinearFVTurbulencePhysics::addTurbulenceOutletBCs()
{
  const auto & outlet_boundaries = _flow_equations_physics->getOutletBoundaries();
  const std::string bc_type = "LinearFVAdvectionDiffusionOutflowBC";

  for (const auto & bdy : outlet_boundaries)
  {
    {
      auto params = getFactory().getValidParams(bc_type);
      params.set<std::vector<BoundaryName>>("boundary") = {bdy};
      params.set<LinearVariableName>("variable") = _tke_name;
      params.set<bool>("use_two_term_expansion") = false;
      getProblem().addLinearFVBC(bc_type, prefix() + "tke_outlet_" + bdy, params);
    }
    {
      auto params = getFactory().getValidParams(bc_type);
      params.set<std::vector<BoundaryName>>("boundary") = {bdy};
      params.set<LinearVariableName>("variable") = _tked_name;
      params.set<bool>("use_two_term_expansion") = false;
      getProblem().addLinearFVBC(bc_type, prefix() + "tked_outlet_" + bdy, params);
    }
  }
}

// ---------------------------------------------------------------------------
// Functor materials
// ---------------------------------------------------------------------------

void
OpenPronghornLinearFVTurbulencePhysics::addFunctorMaterials()
{
  if (_turbulence_model != "k-epsilon")
    return;

  const std::string mat_type = "FunctorEffectiveDynamicViscosity";
  const auto mu_name = _flow_equations_physics->dynamicViscosityName();

  // mu_eff = mu + mu_t  (used by WCNSLinearFVFlowPhysics for momentum flux)
  if (!getProblem().hasFunctor(NS::mu_eff, /*thread_id=*/0))
  {
    InputParameters params = getFactory().getValidParams(mat_type);
    assignBlocks(params, _blocks);
    params.set<MooseFunctorName>("property_name") = NS::mu_eff;
    params.set<MooseFunctorName>(NS::mu) = mu_name;
    params.set<MooseFunctorName>(NS::mu_t) = _turbulent_viscosity_name;
    params.set<MooseFunctorName>(NS::mu_t + "_inverse_factor") = "1";
    getProblem().addMaterial(mat_type, prefix() + "effective_viscosity", params);
  }

  // mu_eff_tke = mu + mu_t / sigma_k  (TKE diffusion coefficient)
  {
    InputParameters params = getFactory().getValidParams(mat_type);
    assignBlocks(params, _blocks);
    params.set<MooseFunctorName>("property_name") = "mu_eff_tke";
    params.set<MooseFunctorName>(NS::mu) = mu_name;
    params.set<MooseFunctorName>(NS::mu_t) = _turbulent_viscosity_name;
    params.set<MooseFunctorName>(NS::mu_t + "_inverse_factor") =
        getParam<MooseFunctorName>("sigma_k");
    getProblem().addMaterial(mat_type, prefix() + "mu_eff_tke", params);
  }

  // mu_eff_tked = mu + mu_t / sigma_eps  (epsilon diffusion coefficient)
  {
    InputParameters params = getFactory().getValidParams(mat_type);
    assignBlocks(params, _blocks);
    params.set<MooseFunctorName>("property_name") = "mu_eff_tked";
    params.set<MooseFunctorName>(NS::mu) = mu_name;
    params.set<MooseFunctorName>(NS::mu_t) = _turbulent_viscosity_name;
    params.set<MooseFunctorName>(NS::mu_t + "_inverse_factor") =
        getParam<MooseFunctorName>("sigma_eps");
    getProblem().addMaterial(mat_type, prefix() + "mu_eff_tked", params);
  }
}
