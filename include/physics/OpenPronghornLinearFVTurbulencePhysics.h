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

#pragma once

#include "WCNSFVTurbulencePhysicsBase.h"

/**
 * Physics action that sets up an OpenPronghorn k-epsilon turbulence model for a
 * linearly-segregated (SIMPLE) finite-volume flow simulation.
 *
 * Extends the MOOSE NavierStokes turbulence framework by substituting OpenPronghorn's
 * custom kernels (kEpsilonTKESourceSink, kEpsilonTKEDSourceSink, kEpsilonViscosity)
 * in place of the standard MOOSE equivalents, enabling k-epsilon model variants
 * (Standard, LowRe, TwoLayer, Realizable) and optional corrections (Yap, buoyancy,
 * compressibility, curvature).
 *
 * Usage: pair with WCNSLinearFVFlowPhysics under Physics/NavierStokes/FlowSegregated,
 * and place this physics under Physics/NavierStokes/TurbulenceSegregated with
 * type = OpenPronghornLinearFVTurbulencePhysics.
 */
class OpenPronghornLinearFVTurbulencePhysics final : public WCNSFVTurbulencePhysicsBase
{
public:
  static InputParameters validParams();

  OpenPronghornLinearFVTurbulencePhysics(const InputParameters & parameters);

protected:
  virtual void initializePhysicsAdditional() override;
  virtual void checkIntegrity() const override;
  unsigned short getNumberAlgebraicGhostingLayersNeeded() const override;

private:
  virtual void addSolverVariables() override;
  virtual void addAuxiliaryVariables() override;
  virtual void addAuxiliaryKernels() override;
  virtual void addFVKernels() override;
  virtual void addFVBCs() override;
  // Overridden to no-op; mu_eff materials are added in addFunctorMaterials() so the
  // linear FV (non-AD) material type is used regardless of the base class dynamic_cast.
  virtual void addMaterials() override {}
  virtual void addFunctorMaterials() override;

  void addKEpsilonTimeDerivatives();
  void addKEpsilonAdvection();
  void addKEpsilonDiffusion();
  void addKEpsilonSourceSink();
  void addTurbulenceInletBCs();
  void addTurbulenceOutletBCs();

  /// Name of the wall-distance auxiliary variable
  const VariableName _wall_distance_name;
};
