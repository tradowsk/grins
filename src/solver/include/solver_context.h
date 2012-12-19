//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_SOLVER_CONTEXT_H
#define GRINS_SOLVER_CONTEXT_H

#include "boost/tr1/memory.hpp"

// libMesh
#include "equation_systems.h"

// GRINS
#include "multiphysics_sys.h"
#include "visualization.h"

namespace GRINS
{
  //! Simple class to hold objects passed to Solver::solve
  /*! Allows some flexibility for adding objects needed to pass to the Solver::solve
      method so that the solver can still be agnostic to creation etc. of those objects,
      but can operate on them. 
   */
  class SolverContext
  {
  public:
    
    SolverContext();
    ~SolverContext();

    GRINS::MultiphysicsSystem* system;
    std::tr1::shared_ptr<libMesh::EquationSystems> equation_system;
    std::tr1::shared_ptr<GRINS::Visualization> vis;
    bool output_vis;
    bool output_residual;

  };
}

#endif // GRINS_SOLVER_CONTEXT_H
