//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

// GRINS
#include "grins/runner.h"
#include "grins/string_utils.h"
#include "grins/absorption_coeff.h"
#include "grins/simulation.h"
#include "grins/multiphysics_sys.h"
#include "grins/composite_qoi.h"
#include "grins/spectroscopic_absorption.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

int main(int argc, char* argv[])
{
#if GRINS_HAVE_ANTIOCH
  GRINS::Runner runner(argc,argv);
  runner.init();
  
  const GetPot & input = runner.get_input_file();
  
  unsigned int gq_order = input.vector_variable_size("SpectroscopyExample/weights");
  
  std::vector<libMesh::Real> weights(gq_order);
  for (unsigned int i=0; i<gq_order; ++i)
    weights[i] = input("SpectroscopyExample/weights",0.0,i);
    
  std::vector<libMesh::Real> intensity_val(gq_order);
  for (unsigned int i=0; i<gq_order; ++i)
    {
      intensity_val[i] = input("SpectroscopyExample/intensity",1.0,i);
    }
  
  libMesh::Real a = input("SpectroscopyExample/a",0.0);
  libMesh::Real b = input("SpectroscopyExample/b",0.0);
  
  libMesh::Real I0_star = 0.0;
  for (unsigned int i=0; i<weights.size(); ++i)
    I0_star += weights[i]*intensity_val[i];
    
  I0_star *= (b-a)/2.0;
  
  runner.run(); // do the rayfire runs

  GRINS::Simulation & sim = runner.get_simulation();
  GRINS::MultiphysicsSystem * system = sim.get_multiphysics_system();
  
  std::string filename = "qoi_vals";

//  std::ofstream output;
//  if (system->get_mesh().comm().rank() == 0)
//    {
//      output.open(filename+".dat",std::ofstream::app);
//    }

  libMesh::Real If_star = 0.0;

  for (unsigned int i=0; i<weights.size(); ++i)
    {
      libMesh::Real val = sim.get_qoi_value(i);
      
      libMesh::Real contrib = val*weights[i]*intensity_val[i];
      
      If_star += contrib;
      
//      if (system->get_mesh().comm().rank() == 0)
//        output <<std::setprecision(16) <<val <<"," <<weights[i] <<"," <<intensity_val[i] <<std::endl;
    }

  libMesh::Real Q_star = 1.0 - ((b-a)/2.0)*(If_star/I0_star);

//  if (system->get_mesh().comm().rank() == 0)
//    output.close();

  if (system->get_mesh().comm().rank() == 0)
        std::cout <<"\n\n===============================\n" <<"total absorption: " <<std::setprecision(16) <<Q_star <<"\n===============================" <<std::endl;

#else
  libmesh_error_msg("ERROR: GRINS must be built with Antioch to use the Spectroscopy example. Please reconfigure your build to include the Antioch library.");
#endif
  return 0;
}

