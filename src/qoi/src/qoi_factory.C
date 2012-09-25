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

#include "qoi_factory.h"

namespace GRINS
{
  QoIFactory::QoIFactory()
  {
    return;
  }
  
  QoIFactory::~QoIFactory()
  {
    return;
  }

  std::tr1::shared_ptr<QoIBase> QoIFactory::build(const GetPot& input)
  {
    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
    std::string qoi_name = input("QoI/enabled_qois", "none" );

    std::tr1::shared_ptr<QoIBase> qoi;
    
    this->add_qoi( input, qoi_name, qoi );

    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
    this->check_qoi_physics_consistency( input, qoi_name );

    if( input( "screen-options/echo_qoi", true ) )
      {
	/*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
	this->echo_qoi_list( qoi_name );
      }

    return qoi;
  }

  void QoIFactory::add_qoi( const GetPot& input, const std::string& qoi_name, std::tr1::shared_ptr<QoIBase>& qoi )
  {
    if( qoi_name == avg_nusselt )
      qoi.reset( new AverageNusseltNumber( input ) );
    else
      {
	 libMesh::err << "Error: Invalid QoI name " << qoi_name << std::endl;
	 libmesh_error();
      }

    return;
  }

  void QoIFactory::check_qoi_physics_consistency( const GetPot& input, 
						  const std::string& qoi_name )
  {
    int num_physics =  input.vector_variable_size("Physics/enabled_physics");

    // This should be checked other places, but let's be double sure.
    libmesh_assert(num_physics > 1);
  
    std::set<std::string> requested_physics;
    
    // Build Physics name set
    for( int i = 0; i < num_physics; i++ )
      {
	requested_physics.insert( input("Physics/enabled_physics", "NULL", i ) );
      }
  
    /* If it's Nusselt, we'd better have HeatTransfer. HeatTransfer implicitly requires fluids,
       so no need to check for those. `*/
    if( qoi_name == avg_nusselt )
      {
	if( requested_physics.find( heat_transfer ) == requested_physics.end() )
	  this->consistency_error_msg( qoi_name, heat_transfer );
      }
      
    return;
  }

  void QoIFactory::echo_qoi_list( const std::string& qoi_name )
  {
    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
    std::cout << "==========================================================" << std::endl
	      << "List of Enabled QoIs:" << std::endl
	      << qoi_name << std::endl
	      <<  "==========================================================" << std::endl;
    return;
  }

  void QoIFactory::consistency_error_msg( const std::string& qoi_name, const std::string& physics_name )
  {
    libMesh::err << "================================================================" << std::endl
	      << "Physics " << physics_name << " could not be found." << std::endl
	      << "It is required for QoI " << qoi_name << std::endl
	      << "================================================================" << std::endl;
    libmesh_error();
  }

} //namespace GRINS