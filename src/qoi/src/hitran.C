
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

// this class
#include "grins/hitran.h"

// libMesh

// GRINS
#include "grins/string_utils.h"

// C++
#include <fstream>
#include <algorithm>

namespace GRINS
{
  HITRAN::HITRAN(const std::vector<std::string> & data_files, const std::vector<std::string> & partition_function_files,
                 libMesh::Real T_min, libMesh::Real T_max, libMesh::Real T_step)
  : _Tmin(T_min),
    _Tmax(T_max),
    _Tstep(T_step),
    _T0(296.0)
  {
    // sanity checks on temperature range specification
    if ( (T_min<0.0) || (T_min>=T_max) || (T_step<=0.0) || (T_min>_T0) || (T_max<_T0) )
      {
        std::stringstream ss;
        ss <<"Invalid specification of temperature range:" <<std::endl;
        ss <<"T_min: " <<T_min <<std::endl;
        ss <<"T_max: " <<T_max <<std::endl;
        ss <<"T_step: " <<T_step <<std::endl;       
        libmesh_error_msg(ss.str());
      }
      
    _data_size.resize(data_files.size());
    _qT.resize(data_files.size());
    _qT0.resize(data_files.size());

    for (unsigned int s=0; s<data_files.size(); s++)
      {
        _isotop.push_back(std::vector<unsigned int>());
        _nu.push_back(std::vector<libMesh::Real>());
        _sw.push_back(std::vector<libMesh::Real>());
        _gamma_air.push_back(std::vector<libMesh::Real>());
        _gamma_self.push_back(std::vector<libMesh::Real>());
        _elower.push_back(std::vector<libMesh::Real>());
        _n.push_back(std::vector<libMesh::Real>());
        _delta_air.push_back(std::vector<libMesh::Real>());      
      
        // open data file
        std::ifstream hitran_file;
        hitran_file.open(data_files[s]);
        if (!hitran_file.is_open())
          {
            std::stringstream ss;
            ss <<"Unable to open provided hitran_file: " <<data_files[s] <<std::endl;    
            libmesh_error_msg(ss.str());
          }
        
        while(!hitran_file.eof())
          {
            std::string line;
            getline(hitran_file,line);
            
            if (line == "")
              continue;
            
            std::vector<libMesh::Real> vals;
            
            GRINS::StringUtilities::split_string_real(line,",",vals);
            
            libmesh_assert_equal_to(vals.size(),8);
            
            // HITRAN gives isotopologue numbers starting at 1
            // shift to zero-based to make indexing easier
            _isotop[s].push_back(static_cast<unsigned int>(vals[0]-1));
            
            _nu[s].push_back(vals[1]);
            _sw[s].push_back(vals[2]);
            _gamma_air[s].push_back(vals[3]);
            _gamma_self[s].push_back(vals[4]);
            _elower[s].push_back(vals[5]);
            _n[s].push_back(vals[6]);
            _delta_air[s].push_back(vals[7]);      
          }

        // sanity checks
        libmesh_assert_equal_to( _isotop[s].size(), _nu[s].size() );
        libmesh_assert_equal_to( _nu[s].size(), _sw[s].size() );
        libmesh_assert_equal_to( _sw[s].size(), _gamma_air[s].size() );
        libmesh_assert_equal_to( _gamma_air[s].size(), _gamma_self[s].size() );
        libmesh_assert_equal_to( _gamma_self[s].size(), _elower[s].size() );
        libmesh_assert_equal_to( _elower[s].size(), _n[s].size() );
        libmesh_assert_equal_to( _n[s].size(), _delta_air[s].size() );

        // save data length and close HITRAN data file
        _data_size[s] = _isotop[s].size(); // all data vectors are the same length
        hitran_file.close();

        // file with partition function values
        std::ifstream qT_file;
        qT_file.open(partition_function_files[s]);
        if (!qT_file.is_open())
          {
            std::stringstream ss;
            ss <<"Unable to open provided partition_function_file: " <<partition_function_files[s] <<std::endl;     
            libmesh_error_msg(ss.str());
          }
        
        // number of temperature values
        unsigned int num_T = (T_max-T_min)/T_step + 1;
        
        // read the partition function values
        unsigned int counter = 0;
        
        while(!qT_file.eof())
          {
            std::string line;
            getline(qT_file,line);
            
            if (line == "")
              continue;
            
            _qT[s].push_back(std::vector<libMesh::Real>());
            
            GRINS::StringUtilities::split_string_real(line,",",_qT[s][counter]);
            
            // we should have a partition function value for each temperature
            libmesh_assert_equal_to(num_T,_qT[s][counter].size());
            
            counter++;   
          }

        // save length and close partition sum file
        _q_size = num_T;
        qT_file.close();
        
        // cache the partition function values at the referece temperature
        for(unsigned int i=0; i<_qT[s].size(); i++)
          _qT0[s].push_back(this->search_partition_function(_T0,s,i));
      }

  }

  unsigned int HITRAN::get_data_size(unsigned int species_idx)
  {
    return _data_size[species_idx];
  }

  unsigned int HITRAN::isotopologue(unsigned int species_idx, unsigned int index)
  {
    libmesh_assert_less(index,_isotop[species_idx].size());
    return _isotop[species_idx][index];
  }

  libMesh::Real HITRAN::nu0(unsigned int species_idx, unsigned int index)
  {
    libmesh_assert_less(index,_nu[species_idx].size());
    return _nu[species_idx][index];
  }

  libMesh::Real HITRAN::sw(unsigned int species_idx, unsigned int index)
  {
    libmesh_assert_less(index,_sw[species_idx].size());
    return _sw[species_idx][index];
  }

  libMesh::Real HITRAN::gamma_air(unsigned int species_idx, unsigned int index)
  {
    libmesh_assert_less(index,_gamma_air[species_idx].size());
    return _gamma_air[species_idx][index];
  }

  libMesh::Real HITRAN::gamma_self(unsigned int species_idx, unsigned int index)
  {
    libmesh_assert_less(index,_gamma_self[species_idx].size());
    return _gamma_self[species_idx][index];
  }

  libMesh::Real HITRAN::elower(unsigned int species_idx, unsigned int index)
  {
    libmesh_assert_less(index,_elower[species_idx].size());
    return _elower[species_idx][index];
  }

  libMesh::Real HITRAN::n_air(unsigned int species_idx, unsigned int index)
  {
    libmesh_assert_less(index,_n[species_idx].size());
    return _n[species_idx][index];
  }

  libMesh::Real HITRAN::delta_air(unsigned int species_idx, unsigned int index)
  {
    libmesh_assert_less(index,_delta_air[species_idx].size());
    return _delta_air[species_idx][index];
  }

  libMesh::Real HITRAN::partition_function(libMesh::Real T, unsigned int species_idx, unsigned int iso)
  {
    libMesh::Real qt;

    if (T==_T0)
      qt = _qT0[species_idx][iso];
    else
      qt = this->search_partition_function(T,species_idx,iso);

    return qt;
  }

  libMesh::Real HITRAN::search_partition_function(libMesh::Real T, unsigned int species_idx, unsigned int iso)
  {
    libMesh::Real retval = -1.0;

    int i = _T_index(T);

    if (i >= 0)
        retval = this->interpolate_values(i,T,_qT[species_idx][iso]);
    else
      {
        std::stringstream ss;
        ss <<"Error: Temperature " <<T <<"K does not exist in the given partition sum data" <<std::endl;
        libmesh_error_msg(ss.str());
      }

    return retval;
  }
  
  unsigned int HITRAN::get_species_idx(unsigned int species_var_index)
  {
    unsigned int retval = libMesh::invalid_uint;

    std::map<unsigned int,unsigned int>::iterator it;
    it = _species_var_map.find(species_var_index);
    if (it != _species_var_map.end())
      retval = it->second;

    return retval;
  }

  void HITRAN::add_species(unsigned int species_var_index, unsigned int species_idx)
  {
    _species_var_map[species_var_index] = species_idx;
  }

  int HITRAN::_T_index(libMesh::Real T)
  {
    unsigned int index = std::ceil((T-_Tmin)/_Tstep) + 1;
    return index;
  }
  
  libMesh::Real HITRAN::interpolate_values( int index_r, libMesh::Real T_star, const std::vector<libMesh::Real> & y) const
  {
    if ( (T_star>_Tmax) || (T_star<_Tmin) )
      {
        std::stringstream ss;
        ss <<"Error, temperature T=" <<T_star <<" is outside the specified range of provided partition function values" <<std::endl;
        libmesh_error_msg(ss.str());
      }
      
    libMesh::Real T = _Tmin + (_Tstep*index_r);
    return y[index_r-1] + ( (T_star-(T-_Tstep))*(y[index_r]-y[index_r-1]) )/(T-(T-_Tstep));
  }

}
