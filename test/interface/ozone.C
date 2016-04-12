//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#include "grins_config.h"

#ifdef GRINS_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include "grins_test_paths.h"
#include "nasa_thermo_test_base.h"
#include "kinetics_test_base.h"
#include "antioch_test_base.h"
#include "cantera_test_base.h"

// GRINS
#include "grins/antioch_evaluator.h"
#include "grins/cantera_evaluator.h"

namespace GRINSTesting
{
  class OzoneTestBase
  {
  public:

    void init_ozone_test( unsigned int& O_idx, unsigned int& O2_idx,
                        unsigned int& O3_idx,
                        std::vector<unsigned int>& active_species )
    {
      O_idx = 0;
      O2_idx = 1;
      O3_idx  = 2;

      active_species.resize(3);
      active_species[0] = O_idx;
      active_species[1] = O2_idx;
      active_species[2] = O3_idx;
    }

  };

  class OzoneNASA9Thermo : public OzoneTestBase,
                         public NASA9ThermoTestBase
  {
  public:

    OzoneNASA9Thermo()
    {
      this->init_ozone_test( _O_idx, _O2_idx, _O3_idx,
                           _active_species );
    }

  };

  class OzoneNASA7Thermo : public OzoneTestBase,
                           public NASA7ThermoTestBase
  {
  public:

    OzoneNASA7Thermo()
    {
      this->init_ozone_test( _O_idx, _O2_idx, _O3_idx,
                             _active_species );
    }

  };
  
  class OzoneKineticsTestBase : public KineticsTestBase,
                              public OzoneTestBase
  {
  public:

    void init_ozone_kinetics()
    {
      this->init_ozone_test( _O_idx, _O2_idx, _O3_idx,
                           _active_species );

      _n_reactions = 3;

      unsigned int n_species = _active_species.size();

      // All in kJ/mol
      _Ea_coeffs.resize(_n_reactions);
      _Ea_coeffs[0] = 0.0;
      _Ea_coeffs[1] = -4.234;
      _Ea_coeffs[2] = 17.38;

      // Convert to J/kmol (since this is what GRINS::Constants::R_universal is in)
      for( unsigned int r = 0; r < _n_reactions; r++ )
        _Ea_coeffs[r] *= 1000.0*1000.0;

      // m^3/kmol
      _preexp_coeffs.resize(_n_reactions);
      _preexp_coeffs[0] = 2.9e+11;
      _preexp_coeffs[1] = 3.427e+7;
      _preexp_coeffs[2] = 5.2e+9;

      _temp_exp_coeffs.resize(_n_reactions);
      _temp_exp_coeffs[0] = -1.0;
      _temp_exp_coeffs[1] = 0.0;
      _temp_exp_coeffs[2] = 0.0;

      _three_body_coeffs.resize(_n_reactions);
      _is_three_body_rxn.resize(_n_reactions,false);
      _is_three_body_rxn[0] = true;
      _three_body_coeffs[0].resize(n_species, 1.0);
      _three_body_coeffs[0][_O_idx] = 1.13;
      _three_body_coeffs[0][_O2_idx] = 0.94;
      _three_body_coeffs[0][_O3_idx] = 0.92;
      
      _is_three_body_rxn[1] = true;
      _three_body_coeffs[1].resize(n_species, 1.0);
      _three_body_coeffs[1][_O_idx] = 1.13;
      _three_body_coeffs[1][_O2_idx] = 0.94;
      _three_body_coeffs[1][_O3_idx] = 0.92;

      _reactant_stoich_coeffs.resize(_n_reactions);
      _product_stoich_coeffs.resize(_n_reactions);

      // O + O + M <=> O2 + M
      _reactant_stoich_coeffs[0].resize(n_species, 0.0);
      _product_stoich_coeffs[0].resize(n_species, 0.0);
      _reactant_stoich_coeffs[0][_O_idx] = 2.0;
      _product_stoich_coeffs[0][_O2_idx] = 1.0;

      // O2 + O + M <=> O3 + M
      _reactant_stoich_coeffs[1].resize(n_species, 0.0);
      _product_stoich_coeffs[1].resize(n_species, 0.0);
      _reactant_stoich_coeffs[1][_O2_idx] = 1.0;
      _product_stoich_coeffs[1][_O_idx] = 1.0;
      _product_stoich_coeffs[1][_O3_idx] = 1.0;

      // O + O3 <=> O2 + O2
      _reactant_stoich_coeffs[2].resize(n_species, 0.0);
      _product_stoich_coeffs[2].resize(n_species, 0.0);
      _reactant_stoich_coeffs[2][_O_idx] = 1.0;
      _reactant_stoich_coeffs[2][_O3_idx] = 1.0;
      _product_stoich_coeffs[2][_O2_idx] = 2.0;
    }

  };


#ifdef GRINS_HAVE_ANTIOCH

  class AntiochOzoneNASA9ThermoTest : public CppUnit::TestCase,
                                    public AntiochTestBase,
                                    public OzoneTestBase,
                                    public NASA9ThermoTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochOzoneNASA9ThermoTest );

    CPPUNIT_TEST( test_cp );
    CPPUNIT_TEST( test_hs );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/antioch_ozone.in";
      this->init_antioch(input_file, "TestMaterial");

      this->init_ozone_test( _O_idx, _O2_idx, _O3_idx,
                           _active_species );

      //this->check_indices(*_antioch_mixture);
    }

    void test_cp()
    {
      std::vector<libMesh::Real> Y(3);
      Y[_O_idx] = 0.0;
      Y[_O2_idx] = 0.8;
      Y[_O3_idx] = 0.2;

      this->test_cp_common<GRINS::AntiochMixture,GRINS::AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> > >
        ( *_antioch_mixture, Y, TestingUtils::epsilon()*100 );
    }

    void test_hs()
    {
      this->test_h_common<GRINS::AntiochMixture,GRINS::AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> > >
        ( *_antioch_mixture, TestingUtils::epsilon()*100 );
    }
  };

  class AntiochOzoneNASA9KineticsTest : public CppUnit::TestCase,
                                      public AntiochTestBase,
                                      public OzoneKineticsTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochOzoneNASA9KineticsTest );

    CPPUNIT_TEST( test_omega_dot );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/antioch_ozone.in";
      this->init_antioch(input_file, "TestMaterial");

      this->init_ozone_kinetics();

      //this->check_indices(*_antioch_mixture);
    }

    void test_omega_dot()
    {
      std::vector<libMesh::Real> Y(3);
      Y[_O_idx] = 0.25;
      Y[_O2_idx] = 0.25;
      Y[_O3_idx] = 0.50;

      OzoneNASA9Thermo thermo;

      this->test_omega_dot_common<GRINS::AntiochMixture,GRINS::AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> > >
        ( *_antioch_mixture, thermo, Y, TestingUtils::epsilon()*1e3 );
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochOzoneNASA9ThermoTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochOzoneNASA9KineticsTest );

#endif // GRINS_HAVE_ANTIOCH


#ifdef GRINS_HAVE_CANTERA

  class CanteraOzoneNASA9ThermoTest : public CppUnit::TestCase,
                                    public CanteraTestBase,
                                    public OzoneTestBase,
                                    public NASA9ThermoTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( CanteraOzoneNASA9ThermoTest );

    CPPUNIT_TEST( test_cp );
    CPPUNIT_TEST( test_hs );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/cantera_ozone.in";
      this->init_cantera(input_file, "TestMaterial");

      this->init_ozone_test( _O_idx, _O2_idx, _O3_idx,
                             _active_species );

      //this->check_indices(*_antioch_mixture);
    }

    void test_cp()
    {
      std::vector<libMesh::Real> Y(3);
      Y[_O_idx] = 1.0;
      Y[_O2_idx] = 0.0;
      Y[_O3_idx] = 0.0;

      this->test_cp_common<GRINS::CanteraMixture,GRINS::CanteraEvaluator>
        ( *_cantera_mixture, Y, 1.0e-4 );
    }

    void test_hs()
    {
      this->test_h_common<GRINS::CanteraMixture,GRINS::CanteraEvaluator>
        ( *_cantera_mixture, 1.0e-4 );
    }
  };

  class CanteraOzoneNASA9KineticsTest : public CppUnit::TestCase,
                                      public CanteraTestBase,
                                      public OzoneKineticsTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( CanteraOzoneNASA9KineticsTest );

    CPPUNIT_TEST( test_omega_dot );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/cantera_ozone.in";
      this->init_cantera(input_file, "TestMaterial");

      this->init_ozone_kinetics();

      //this->check_indices(*_antioch_mixture);
    }

    void test_omega_dot()
    {
      std::cout << "\n******************************\n\nRunning Cantera test_omega_dot()" << std::endl;

      std::vector<libMesh::Real> Y(3);
      Y[_O_idx] = 1.0;
      Y[_O2_idx] = 0.0;
      Y[_O3_idx] = 0.0;

      OzoneNASA9Thermo thermo;

      this->test_omega_dot_common<GRINS::CanteraMixture,GRINS::CanteraEvaluator>
        ( *_cantera_mixture, thermo, Y, 3.0e-1 );
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( CanteraOzoneNASA9ThermoTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( CanteraOzoneNASA9KineticsTest );

#endif // GRINS_HAVE_CANTERA

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
