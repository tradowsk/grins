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


// This class
#include "grins/absorption_coeff.h"

// GRINS
#include "grins/variable_warehouse.h"
#include "grins/physical_constants.h"
#include "grins/math_constants.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

#if GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

// libMesh
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace GRINS
{
  template<typename Chemistry>
  AbsorptionCoeff<Chemistry>::AbsorptionCoeff(SharedPtr<Chemistry> & chem, SharedPtr<HITRAN> & hitran,
                                              libMesh::Real nu_min, libMesh::Real nu_max,
                                              libMesh::Real desired_nu, const std::string & species,
                                              libMesh::Real thermo_pressure)
    : _chemistry(chem),
      _hitran(hitran),
      _nu(desired_nu),
      _T_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>("Temperature")),
      _P_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>("Pressure")),
      _Y_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>("SpeciesMassFractions")),
      _T0(296), // [K]
      _Pref(1), // [atm]
      _rad_coeff(1.43887752) // [cm K]
  {
    // sanity checks
    if ( (nu_min>nu_max) || (desired_nu>nu_max) || (desired_nu<nu_min) )
      {
        std::stringstream ss;
        ss <<"Invalid specification of wavenumber range:" <<std::endl
           <<"nu_min: " <<nu_min <<std::endl
           <<"nu_max: " <<nu_max <<std::endl
           <<"desired_nu: " <<desired_nu <<std::endl;
        libmesh_error_msg(ss.str());
      }

    _species_idx = _chemistry->species_index(species);
    unsigned int data_size = _hitran->get_data_size();

    bool min_flag = false;

    for (unsigned int i=0; i<data_size; i++)
      {
        if (_hitran->nu0(i) > nu_min)
          {
            _min_index = i;
            min_flag = true;
            break;
          }
      }

    if (!min_flag)
      {
        std::stringstream ss;
        ss <<"Minimum wavenumber " <<nu_min <<" is greater than the maximum wavenumber in provided HITRAN data";
        libmesh_error_msg(ss.str());
      }

    bool max_flag = false;

    for (unsigned int i=data_size-1; i>=0; i--)
      {
        if (_hitran->nu0(i) < nu_max)
          {
            _max_index = i;
            max_flag = true;
            break;
          }
      }

    if (!max_flag)
      _max_index = data_size-1;

    if (thermo_pressure == -1.0) {
      _calc_thermo_pressure = true;
      _thermo_pressure = -1.0; // not needed in this case
    } else {
      _calc_thermo_pressure = false;
      _thermo_pressure = thermo_pressure;
    }

    this->init_voigt();
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::operator()(const libMesh::FEMContext& context, const libMesh::Point& qp_xyz, const libMesh::Real /*t*/)
  {
    START_LOG("operator()","AbsorptionCoeff");
    libMesh::Real T,p,thermo_p; // temperature, hydrostatic pressure, thermodynamic pressure
    std::vector<libMesh::Real> Y(_chemistry->n_species()); // mass fractions

    context.point_value(_T_var.T(), qp_xyz, T); // [K]

    if (_calc_thermo_pressure) {
      libmesh_not_implemented();
    } else {
      thermo_p = _thermo_pressure;
    }

    context.point_value(_P_var.p(), qp_xyz, p); // [Pa]

    libMesh::Real P = p + thermo_p; // total pressure [Pa]
    libmesh_assert_greater(P,0.0);

    // all mass fractions needed to get M_mix
    for (unsigned int s=0; s<_chemistry->n_species(); s++)
      context.point_value(_Y_var.species(s), qp_xyz, Y[s]);

    libMesh::Real M = _chemistry->M(_species_idx); // [kg/mol]
    libMesh::Real M_mix = _chemistry->M_mix(Y); // [kg/mol]
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]);

    P /= 101325.0; // convert to [atm]

    return this->kv(P,T,X,M);
    STOP_LOG("operator()","AbsorptionCoeff");
  }

  template<typename Chemistry>
  void AbsorptionCoeff<Chemistry>::derivatives( libMesh::FEMContext & context,
                                                const libMesh::Point & qp_xyz,
                                                const libMesh::Real & JxW,
                                                const unsigned int qoi_index,
                                                const libMesh::Real /*time*/)
  {
    START_LOG("derivatives()","AbsorptionCoeff");
    std::vector<libMesh::Point> qp(1);
    qp[0] = qp_xyz;

    libMesh::DenseSubVector<libMesh::Number> & dT = context.get_qoi_derivatives(qoi_index, _T_var.T());
    libMesh::DenseSubVector<libMesh::Number> & dP = context.get_qoi_derivatives(qoi_index, _P_var.p());
    libMesh::DenseSubVector<libMesh::Number> & dY = context.get_qoi_derivatives(qoi_index, _Y_var.species(_species_idx));


    libMesh::UniquePtr< libMesh::FEBase > T_fe = libMesh::FEBase::build(context.get_elem().dim(), context.get_element_fe(_T_var.T())->get_fe_type() );
    T_fe->get_phi();
    T_fe->reinit(&(context.get_elem()),&qp);
    const std::vector<std::vector<libMesh::Real> > & T_phi = T_fe->get_phi();

    libMesh::UniquePtr< libMesh::FEBase > P_fe = libMesh::FEBase::build(context.get_elem().dim(), context.get_element_fe(_P_var.p())->get_fe_type() );
    P_fe->get_phi();
    P_fe->reinit(&(context.get_elem()),&qp);
    const std::vector<std::vector<libMesh::Real> > & P_phi = P_fe->get_phi();
 
    libMesh::UniquePtr< libMesh::FEBase > Y_fe = libMesh::FEBase::build(context.get_elem().dim(), context.get_element_fe(_Y_var.species(_species_idx))->get_fe_type() );
    Y_fe->get_phi();
    Y_fe->reinit(&(context.get_elem()),&qp);
    const std::vector<std::vector<libMesh::Real> > & Y_phi = Y_fe->get_phi();

    libMesh::Real T,p,thermo_p; // temperature, hydrostatic pressure, thermodynamic pressure
    std::vector<libMesh::Real> Y(_chemistry->n_species()); // mass fractions

    context.point_value(_T_var.T(), qp_xyz, T); // [K]

    if (_calc_thermo_pressure) {
      libmesh_not_implemented();
    } else {
      thermo_p = _thermo_pressure;
    }

    context.point_value(_P_var.p(), qp_xyz, p); // [Pa]

    libMesh::Real P = p + thermo_p; // total pressure [Pa]
    libmesh_assert_greater(P,0.0);

    // all mass fractions needed to get M_mix
    for (unsigned int s=0; s<_chemistry->n_species(); s++)
      context.point_value(_Y_var.species(s), qp_xyz, Y[s]);

    libMesh::Real MW = _chemistry->M(_species_idx); // [kg/mol]
    libMesh::Real MW_mix = _chemistry->M_mix(Y); // [kg/mol]
    libMesh::Real X = _chemistry->X(_species_idx,MW_mix,Y[_species_idx]);

    P /= 101325.0; // convert to [atm]

    for (unsigned int i=_min_index; i<=_max_index; i++)
      {
        // linecenter wavenumber
        libMesh::Real nu0 = _hitran->nu0(i);

        // air pressure-induced line shift
        libMesh::Real d_air = _hitran->delta_air(i);

        // pressure shift of the linecenter wavenumber
        libMesh::Real nu = nu0 + d_air*(P/_Pref);

        // linestrength
        libMesh::Real S = Sw(T,nu,i);

        // collisional FWHM [cm^-1]
        libMesh::Real nu_c = this->nu_C(T,X,P,i);

        // Doppler FWHM [cm^-1]
        libMesh::Real nu_D = this->nu_D(nu,T,MW);

        // Voigt profile [cm^-1]
        libMesh::Real phi_V = this->voigt(nu_D,nu_c,nu);

        // no velocity dependence

        // temperature deriv
        for (unsigned int j=0; j<dT.size(); j++)
          dT(j) += d_kv_dT(nu_D,nu_c,nu,P,T,T_phi[j][0],MW,X,S,phi_V,i)*JxW;

        // pressure deriv
        for (unsigned int j=0; j<dP.size(); j++)
          dP(j) += d_kv_dP(nu_D,nu_c,nu,P,P_phi[j][0],T,MW,X,S,phi_V,i)*JxW;

        // mass fraction deriv
        for (unsigned int j=0; j<dY.size(); j++)
          dY(j) += d_kv_dY(nu_D,nu_c,nu,P,T,MW,MW_mix,X,Y_phi[j][0],S,phi_V,i)*JxW;
      }
    STOP_LOG("derivatives()","AbsorptionCoeff");
  }

  template<typename Chemistry>
  void AbsorptionCoeff<Chemistry>::operator()( const libMesh::FEMContext & /*context*/,
                                               const libMesh::Point & /*p*/,
                                               const libMesh::Real /*time*/,
                                               libMesh::DenseVector<libMesh::Real> & /*output*/)
  {
    libmesh_not_implemented();
  }

  template<typename Chemistry>
  libMesh::UniquePtr<libMesh::FEMFunctionBase<libMesh::Real> > AbsorptionCoeff<Chemistry>::clone() const
  {
    libmesh_not_implemented();
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::kv(libMesh::Real P,libMesh::Real T,libMesh::Real X,libMesh::Real M)
  {
    libMesh::Real kv = 0.0;

    for (unsigned int i=_min_index; i<=_max_index; i++)
      {
        // linecenter wavenumber
        libMesh::Real nu0 = _hitran->nu0(i);

        // air pressure-induced line shift
        libMesh::Real d_air = _hitran->delta_air(i);

        // pressure shift of the linecenter wavenumber
        libMesh::Real nu = nu0 + d_air*(P/_Pref);

        // linestrength
        libMesh::Real S = Sw(T,nu,i);

        // collisional FWHM [cm^-1]
        libMesh::Real nu_c = this->nu_C(T,X,P,i);

        // Doppler FWHM [cm^-1]
        libMesh::Real nu_D = this->nu_D(nu,T,M);

        // Voigt profile [cm^-1]
        libMesh::Real phi_V = this->voigt(nu_D,nu_c,nu);

        // absorption coefficient [cm^-1]
        kv += S*P*X*phi_V;
      }

    return kv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::Sw(libMesh::Real T, libMesh::Real nu, unsigned int i)
  {
    // isotopologue
    unsigned int iso = _hitran->isotopologue(i);

    // linestrength
    libMesh::Real sw0 = _hitran->sw(i);

    // lower state energy of transition
    libMesh::Real E = _hitran->elower(i);

    // partition function at reference temp
    libMesh::Real QT0 = _hitran->partition_function(_T0,iso);

    // partition function at current temp
    libMesh::Real QT = _hitran->partition_function(T,iso);

    libMesh::Real S = sw0 * (QT0/QT) * std::exp(-E*_rad_coeff*( (1.0/T) - (1.0/_T0) )) * ( 1.0-std::exp(-_rad_coeff*nu/T) ) * pow(1.0-std::exp(-_rad_coeff*nu/_T0),-1.0);

    // convert linestrength units to [cm^-2 atm^-1]
    libMesh::Real loschmidt = (101325.0*1.0e-6)/(T*Constants::Boltzmann);
    S = S*loschmidt;

    return S;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::nu_D(libMesh::Real nu,libMesh::Real T,libMesh::Real M)
  {
    libMesh::Real k = Constants::Boltzmann;
    libMesh::Real c = Constants::c_vacuum;
    libMesh::Real NA = Constants::Avogadro;

    return (nu/c)*std::sqrt( ( 8.0*k*T*std::log(2.0) )/( M/NA ) );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::nu_C(libMesh::Real T,libMesh::Real X,libMesh::Real P, unsigned int index)
  {
    libMesh::Real g_self = _hitran->gamma_self(index);
    libMesh::Real g_air = _hitran->gamma_air(index);
    libMesh::Real n = _hitran->n_air(index);

    return 2.0*P*pow(_T0/T,n) * ( X*g_self + (1.0-X)*g_air );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::voigt(libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu)
  {
    //repeated term
    libMesh::Real root_ln2 = std::sqrt(std::log(2.0));

    libMesh::Real a = root_ln2*nu_c/nu_D;
    libMesh::Real w = 2*root_ln2*(_nu-nu)/nu_D;

    // Voigt coefficient
    libMesh::Real V = 0.0;

    for(int i=0; i<4; i++)
      {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        V += ( Ci*(a-Ai) + Di*(w-Bi) )/( (a-Ai)*(a-Ai) + (w-Bi)*(w-Bi) );
      }

    libMesh::Real phi_V = (2.0*root_ln2)/(std::sqrt(Constants::pi)*nu_D)*V;

    return phi_V;
  }

  template<typename Chemistry>
  void AbsorptionCoeff<Chemistry>::init_voigt() {
    _voigt_coeffs.resize(4);
    for (int i=0; i<4; i++)
      _voigt_coeffs[i].resize(4);

    _voigt_coeffs[0][0] = -1.2150;
    _voigt_coeffs[0][1] = -1.3509;
    _voigt_coeffs[0][2] = -1.2150;
    _voigt_coeffs[0][3] = -1.3509;
    _voigt_coeffs[1][0] =  1.2359;
    _voigt_coeffs[1][1] =  0.3786;
    _voigt_coeffs[1][2] = -1.2359;
    _voigt_coeffs[1][3] = -0.3786;
    _voigt_coeffs[2][0] = -0.3085;
    _voigt_coeffs[2][1] =  0.5906;
    _voigt_coeffs[2][2] = -0.3085;
    _voigt_coeffs[2][3] =  0.5906;
    _voigt_coeffs[3][0] =  0.0210;
    _voigt_coeffs[3][1] = -1.1858;
    _voigt_coeffs[3][2] = -0.0210;
    _voigt_coeffs[3][3] =  1.1858;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_kv_dT(libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu,
                                                    libMesh::Real P, libMesh::Real T, libMesh::Real T_phi, libMesh::Real MW,
                                                    libMesh::Real X, libMesh::Real S, libMesh::Real V, unsigned int i)
  {
    libMesh::Real dV_dT = this->d_voigt_dT(nu_D,nu_c,nu,T,T_phi,X,P,MW,i);
    libMesh::Real dS_dT = this->dS_dT(T,T_phi,nu,i);

    libMesh::Real deriv = X*P*( S*dV_dT + V*dS_dT );
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_kv_dP(libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu,
                                                    libMesh::Real P, libMesh::Real P_phi, libMesh::Real T, libMesh::Real MW,
                                                    libMesh::Real X, libMesh::Real S, libMesh::Real V, unsigned int i)
  {
    libMesh::Real dV_dP = this->d_voigt_dP(nu_D,nu_c,nu,T,MW,X,P_phi,i);
    libMesh::Real dS_dP = this->dS_dP(T,nu,P_phi,i);

    libMesh::Real deriv = X*( S*(P*dV_dP + V*P_phi) + P*V*dS_dP );
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_kv_dY(libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu,
                                                    libMesh::Real P, libMesh::Real T, libMesh::Real MW, libMesh::Real MW_mix,
                                                    libMesh::Real X, libMesh::Real Y_phi, libMesh::Real S, libMesh::Real V, unsigned int i)
  {
    libMesh::Real Y = X*MW/MW_mix;

    libMesh::Real dV_dY = this->d_voigt_dY(nu_D,nu_c,nu,T,P,MW,MW_mix,Y_phi,i);
    libMesh::Real dX = (MW_mix*MW_mix)/MW * ( (1.0/MW_mix) - (Y/MW) ) * Y_phi;

    libMesh::Real deriv = P*S*(X*dV_dY + V*dX);
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::dS_dT(libMesh::Real T, libMesh::Real T_phi, libMesh::Real nu, unsigned int i)
  {
    unsigned int iso = _hitran->isotopologue(i);
    libMesh::Real sw0 = _hitran->sw(i);
    libMesh::Real E = _hitran->elower(i);

    libMesh::Real QT0 = _hitran->partition_function(_T0,iso);
    libMesh::Real QT = _hitran->partition_function(T,iso);

    libMesh::Real const_terms = sw0 * _Pref/(Constants::Boltzmann) * QT0 * std::pow(1.0-std::exp(-_rad_coeff*nu/_T0),-1.0);

    libMesh::Real A = 1.0/T;
    libMesh::Real dA = (-1.0/(T*T))*T_phi;

    libMesh::Real B = 1.0/QT;
    libMesh::Real dB = -1.0/(QT*QT) * dQ_dT(T,iso) * T_phi;

    libMesh::Real C = std::exp( -_rad_coeff*E*(1.0/T - 1.0/_T0) );
    libMesh::Real dC = std::exp( -_rad_coeff*E*(1.0/T - 1.0/_T0) )*(_rad_coeff*E/(T*T))*T_phi;

    libMesh::Real D = 1.0-std::exp(-_rad_coeff*nu/T);
    libMesh::Real dD = -std::exp(-_rad_coeff*nu/T)*(_rad_coeff*nu/(T*T))*T_phi;

    libMesh::Real deriv = const_terms * ( A*(B*(C*dD+D*dC) + C*D*dB) + B*C*D*dA );
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::dS_dP(libMesh::Real T, libMesh::Real nu, libMesh::Real P_phi, unsigned int i)
  {
    unsigned int iso = _hitran->isotopologue(i);
    libMesh::Real sw0 = _hitran->sw(i);
    libMesh::Real E = _hitran->elower(i);
    libMesh::Real d_air = _hitran->delta_air(i);

    libMesh::Real QT0 = _hitran->partition_function(_T0,iso);
    libMesh::Real QT = _hitran->partition_function(T,iso);

    libMesh::Real A = 1.0-std::exp(-_rad_coeff*nu/T);
    libMesh::Real dA = -std::exp(-_rad_coeff*nu/T)*(-_rad_coeff*d_air/(_Pref*T))*P_phi;

    libMesh::Real B = std::pow(1.0-std::exp(-_rad_coeff*nu/_T0),-1.0);
    libMesh::Real dB = -1.0*std::pow(1.0-std::exp(-_rad_coeff*nu/_T0),-2.0)*-std::exp(-_rad_coeff*nu/_T0)*(-_rad_coeff*d_air/(_Pref*_T0))*P_phi;

    libMesh::Real const_terms = sw0 * _Pref/(T*Constants::Boltzmann) * QT0/QT * std::exp(-E*_rad_coeff*( (1.0/T) - (1.0/_T0) ));

    libMesh::Real deriv = const_terms * ( A*dB + B*dA );
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuD_dT(libMesh::Real nu, libMesh::Real T, libMesh::Real T_phi, libMesh::Real M)
  {
    libMesh::Real k = Constants::Boltzmann;
    libMesh::Real c = Constants::c_vacuum;
    libMesh::Real NA = Constants::Avogadro;

    libMesh::Real deriv = (nu/c) * std::sqrt(8.0*std::log(2.0)*k/(M/NA)) * 0.5 * 1.0/std::sqrt(T)*T_phi;
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuD_dP(libMesh::Real T, libMesh::Real M, libMesh::Real P_phi, unsigned int i)
  {
    libMesh::Real k = Constants::Boltzmann;
    libMesh::Real c = Constants::c_vacuum;
    libMesh::Real NA = Constants::Avogadro;

    libMesh::Real d_air = _hitran->delta_air(i);

    libMesh::Real deriv = (d_air/(c*_Pref))*std::sqrt( ( 8.0*k*T*std::log(2.0) )/( M/NA ) )*P_phi;
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuC_dT( libMesh::Real T, libMesh::Real T_phi,
                                                      libMesh::Real X, libMesh::Real P, unsigned int i)
  {
    libMesh::Real n = _hitran->n_air(i);
    libMesh::Real g_air = _hitran->gamma_air(i);
    libMesh::Real g_self = _hitran->gamma_self(i);

    libMesh::Real deriv = 2.0*P*n*std::pow(T/_T0,n-1)*(1.0/_T0)*T_phi*(X*g_self + (1.0-X)*g_air);
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuC_dP( libMesh::Real T, libMesh::Real X,
                                                      libMesh::Real P_phi, unsigned int i)
  {
    libMesh::Real n = _hitran->n_air(i);
    libMesh::Real g_air = _hitran->gamma_air(i);
    libMesh::Real g_self = _hitran->gamma_self(i);

    libMesh::Real deriv = 2.0*std::pow(T/_T0,n)*(X*g_self + (1.0-X)*g_air)*P_phi;
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuC_dY( libMesh::Real P, libMesh::Real T, libMesh::Real MW,
                                                      libMesh::Real MW_mix, libMesh::Real Y_phi, unsigned int i)
  {
    libMesh::Real n = _hitran->n_air(i);
    libMesh::Real g_air = _hitran->gamma_air(i);
    libMesh::Real g_self = _hitran->gamma_self(i);

    libMesh::Real dX = MW_mix/MW;

    libMesh::Real deriv = 2.0*P*std::pow(T/_T0,n)*(g_self - g_air)*dX*Y_phi;
    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_dT( libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu,
                                                        libMesh::Real T, libMesh::Real T_phi, libMesh::Real X,
                                                        libMesh::Real P, libMesh::Real M, unsigned int i)
  {
    //repeated term
    libMesh::Real root_ln2 = std::sqrt(std::log(2.0));

    libMesh::Real a = root_ln2*nu_c/nu_D;
    libMesh::Real w = 2*root_ln2*(_nu-nu)/nu_D;

    // Broadening derivatives
    libMesh::Real dnuC_dT = d_nuC_dT(T,T_phi,X,P,i);
    libMesh::Real dnuD_dT = d_nuD_dT(nu,T,T_phi,M);

    // Voigt coefficients
    libMesh::Real V0 = 0.0;
    libMesh::Real V1 = 0.0;

    // Voigt parameter derivatives
    libMesh::Real da_dT = root_ln2*(nu_D*dnuC_dT - nu_c*dnuD_dT)/(nu_D*nu_D);
    libMesh::Real dw_dT = 2.0*root_ln2*(_nu-nu)/(nu_D*nu_D) * (-dnuD_dT);

    for(int i=0; i<4; i++)
     {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        libMesh::Real aAi = (a-Ai);
        libMesh::Real wBi = (w-Bi);

        libMesh::Real aAi2 = aAi*aAi;
        libMesh::Real wBi2 = wBi*wBi;

        V0 += ( (aAi2+wBi2)*(Ci*da_dT + Di*dw_dT) - (Ci*aAi + Di*wBi)*( 2.0*aAi*da_dT + 2.0*wBi*dw_dT ) )/std::pow(aAi2+wBi2, 2.0);
        V1 += ( Ci*aAi + Di*wBi )/( aAi2+wBi2 );
     }

    libMesh::Real deriv = ( (2.0*root_ln2)/(std::sqrt(Constants::pi)*nu_D) )*V0 + V1*( (2.0*root_ln2*(-dnuD_dT))/(std::sqrt(Constants::pi)*nu_D*nu_D) );

    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_dP( libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu,
                                                        libMesh::Real T, libMesh::Real M, libMesh::Real X,
                                                        libMesh::Real P_phi, unsigned int i)
  {
    //repeated term
    libMesh::Real root_ln2 = std::sqrt(std::log(2.0));

    libMesh::Real d_air = _hitran->delta_air(i);

    libMesh::Real a = root_ln2*nu_c/nu_D;
    libMesh::Real w = 2*root_ln2*(_nu-nu)/nu_D;

    // Broadening derivatives
    libMesh::Real dnuC_dP = d_nuC_dP(T,X,P_phi,i);
    libMesh::Real dnuD_dP = d_nuD_dP(T,M,P_phi,i);

    // Voigt coefficients
    libMesh::Real V0 = 0.0;
    libMesh::Real V1 = 0.0;

    // Voigt parameter derivatives
    libMesh::Real da_dP = root_ln2*(nu_D*dnuC_dP - nu_c*dnuD_dP)/(nu_D*nu_D);
    libMesh::Real dw_dP = 2.0*root_ln2*( nu_D*(-d_air/_Pref) - (_nu-nu)*dnuD_dP )/(nu_D*nu_D);

    for(int i=0; i<4; i++)
     {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        libMesh::Real aAi = (a-Ai);
        libMesh::Real wBi = (w-Bi);

        libMesh::Real aAi2 = aAi*aAi;
        libMesh::Real wBi2 = wBi*wBi;

        V0 += ( (aAi2+wBi2)*(Ci*da_dP + Di*dw_dP) - (Ci*aAi + Di*wBi)*( 2.0*aAi*da_dP + 2.0*wBi*dw_dP ) )/std::pow(aAi2+wBi2, 2.0);
        V1 += ( Ci*aAi + Di*wBi )/( aAi2+wBi2 );
     }

    libMesh::Real deriv = ( (2.0*root_ln2)/(std::sqrt(Constants::pi)*nu_D) )*V0 + V1*( (2.0*root_ln2*(-dnuD_dP))/(std::sqrt(Constants::pi)*nu_D*nu_D) );

    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_dY( libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu,
                                                        libMesh::Real T, libMesh::Real P, libMesh::Real MW,
                                                        libMesh::Real MW_mix, libMesh::Real Y_phi, unsigned int i)
  {
    //repeated term
    libMesh::Real root_ln2 = std::sqrt(std::log(2.0));

    libMesh::Real a = root_ln2*nu_c/nu_D;
    libMesh::Real w = 2*root_ln2*(_nu-nu)/nu_D;

    // Voigt coefficient
    libMesh::Real V0 = 0.0;

    // Voigt parameter derivative
    libMesh::Real da_dY = (root_ln2/nu_D)*d_nuC_dY(P,T,MW,MW_mix,Y_phi,i);

    for(int i=0; i<4; i++)
     {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        libMesh::Real aAi = (a-Ai);
        libMesh::Real wBi = (w-Bi);

        libMesh::Real aAi2 = aAi*aAi;
        libMesh::Real wBi2 = wBi*wBi;

        V0 += ( (aAi2+wBi2)*(Ci*da_dY) - (Ci*aAi + Di*wBi)*( 2.0*aAi*da_dY ) )/std::pow(aAi2+wBi2, 2.0);
     }

    libMesh::Real deriv = ( (2.0*root_ln2)/(std::sqrt(Constants::pi)*nu_D) )*V0;

    return deriv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::dQ_dT(libMesh::Real T, unsigned int iso)
  {
    libMesh::Real deriv = _hitran->partition_function_derivative(T,iso);
    return deriv;
  }

#if GRINS_HAVE_ANTIOCH
  template class AbsorptionCoeff<AntiochChemistry>;
#endif

#if GRINS_HAVE_CANTERA
  template class AbsorptionCoeff<CanteraMixture>;
#endif
} //namespace GRINS
