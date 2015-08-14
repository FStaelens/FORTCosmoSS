```
* main.f90
    |
    - program main
         |
         - subroutine write_output
         |
         - subroutine write_in_log 

* common.f90
    |
    - module gauge_choice
         |
         - function f_alpha
         |
    - module matter
    |
    - module BSSN_var
         |
         - subroutine update_derivatives
         |
    - module cosmov
         |
         - subroutine update_cosmov
         |
    - module metric
         |
         - subroutine update_metric
         
* constraints.f90
    |
    - module constraints
         |
         - subroutine H_cons
         |
         - subroutine M_cons
         |
         - subroutine L2_norm
         |
         - subroutine H_cons_cosm
         |
         
* evolution.f90
    |
    - module evolution
         |
         - subroutine PIRK
         |
         - subroutine make_L1_cosmo
         |
         - subroutine make_L2_cosmo
         |
         - subroutine make_L3_cosmo
         |
         - subroutine make_L1
         |
         - subroutine make_L2
         |
         - subroutine make_L3
         |
         - subroutine make_L2bar
         |
         - subroutine make_L3bar
         
 * fparser.f90	(public parts only)
    |
    - module fparser
         |
         - subroutine initf
         |
         - subroutine parsef
         |
         - subroutine evalf
         |
         - subroutine EvalErrMsg
         |
         - subroutine EvalErrType

 * grid.f90	
    |
    - module grid 
         
 * hydro.f90
    |
    - module hydro
         |
         - subroutine buildlr
         |
         - subroutine hlle
         |
         - subroutine eigenvalue
         |
         - subroutine cons2prim

 * inifile.f90
    |
    - module IniFile
         |
         - subroutine read_input
         |
         - subroutine assign_real
         |
         - subroutine assign_integer
         |
         - subroutine assign_str
         |
         - function isthere
      
 * init.f90
    |
    - module init
         |
         - subroutine init_data
         |
         - subroutine fsub
         |
         - subroutine bcsub
         
 * initial_profiles.f90
    |
    - module initial_profiles
         |
         - function rho_mix
         |
         - function  phi_ix
         |
         - function Psi_phiix
         |
         - function V_phi_ix

* input.f90
    |
    - module input_param
         |
         - subroutine input_list
         |
         - subroutine input_potential_list
         
* math_lib.f90
    |
    - module profiles
         |
         - function logistic
         |
         - function bump
         |
         - function dxbump
         |
         - function sym_gaussian
         |
         - function dxsym_gaussian
         |
    - module root_finding
         |
         - function root
         |
    - module derivatives_fcn
         |
         - function dxf
         |
         - function d2xf
         |
         - function Delta4_x
         |
         - function Delta4_x_i
         |
    - module boundaries
         |
         - subroutine symmetrise
         |
         - subroutine anti_symmetrise
         |
         - subroutine symmetrise_centred
         |
         - subroutine sommerfeld
         |
    - module num_integration
         |
         - function integral
         |
    - module find_value
         |
         - subroutine find_closest
         |
    - module interpolation
         |
         - subroutine POLINT
         
* output.f90
    |
    - module outputs
         |
         - subroutine output_list
         |
         - subroutine append_to_list
         |
         - subroutine print_list
         |
         - subroutine progress
         
* potential.f90
    |
    - module potential
         |
         - subroutine init_potential
         |
         - subroutine V_phi
         |
         - function phi_V
         
* scales.f90
    |
    - module constants
    
* sources.f90
    |
    - module sources
         |
         - subroutine build_hydro_sources
         |
         - subroutine build_matter_sources
```

