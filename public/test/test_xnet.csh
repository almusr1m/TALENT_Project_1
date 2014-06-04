#!/bin/csh -f
# Regression Testing for XNet
# Selection of test problem controlled by argument.
# All serial test problems = 0
# Thermonuclear SN     with alpha (= 1),  150 species (= 2) or 300 species (= 3)
# Tidally-Crushed WD   with alpha (= 4),  150 species (= 5) or 300 species (= 6)
# Nova        with 16 species CNO (= 7),  169 species (= 8)
# X-ray burst with 16 species CNO (= 9),  304 species (=10)
# Core Collapse SN     with alpha (=11), 1072 network(=12) 
# Neutron Star Wind    with alpha (=13), 3837 species(=14)
# All parallel test problems = 30
# 4 different SN       with alpha (=31),  150 species (=32)

set xnet = ../source/xnetp
echo 'Testing ' $xnet

# TN SN tracer, from Ed Brown, with alpha network 
if ($argv[1] == 0 || $argv[1] == 1) then
  cat test_settings_small Test_Problems/setup_tnsn_alpha >! control
  echo 'Test: Thermonuclear SN with alpha network'
  $xnet
  mv -f net_diag1 net_diag_tnsn_alpha
endif

# TN SN tracer, from Ed Brown, with 150 species network 
if ($argv[1] == 0 || $argv[1] == 2) then
  cat test_settings Test_Problems/setup_tnsn_sn150 >! control
  echo 'Test: Thermonuclear SN with 150 species network'
  $xnet
  mv -f net_diag1 net_diag_tnsn_sn150
endif

# TN SN tracer, from Ed Brown, with 300 species network 
if ($argv[1] == 0 || $argv[1] == 3) then
  cat test_settings Test_Problems/setup_tnsn_si2002 >! control
  echo 'Test: Thermonuclear SN with 200 species network'
  $xnet
  mv -f net_diag1 net_diag_tnsn_si2002
endif

# TI SN tracer, from Stephan Rosswog, with alpha network 
if ($argv[1] == 0 || $argv[1] == 4) then
  cat test_settings_small Test_Problems/setup_tisn_alpha >! control
  echo 'Test: Tidally Induced SN with alpha network'
  $xnet
  mv -f net_diag1 net_diag_tisn_alpha
endif

# TI SN tracer, from Stephan Rosswog, with 150 species network 
if ($argv[1] == 0 || $argv[1] == 5) then
  cat test_settings Test_Problems/setup_tisn_sn150 >! control
  echo 'Test: Tidally Induced SN with 150 species network'
  $xnet
  mv -f net_diag1 net_diag_tisn_sn150
endif

# TI SN tracer, from Stephan Rosswog, with 300 species network 
if ($argv[1] == 0 || $argv[1] == 6) then
  cat test_settings Test_Problems/setup_tisn_si2002 >! control
  echo 'Test: Tidally Induced SN with 300 species network'
  $xnet
  mv -f net_diag1 net_diag_tisn_si2002
endif

# Nova zone, from Sumner Starrfield, with 16 species CNO network 
if ($argv[1] == 0 || $argv[1] == 7) then
  cat test_settings_small Test_Problems/setup_nova_cno >! control
  echo 'Test: Nova with minimal CNO network'
  $xnet
  mv -f net_diag1 net_diag_nova_cno
endif

# Nova zone, from Sumner Starrfield, with 189 species network from Chritian Iliadis
if ($argv[1] == 0 || $argv[1] == 8) then
  cat test_settings Test_Problems/setup_nova_Iliadis >! control
  echo 'Test: Nova with 189 species network'
  $xnet
  mv -f net_diag1 net_diag_nova_Iliadis
endif

# XRB tracer, from Jacob Fisker, with minimal CNO network 
if ($argv[1] == 0 || $argv[1] == 9) then
  cat test_settings_small Test_Problems/setup_xrb_cno >! control
  echo 'Test: X-ray burst with CNO network'
  $xnet
  mv -f net_diag1 net_diag_xrb_cno
endif

# XRB tracer, from Jacob Fisker, with 304 species rp-process network 
if ($argv[1] == 0 || $argv[1] == 10) then
  cat test_settings Test_Problems/setup_xrb_fisker >! control
  echo 'Test: X-ray burst with Fiskers network'
  $xnet
  mv -f net_diag1 net_diag_xrb_fisker
endif

# CC SN zone, from Carla Froehlich, with alpha network 
if ($argv[1] == 0 || $argv[1] == 11) then
  cat test_settings_small Test_Problems/setup_ccsn_alpha >! control
  echo 'Test: Core-Collapse SN with alpha network'
  $xnet
  mv -f net_diag1 net_diag_ccsn_alpha
endif

#  CC SN zone, from Carla Froehlich with 150 species network
if ($argv[1] == 0 || $argv[1] == 12) then
  cat test_settings Test_Problems/setup_ccsn_sn150 >! control
  echo 'Test: Core-Collapse SN with 150 species network'
  $xnet
  mv -f net_diag1 net_diag_ccsn_sn150
endif

#  CC SN zone, from Carla Froehlich with 300 species network
if ($argv[1] == 0 || $argv[1] == 13) then
  cat test_settings Test_Problems/setup_ccsn_si2002 >! control
  echo 'Test: Core-Collapse SN with 300 species network'
  $xnet
  mv -f net_diag1 net_diag_ccsn_si2002
endif

#  CC SN zone, from Carla Froehlich with nu p-process network 
if ($argv[1] == 0 || $argv[1] == 14) then
  cat test_settings Test_Problems/setup_ccsn_nup >! control
  echo 'Test: Core-Collapse SN with nu-p process network'
  $xnet
  mv -f net_diag1 net_diag_ccsn_nup
endif

# Neutino driven wind example, from Josh Beun, with alpha network 
if ($argv[1] == 0 || $argv[1] == 15) then
  cat test_settings_small Test_Problems/setup_nuwind_alpha >! control
  echo 'Test: Neutrino-driven wind with alpha network'
  $xnet
  mv -f net_diag1 net_diag_nuwind_alpha
endif

# Neutino driven wind (Meyer & Brown) with large network for JINA REACLIB v2.0 
if ($argv[1] == 0 || $argv[1] == 16) then
  cat test_settings Test_Problems/setup_nuwind_rprocess >! control
  echo 'Test: Neutrino-driven wind with 4510 species network'
  $xnet
  mv -f net_diag1 net_diag_nuwind_rprocess
endif

# Parallel test, runs 4 different zones using SN150 
if ($argv[1] == 30 || $argv[1] == 31) then
  cat test_settings_parallel Test_Problems/setup_parallel_sn150 >! control
  echo 'Test: 4 different zones in parallel'
  mpirun -n 4 ../source/xnetp_mpi
  mv -f mp_diag0 net_diag_parallel_sn150_1
  mv -f mp_diag1 net_diag_parallel_sn150_2
  mv -f mp_diag2 net_diag_parallel_sn150_3
  mv -f mp_diag3 net_diag_parallel_sn150_4
endif

# Parallel test, runs 4 different zones using alpha
if ($argv[1] == 30 || $argv[1] == 32) then
  cat test_settings_parallel Test_Problems/setup_parallel_alpha >! control
  echo 'Test: 4 different zones in parallel'
  mpirun -n 4 ../source/xnetp_mpi
  mv -f mp_diag0 net_diag_parallel_alpha_1
  mv -f mp_diag1 net_diag_parallel_alpha_2
  mv -f mp_diag2 net_diag_parallel_alpha_3
  mv -f mp_diag3 net_diag_parallel_alpha_4
endif

# NSE initial abundance test
if ($argv[1] == 40 || $argv[1] == 41) then
  cat test_settings_nse Test_Problems/setup_nse_nup >! control
  echo 'Test NSE: Core-Collapse SN with nu-p process network'
  ../source/xnet_nse
  mv -f net_diag1 net_diag_nse_ccsn_nup
endif
