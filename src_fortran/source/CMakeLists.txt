cmake_minimum_required(VERSION 3.12)
project(KrakenFortran LANGUAGES Fortran)

# Compiler flags
set(CMAKE_Fortran_FLAGS "-cpp -fPIC -O3 -finline-functions -fno-strength-reduce -fomit-frame-pointer -falign-functions=2")

# Define the shared library and its source files
add_library(kraken SHARED
    ../krakmod.f90
    ../sdrdrmod.f90
    ../zsecx.f90
    ../zbrentx.f90
    ../bcimp.f90
    ../kuping.f90
    ../sinvitd.f90
    ../mergev.f90
    ../weight.f90
    ../twersk.f90
    ../subtab.f90
    ../readin.f90
    ../errout.f90
    ../sorti.f90
    ../splinec.f90
    ../kraken.f90
)

# Standard installation rules (handles Windows 'bin' vs Unix 'lib' automatically)
include(GNUInstallDirs)
install(TARGETS kraken
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
