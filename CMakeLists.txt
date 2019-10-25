# libaatm

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# foreach(policy
#     CMP0048
#     CMP0074
#     CMP0077
#   )
#   if(POLICY ${policy})
#     cmake_policy(SET ${policy} NEW)
#   endif()
# endforeach()

project(libaatm VERSION 1.0.0 LANGUAGES C CXX)

# Auxiliary files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Use GNUInstallDirs to install libraries into correct locations on all
# platforms.
include(GNUInstallDirs)

# Build defaults
include(BuildType)

# External packages


# Internal products

enable_testing()
add_subdirectory(src)
