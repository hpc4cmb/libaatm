
# Name of the static library
set(AATM_STATIC aatm_static)

# Name of loadable module / shared library
set(AATM_MOD aatm)

# Library sources
set(AATM_SOURCES
    src/ATMAngle.cpp
    src/ATMError.cpp
    src/ATMException.cpp
    src/ATMFrequency.cpp
    src/ATMInverseLength.cpp
    src/ATMLength.cpp
    src/ATMMassDensity.cpp
    src/ATMNumberDensity.cpp
    src/ATMOpacity.cpp
    src/ATMPercent.cpp
    src/ATMPressure.cpp
    src/ATMProfile.cpp
    src/ATMRefractiveIndex.cpp
    src/ATMRefractiveIndexProfile.cpp
    src/ATMSkyStatus.cpp
    src/ATMSpectralGrid.cpp
    src/ATMTemperature.cpp
    src/ATMVersion.cpp
    src/ATMWaterVaporRadiometer.cpp
    src/ATMWVRMeasurement.cpp
)

# Add the internal object library target

add_library(aatmobj OBJECT ${AATM_SOURCES})
set_property(TARGET aatmobj PROPERTY POSITION_INDEPENDENT_CODE 1)

target_include_directories(aatmobj PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    src
)

# Static library

add_library(${AATM_STATIC} STATIC $<TARGET_OBJECTS:aatmobj>)

target_include_directories(${AATM_STATIC} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

# Now define rules for the shared library

add_library(${AATM_MOD} SHARED $<TARGET_OBJECTS:aatmobj>)

target_include_directories(${AATM_MOD} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
    TARGETS
    ${AATM_MOD}
    ${AATM_STATIC}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
