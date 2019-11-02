
# Name of the internal static library
set(AATM_LIB aatm_)

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

# Add the internal library target

add_library(${AATM_LIB} STATIC ${AATM_SOURCES})

target_include_directories(${AATM_LIB} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

# Now define rules for the shared library

add_library(${AATM_MOD} SHARED ATM.cpp)

target_link_libraries(${AATM_MOD} PRIVATE ${AATM_LIB})

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ${AATM_MOD} DESTINATION ${CMAKE_INSTALL_LIBDIR})
