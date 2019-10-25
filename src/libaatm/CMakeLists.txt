
# Add the library target

add_library(aatm
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

target_include_directories(aatm PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Hardcode this to "lib" for now, since lib32/lib64 is a pain
# to use on many systems.
#install(TARGETS aatm DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS aatm DESTINATION lib)
