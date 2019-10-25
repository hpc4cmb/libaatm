
#======================================================

add_executable(aatm_test_basic
    AtmBasicTest.cpp
)

target_include_directories(aatm_test_basic PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_basic aatm)

install(TARGETS aatm_test_basic DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_basic COMMAND aatm_test_basic)

#======================================================

add_executable(aatm_test_version
    ATMVersionTest.cpp
)

target_include_directories(aatm_test_version PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_version aatm)

install(TARGETS aatm_test_version DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_version COMMAND aatm_test_version)

#======================================================

add_executable(aatm_test_abinitio
    AbInitioTest.cpp
)

target_include_directories(aatm_test_abinitio PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_abinitio aatm)

install(TARGETS aatm_test_abinitio DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_abinitio COMMAND aatm_test_abinitio)

add_custom_command(TARGET aatm_test_abinitio POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:aatm_test_abinitio>/WVR_MAUNA_KEA
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/WVR_MAUNA_KEA/radiometer_data.dat
    $<TARGET_FILE_DIR:aatm_test_abinitio>/WVR_MAUNA_KEA/
)

#======================================================

add_executable(aatm_test_apex
    ApexTest.cpp
)

target_include_directories(aatm_test_apex PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_apex aatm)

install(TARGETS aatm_test_apex DESTINATION ${CMAKE_INSTALL_BINDIR})

# This executable reliably segfaults...
# add_test(NAME test_apex COMMAND aatm_test_apex)

add_custom_command(TARGET aatm_test_apex POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:aatm_test_apex>/APEX
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/APEX/apex_radiometer_data_MJD53965_skydip.dat
    $<TARGET_FILE_DIR:aatm_test_apex>/APEX/
)

#======================================================

add_executable(aatm_test_length
    LengthTest.cpp
)

target_include_directories(aatm_test_length PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_length aatm)

install(TARGETS aatm_test_length DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_length COMMAND aatm_test_length)

#======================================================

add_executable(aatm_test_profile
    AtmProfileTest.cpp
)

target_include_directories(aatm_test_profile PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_profile aatm)

install(TARGETS aatm_test_profile DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_profile COMMAND aatm_test_profile)

#======================================================

add_executable(aatm_test_skystatus
    SkyStatusTest.cpp
)

target_include_directories(aatm_test_skystatus PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_skystatus aatm)

install(TARGETS aatm_test_skystatus DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_skystatus COMMAND aatm_test_skystatus)

#======================================================

add_executable(aatm_test_skystatusinv
    SkyStatusTest_TINV.cpp
)

target_include_directories(aatm_test_skystatusinv PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_skystatusinv aatm)

install(TARGETS aatm_test_skystatusinv DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_skystatusinv COMMAND aatm_test_skystatusinv)

add_custom_command(TARGET aatm_test_skystatusinv POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:aatm_test_skystatusinv>/T-INVERSION
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/T-INVERSION/invt.dat
    $<TARGET_FILE_DIR:aatm_test_skystatusinv>/T-INVERSION/
)

#======================================================

add_executable(aatm_test_refindex
    RefractiveIndexProfileTest.cpp
)

target_include_directories(aatm_test_refindex PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_refindex aatm)

install(TARGETS aatm_test_refindex DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_refindex COMMAND aatm_test_refindex)

#======================================================

add_executable(aatm_test_specgrid
    SpectralGridTest.cpp
)

target_include_directories(aatm_test_specgrid PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_specgrid aatm)

install(TARGETS aatm_test_specgrid DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_specgrid COMMAND aatm_test_specgrid)

#======================================================

add_executable(aatm_test_sma
    SMATest.cpp
)

target_include_directories(aatm_test_sma PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_sma aatm)

install(TARGETS aatm_test_sma DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_sma COMMAND aatm_test_sma)

add_custom_command(TARGET aatm_test_sma POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:aatm_test_sma>/SMA
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/SMA/SMA_17JUL2006_skydip1.dat
    $<TARGET_FILE_DIR:aatm_test_sma>/SMA/
)

#======================================================

add_executable(aatm_test_atm651
    atm651.cpp
)

target_include_directories(aatm_test_atm651 PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libaatm/src"
)

target_link_libraries(aatm_test_atm651 aatm)

install(TARGETS aatm_test_atm651 DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_atm651 COMMAND aatm_test_atm651)
