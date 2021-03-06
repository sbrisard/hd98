add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)

find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(${PROJECT_NAME} INTERFACE ${MATH_LIBRARY})
endif()

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}
  DESTINATION include)

install(
  TARGETS ${PROJECT_NAME}
  DESTINATION lib
  EXPORT ${PROJECT_NAME}_targets)

install(
  EXPORT ${PROJECT_NAME}_targets
  FILE ${PROJECT_NAME}-targets.cmake
  DESTINATION lib/cmake/${PROJECT_NAME})

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
        INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
        NO_SET_AND_CHECK_MACRO
        NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# write_basic_package_version_file(
#   "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
#   VERSION "${Tutorial_VERSION_MAJOR}.${Tutorial_VERSION_MINOR}"
#   COMPATIBILITY AnyNewerVersion
# )

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
        DESTINATION lib/cmake/${PROJECT_NAME})
