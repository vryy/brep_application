set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosBRepApplication ****")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/custom_cmake_modules")

######### add include directory here
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
#######################################################################

FIND_PACKAGE(MathPresso REQUIRED)
IF(${MATHPRESSO_FOUND} MATCHES TRUE)
    MESSAGE("mathpresso is found at " ${MATHPRESSO_DIR})
    ADD_DEFINITIONS(-DBREP_APPLICATION_USE_MASHPRESSO)
ELSE()
    MESSAGE("mathpresso is not found")
ENDIF()

FIND_PACKAGE(OpenCasCade REQUIRED)
IF(${OCC_FOUND} MATCHES TRUE)
    MESSAGE("OpenCasCade libraries:" ${OCC_LIBRARIES})
    include_directories(${OCC_INCLUDE_DIR})
    ADD_DEFINITIONS(-DBREP_APPLICATION_USE_OPENCASCADE)
ELSE()
    MESSAGE("OpenCasCade is not found")
ENDIF()

######### add compilation units here
set( KRATOS_BREP_APPLICATION_SOURCES
    # main application compilation unit
    ${CMAKE_CURRENT_SOURCE_DIR}/brep_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_algebra/brep.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_mesh_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_algebra_1_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_algebra_2_to_python.cpp
)

IF(${OCC_FOUND} MATCHES TRUE)
set( KRATOS_BREP_APPLICATION_SOURCES ${KRATOS_BREP_APPLICATION_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/MakeSphere.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/MakeBottle.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/ReadSTEP.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/WriteSTEP.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_occ_to_python.cpp
)
endif()

set( KRATOS_BREP_APPLICATION_SOURCES ${KRATOS_BREP_APPLICATION_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/brep_python_application.cpp
)
#######################################################################

## define library Kratos which defines the basic python interface
add_library(KratosBRepApplication SHARED ${KRATOS_BREP_APPLICATION_SOURCES})
target_link_libraries(KratosBRepApplication KratosCore)
IF(${MATHPRESSO_FOUND} MATCHES TRUE)
    target_link_libraries(KratosBRepApplication ${MATHPRESSO_LIBRARIES})
ENDIF()
IF(${OCC_FOUND} MATCHES TRUE)
    target_link_libraries(KratosBRepApplication ${OCC_LIBRARIES})
ENDIF()
set_target_properties(KratosBRepApplication PROPERTIES PREFIX "")
install(TARGETS KratosBRepApplication DESTINATION libs )
#######################################################################

## changing the .dll suffix to .pyd for Windows compatibility
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosBRepApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#######################################################################

## install application's python script if required
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
#######################################################################

## Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/BRepApplication.py" DESTINATION KratosMultiphysics )
install(FILES ${MATHPRESSO_LIBRARIES} DESTINATION libs )
#######################################################################
