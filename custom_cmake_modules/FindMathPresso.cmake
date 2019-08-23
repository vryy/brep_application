#
# Looks for MathPresso packages
#

message("MATHPRESSO_DIR:" ${MATHPRESSO_DIR})
IF(MATHPRESSO_DIR)
    SET(MATHPRESSO_INC_DIRS ${MATHPRESSO_DIR}/include)
    SET(MATHPRESSO_LIB_DIRS ${MATHPRESSO_DIR}/lib)

    INCLUDE_DIRECTORIES(${MATHPRESSO_INC_DIRS})

    FIND_LIBRARY(AUX_MATHPRESSO mathpresso ${MATHPRESSO_LIB_DIRS} NO_DEFAULT_PATH)
    
    IF(AUX_MATHPRESSO)
        SET(MATHPRESSO_LIBRARIES ${AUX_MATHPRESSO})
        SET(MATHPRESSO_FOUND TRUE)
    ELSE()
        SET(MATHPRESSO_FOUND FALSE)
    ENDIF()
ELSE()
    SET(MATHPRESSO_FOUND FALSE)
ENDIF()

