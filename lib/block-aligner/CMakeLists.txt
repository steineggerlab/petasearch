project(block-aligner)

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/c/libblock_aligner.a
        ${CMAKE_CURRENT_SOURCE_DIR}/c/block_aligner.h
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/c
        COMMAND make block_aligner
        COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/target/release/libblock_aligner.a ${CMAKE_CURRENT_SOURCE_DIR}/c/
)

# https://stackoverflow.com/questions/31274577/custom-target-as-a-target-library-in-cmake
add_custom_target(block_aligner_target DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/c/libblock_aligner.a)

add_library(block_aligner STATIC IMPORTED GLOBAL)
add_dependencies(block_aligner block_aligner_target)
set_target_properties(block_aligner PROPERTIES
        IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/c/libblock_aligner.a
        LINKER_LANGUAGE CXX)