include(ExternalProject)
include_directories(commons)
add_subdirectory(commons)
add_subdirectory(sra)
add_subdirectory(version)
add_subdirectory(workflow)
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_SOURCE_DIR}/lib/block-aligner)

add_executable(srasearch
        ${commons_source_files}
        ${sra_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        srasearch.cpp
        commons/BitManipulateMacros.h)

ExternalProject_Add(
        blockaligner
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND cargo build COMMAND cargo build --release
        BINARY_DIR "${CMAKE_SOURCE_DIR}/lib/block-aligner"
        INSTALL_COMMAND ""
        LOG_BUILD ON)

mmseqs_setup_derived_target(srasearch)
add_dependencies(srasearch local-generated blockaligner)

find_library(BLOCK_ALIGNER block_aligner HINTS "${CMAKE_SOURCE_DIR}/lib/block-aligner/target/release")
target_link_libraries(srasearch version atomic pthread m dl ${BLOCK_ALIGNER})

install(TARGETS srasearch DESTINATION bin)