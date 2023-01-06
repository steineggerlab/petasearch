include_directories(commons)
add_subdirectory(commons)
add_subdirectory(sra)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(srasearch
        ${commons_source_files}
        ${sra_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        srasearch.cpp
        )

mmseqs_setup_derived_target(srasearch)
add_dependencies(srasearch local-generated)
target_link_libraries(srasearch version block-aligner)

install(TARGETS srasearch DESTINATION bin)
