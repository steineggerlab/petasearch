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
        commons/BitManipulateMacros.h)
mmseqs_setup_derived_target(srasearch)
target_link_libraries(srasearch version)

add_dependencies(srasearch local-generated)

install(TARGETS srasearch DESTINATION bin)
