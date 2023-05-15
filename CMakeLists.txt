cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
# keep old policy for setting implicit link libraries
# zlib causes issues in static builds otherwise
cmake_policy(SET CMP0060 OLD)

project(srasearch C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake")

set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(lib/mmseqs)

add_subdirectory(lib/corrosion)
# don't try to link to gcc_s, its not needed anyway
list(REMOVE_ITEM Rust_CARGO_TARGET_LINK_NATIVE_LIBS "gcc_s")

if(Rust_VERSION VERSION_LESS "1.68.0" AND NOT IGNORE_RUST_VERSION)
    message(FATAL_ERROR "Rust >= 1.68.0 is required")
endif()

set(RUST_FEATURE "")
if(HAVE_AVX2)
    set(RUST_FEATURE simd_avx2)
elseif(HAVE_SSE4_1 OR HAVE_SSE2)
    set(RUST_FEATURE simd_sse2)
elseif(HAVE_ARM8)
    set(RUST_FEATURE simd_neon)
endif()

if (NATIVE_ARCH AND (RUST_FEATURE STREQUAL ""))
    if(ARM)
        set(RUST_FEATURE simd_neon)
    elseif(X86 OR X64)
        if (HAVE_AVX2_EXTENSIONS)
            set(RUST_FEATURE simd_avx2)
        elseif (HAVE_SSE2_EXTENSIONS)
            set(RUST_FEATURE simd_sse2)
        endif()
    endif()
endif()

if(RUST_FEATURE STREQUAL "")
    message(FATAL_ERROR "No SIMD support for this architecture")
endif()
message("-- Rust Feature: ${RUST_FEATURE}")
corrosion_import_crate(MANIFEST_PATH lib/block-aligner/c/Cargo.toml FEATURES "${RUST_FEATURE}")
include_directories(lib/block-aligner/c)

add_subdirectory(data)
add_subdirectory(src)
