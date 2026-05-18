cmake_minimum_required(VERSION 3.25)

project(Burnham VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BURNHAM_BUILD_TESTS "Build Burnham tests" ON)
option(BURNHAM_ENABLE_HTSLIB "Enable htslib-backed binary alignment probing" OFF)

add_library(burnham_core
    src/chain_index.cpp
    src/commands.cpp
    src/cigar.cpp
    src/fs.cpp
    src/hts.cpp
    src/sample_index.cpp
    src/sam.cpp
    src/text.cpp
    src/vcf.cpp
)

target_include_directories(burnham_core PUBLIC include)

if(BURNHAM_ENABLE_HTSLIB)
    find_package(htslib CONFIG REQUIRED)
    set(BURNHAM_HTSLIB_TARGET "")
    foreach(candidate htslib::htslib HTSlib::HTSlib unofficial::htslib::htslib)
        if(TARGET ${candidate})
            set(BURNHAM_HTSLIB_TARGET ${candidate})
            break()
        endif()
    endforeach()
    if(NOT BURNHAM_HTSLIB_TARGET)
        message(FATAL_ERROR "htslib package was found, but no known CMake target was exported")
    endif()
    target_link_libraries(burnham_core PUBLIC ${BURNHAM_HTSLIB_TARGET})
    target_compile_definitions(burnham_core PUBLIC BURNHAM_HAS_HTSLIB=1)
else()
    target_compile_definitions(burnham_core PUBLIC BURNHAM_HAS_HTSLIB=0)
endif()

if(MSVC)
    target_compile_options(burnham_core PRIVATE /W4 /permissive-)
else()
    target_compile_options(burnham_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(burnham src/main.cpp)
target_link_libraries(burnham PRIVATE burnham_core)

install(TARGETS burnham RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES README.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
install(DIRECTORY docs/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/docs FILES_MATCHING PATTERN "*.md")

if(BURNHAM_BUILD_TESTS)
    enable_testing()
    add_executable(burnham_tests tests/unit_tests.cpp)
    target_link_libraries(burnham_tests PRIVATE burnham_core)
    add_test(NAME burnham_unit_tests COMMAND burnham_tests)
endif()
