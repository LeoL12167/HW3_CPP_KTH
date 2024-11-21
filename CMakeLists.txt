# CMakeList.txt: CMake-Projekt für "HW3_CPP_KTH". Schließen Sie die Quelle ein, und definieren Sie
# projektspezifische Logik hier.
#
cmake_minimum_required (VERSION 3.8)

project ("HW3_CPP_KTH" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(UTEST_HEADER_DIR ${CMAKE_BINARY_DIR}/external/utest)
file(DOWNLOAD
        https://raw.githubusercontent.com/sheredom/utest.h/master/utest.h
        ${UTEST_HEADER_DIR}/utest.h
)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})



# Fügen Sie der ausführbaren Datei dieses Projekts eine Quelle hinzu.
add_executable (HW3_CPP_KTH "HW3_CPP_KTH.cpp" "HW3_CPP_KTH.h")
add_executable (with_cache src/main.cpp src/geometry.hpp)
add_executable (no_cache src/main.cpp src/geometry.hpp)

target_compile_definitions(with_cache PRIVATE USE_CACHE)

target_link_libraries(with_cache ${Boost_LIBRARIES})
target_link_libraries(no_cache ${Boost_LIBRARIES})
