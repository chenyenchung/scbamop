cmake_minimum_required(VERSION 3.18)
project(scbamop VERSION 0.5.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Sanitizer options
option(ENABLE_ASAN "Enable AddressSanitizer" ON)
option(ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" ON)
option(ENABLE_LSAN "Enable LeakSanitizer" OFF)
option(ENABLE_TSAN "Enable ThreadSanitizer" OFF)

# Find HTSlib
find_package(PkgConfig QUIET)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(HTSLIB QUIET htslib)
endif()

set(HTSLIB_SEARCH_DIRS
    $ENV{HTSLIB_ROOT}
    ${HTSLIB_ROOT}
    /usr
    /usr/local
    /opt/homebrew
)

# Include dir
find_path(HTSlib_INCLUDE_DIR
    NAMES htslib/sam.h
    PATHS ${HTSLIB_SEARCH_DIRS} ${HTSLIB_PKGCONF_INCLUDE_DIRS}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

# Library
find_library(HTSlib_LIBRARY
    NAMES libhts.so.3 libhts.so.2 libhts.so libhts.3.dylib libhts.2.dylib libhts.dylib
    PATHS ${HTSLIB_SEARCH_DIRS} ${HTSLIB_PKGCONF_LIBRARY_DIRS}
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH
)

if(NOT HTSlib_INCLUDE_DIR OR NOT HTSlib_LIBRARY)
    message(FATAL_ERROR "HTSlib not found. Please install HTSlib development package.")
endif()

message(STATUS "HTSlib include: ${HTSlib_INCLUDE_DIR}")
message(STATUS "HTSlib library: ${HTSlib_LIBRARY}")

# Create executable
add_executable(${PROJECT_NAME}
    src/main.c
    src/hash.c
    src/utils.c
    src/sort.c
    src/dedup_3pass.c
)

target_include_directories(${PROJECT_NAME} PUBLIC ${HTSlib_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${HTSlib_LIBRARY})

# Find threads
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME} Threads::Threads)

# Configure sanitizers
set(SANITIZER_FLAGS "")
set(SANITIZER_LINK_FLAGS "")

if(ENABLE_ASAN)
    list(APPEND SANITIZER_FLAGS "-fsanitize=address")
    list(APPEND SANITIZER_LINK_FLAGS "-fsanitize=address")
    message(STATUS "AddressSanitizer enabled")
endif()

if(ENABLE_UBSAN)
    list(APPEND SANITIZER_FLAGS "-fsanitize=undefined")
    list(APPEND SANITIZER_LINK_FLAGS "-fsanitize=undefined")
    message(STATUS "UndefinedBehaviorSanitizer enabled")
endif()

if(ENABLE_LSAN)
    list(APPEND SANITIZER_FLAGS "-fsanitize=leak")
    list(APPEND SANITIZER_LINK_FLAGS "-fsanitize=leak")
    message(STATUS "LeakSanitizer enabled")
endif()

if(ENABLE_TSAN)
    list(APPEND SANITIZER_FLAGS "-fsanitize=thread")
    list(APPEND SANITIZER_LINK_FLAGS "-fsanitize=thread")
    message(STATUS "ThreadSanitizer enabled")
endif()

# Apply sanitizer flags
if(SANITIZER_FLAGS)
    list(APPEND SANITIZER_FLAGS "-fno-omit-frame-pointer" "-g")
    target_compile_options(${PROJECT_NAME} PRIVATE ${SANITIZER_FLAGS})
    target_link_options(${PROJECT_NAME} PRIVATE ${SANITIZER_LINK_FLAGS})
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
