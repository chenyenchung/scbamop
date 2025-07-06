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

# Build HTSlib from submodule using ExternalProject
include(ExternalProject)

# Find Git
find_package(Git QUIET)
if(NOT GIT_FOUND)
    message(FATAL_ERROR "Git not found. Git is required to fetch HTSlib submodule.")
endif()

# Check for required build tools
find_program(MAKE_EXECUTABLE make)
if(NOT MAKE_EXECUTABLE)
    message(FATAL_ERROR "make not found. On Debian/Ubuntu: apt-get install build-essential")
endif()

find_program(GCC_EXECUTABLE gcc)
if(NOT GCC_EXECUTABLE)
    message(FATAL_ERROR "gcc not found. On Debian/Ubuntu: apt-get install build-essential")
endif()

find_program(AR_EXECUTABLE ar)
if(NOT AR_EXECUTABLE)
    message(FATAL_ERROR "ar not found. On Debian/Ubuntu: apt-get install build-essential")
endif()

find_program(RANLIB_EXECUTABLE ranlib)
if(NOT RANLIB_EXECUTABLE)
    message(FATAL_ERROR "ranlib not found. On Debian/Ubuntu: apt-get install build-essential")
endif()

message(STATUS "Build tools found: make=${MAKE_EXECUTABLE}, gcc=${GCC_EXECUTABLE}")

# Check for HTSlib runtime dependencies
find_package(ZLIB)
if(NOT ZLIB_FOUND)
    message(FATAL_ERROR "zlib not found. On Debian/Ubuntu: apt-get install zlib1g-dev")
endif()

find_package(BZip2)
if(NOT BZIP2_FOUND)
    message(FATAL_ERROR "bzip2 not found. On Debian/Ubuntu: apt-get install libbz2-dev")
endif()

find_package(LibLZMA)
if(NOT LIBLZMA_FOUND)
    message(FATAL_ERROR "liblzma not found. On Debian/Ubuntu: apt-get install liblzma-dev")
endif()

find_package(CURL)
if(NOT CURL_FOUND)
    message(FATAL_ERROR "libcurl not found. On Debian/Ubuntu: apt-get install libcurl4-openssl-dev")
endif()

message(STATUS "HTSlib dependencies found: zlib=${ZLIB_VERSION_STRING}, bzip2=${BZIP2_VERSION_STRING}, lzma=${LIBLZMA_VERSION_STRING}, curl=${CURL_VERSION_STRING}")

set(HTSLIB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/htslib)
set(HTSLIB_INCLUDE_DIR ${HTSLIB_SOURCE_DIR})
set(HTSLIB_LIBRARY ${HTSLIB_SOURCE_DIR}/libhts.a)

# Initialize and update submodule if needed
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/htslib/.git")
    message(STATUS "Initializing HTSlib submodule...")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive third_party/htslib
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT
    )
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "Failed to initialize HTSlib submodule")
    endif()
endif()

# Ensure we're on the correct version (1.22)
set(HTSLIB_REQUIRED_VERSION "1.22")
execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --exact-match
    WORKING_DIRECTORY ${HTSLIB_SOURCE_DIR}
    OUTPUT_VARIABLE HTSLIB_CURRENT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

if(NOT "${HTSLIB_CURRENT_VERSION}" STREQUAL "${HTSLIB_REQUIRED_VERSION}")
    message(STATUS "Switching HTSlib to version ${HTSLIB_REQUIRED_VERSION} (current: ${HTSLIB_CURRENT_VERSION})")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} checkout ${HTSLIB_REQUIRED_VERSION}
        WORKING_DIRECTORY ${HTSLIB_SOURCE_DIR}
        RESULT_VARIABLE GIT_CHECKOUT_RESULT
    )
    if(NOT GIT_CHECKOUT_RESULT EQUAL "0")
        message(FATAL_ERROR "Failed to checkout HTSlib version ${HTSLIB_REQUIRED_VERSION}")
    endif()
else()
    message(STATUS "HTSlib is at correct version: ${HTSLIB_REQUIRED_VERSION}")
endif()

# Validate HTSlib source directory and key files
if(NOT EXISTS ${HTSLIB_SOURCE_DIR})
    message(FATAL_ERROR "HTSlib source directory not found: ${HTSLIB_SOURCE_DIR}. CMake should have initialized this automatically - please check for errors above.")
endif()

if(NOT EXISTS ${HTSLIB_SOURCE_DIR}/Makefile)
    message(FATAL_ERROR "HTSlib Makefile not found: ${HTSLIB_SOURCE_DIR}/Makefile. The submodule may not be properly initialized.")
endif()

if(NOT EXISTS ${HTSLIB_SOURCE_DIR}/htslib/sam.h)
    message(FATAL_ERROR "HTSlib headers not found: ${HTSLIB_SOURCE_DIR}/htslib/sam.h. The submodule may be incomplete.")
endif()

message(STATUS "HTSlib source validated: ${HTSLIB_SOURCE_DIR}")

ExternalProject_Add(htslib_external
    SOURCE_DIR ${HTSLIB_SOURCE_DIR}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E chdir ${HTSLIB_SOURCE_DIR} ${MAKE_EXECUTABLE} clean
    BUILD_COMMAND ${MAKE_EXECUTABLE} -C ${HTSLIB_SOURCE_DIR} -j4 lib-static
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    LOG_BUILD 1
)

# Create imported target for HTSlib
add_library(hts STATIC IMPORTED)
set_target_properties(hts PROPERTIES
    IMPORTED_LOCATION ${HTSLIB_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES ${HTSLIB_INCLUDE_DIR}
)
add_dependencies(hts htslib_external)

message(STATUS "HTSlib include: ${HTSLIB_INCLUDE_DIR}")
message(STATUS "HTSlib library: ${HTSLIB_LIBRARY}")
message(STATUS "HTSlib target: hts")

# Create executable
add_executable(${PROJECT_NAME}
    src/main.c
    src/hash.c
    src/utils.c
    src/sort.c
    src/dedup_3pass.c
)

add_dependencies(${PROJECT_NAME} hts)
target_link_libraries(${PROJECT_NAME} hts z m bz2 lzma curl pthread)

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
    
    # Force ASan to be linked first to avoid runtime warnings
    if(ENABLE_ASAN)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU")
            # For GCC, use -static-libasan to avoid LD_PRELOAD issues
            target_link_options(${PROJECT_NAME} PRIVATE -static-libasan)
            # Disable PIE to avoid shadow memory conflicts
            target_compile_options(${PROJECT_NAME} PRIVATE -fno-pie)
            target_link_options(${PROJECT_NAME} PRIVATE -no-pie)
        elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
            # For Clang, add ASan runtime explicitly
            target_link_options(${PROJECT_NAME} PRIVATE -shared-libasan)
            # Disable PIE to avoid shadow memory conflicts
            target_compile_options(${PROJECT_NAME} PRIVATE -fno-pie)
            target_link_options(${PROJECT_NAME} PRIVATE -no-pie)
        endif()
    endif()
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
