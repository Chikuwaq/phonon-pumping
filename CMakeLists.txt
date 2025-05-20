cmake_minimum_required(VERSION 3.16)

project(PhononPumpingPlotter LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# switch between Debug and Release modes
# 'Release' disables assert().
set(DEBUG_MODE OFF)

if(${DEBUG_MODE})
    # include debug info in the Debug mode
    add_compile_options(-g)
else()
    add_compile_definitions(NDEBUG)
endif()


# Explicitly set architecture to arm64 on Apple Silicon
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "" FORCE)
endif()

# Set compiler flags for Apple Clang on ARM64
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    add_compile_options(-mstrict-align -ftree-vectorize)
endif()

file(GLOB SOURCES "src/*.cpp")
add_executable(PhononPumpingPlotter ${SOURCES})
target_include_directories(PhononPumpingPlotter PUBLIC include)

if(${DEBUG_MODE})
    set_target_properties(PhononPumpingPlotter PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/debug"
    )
else()
    set_target_properties(PhononPumpingPlotter PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/release"
    )
endif()
