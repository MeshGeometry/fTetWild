################################################################################
# Prepare dependencies
################################################################################

# Download external dependencies
include(FloatTetwildDownloadExternal)

################################################################################
# Required libraries
################################################################################

# Sanitizers
if(FLOAT_TETWILD_WITH_SANITIZERS)
	float_tetwild_download_sanitizers()
	find_package(Sanitizers)
endif()

# CL11
if(FLOAT_TETWILD_TOPLEVEL_PROJECT AND NOT TARGET CLI11::CLI11)
	float_tetwild_download_cli11()
	add_subdirectory(${FLOAT_TETWILD_EXTERNAL}/cli11)
endif()

# fmt
if(NOT TARGET fmt::fmt)
	float_tetwild_download_fmt()
	add_subdirectory(${FLOAT_TETWILD_EXTERNAL}/fmt)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
	float_tetwild_download_spdlog()

	# Create interface target
	add_library(spdlog INTERFACE)
	add_library(spdlog::spdlog ALIAS spdlog)
	target_include_directories(spdlog INTERFACE ${FLOAT_TETWILD_EXTERNAL}/spdlog/include)
	target_link_libraries(spdlog INTERFACE fmt::fmt)
	target_compile_definitions(spdlog INTERFACE SPDLOG_FMT_EXTERNAL)
endif()


# Geogram
if(NOT TARGET geogram::geogram)
	float_tetwild_download_geogram()
	include(geogram)
endif()

# TBB
if(FLOAT_TETWILD_ENABLE_TBB AND NOT TARGET tbb::tbb)
	float_tetwild_download_tbb()

	set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
	set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
	set(TBB_NO_DATE ON CACHE BOOL " " FORCE)

	add_subdirectory(${FLOAT_TETWILD_EXTERNAL}/tbb tbb)
	set_target_properties(tbb_static PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES "${FLOAT_TETWILD_EXTERNAL}/tbb/include"
	)
	if(NOT MSVC)
		set_target_properties(tbb_static PROPERTIES
			COMPILE_FLAGS "-Wno-implicit-fallthrough -Wno-missing-field-initializers -Wno-unused-parameter -Wno-keyword-macro"
		)
		set_target_properties(tbb_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
	endif()
	add_library(tbb::tbb ALIAS tbb_static)
endif()

# C++11 threads
find_package(Threads REQUIRED)
