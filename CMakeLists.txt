cmake_minimum_required(VERSION 3.18)
project(fpocket-python LANGUAGES C CXX)

# Detect platform for molfile plugin path
if(APPLE)
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        set(FPOCKET_ARCH "MACOSXARM64")
    else()
        set(FPOCKET_ARCH "MACOSXX86_64")
    endif()
    set(FPOCKET_OS_FLAG "M_OS_MAC")
elseif(UNIX)
    set(FPOCKET_ARCH "LINUXAMD64")
    set(FPOCKET_OS_FLAG "M_OS_LINUX")
endif()

set(MOLFILE_DIR "${CMAKE_SOURCE_DIR}/plugins/${FPOCKET_ARCH}/molfile")
set(PLUGIN_INCLUDE "${CMAKE_SOURCE_DIR}/plugins/include")

# Qhull sources
file(GLOB QHULL_SOURCES
    src/qhull/src/libqhull/geom2.c
    src/qhull/src/libqhull/geom.c
    src/qhull/src/libqhull/global.c
    src/qhull/src/libqhull/io.c
    src/qhull/src/libqhull/libqhull.c
    src/qhull/src/libqhull/mem.c
    src/qhull/src/libqhull/merge.c
    src/qhull/src/libqhull/poly2.c
    src/qhull/src/libqhull/poly.c
    src/qhull/src/libqhull/qset.c
    src/qhull/src/libqhull/random.c
    src/qhull/src/libqhull/rboxlib.c
    src/qhull/src/libqhull/stat.c
    src/qhull/src/libqhull/user.c
    src/qhull/src/libqhull/usermem.c
    src/qhull/src/libqhull/userprintf.c
    src/qhull/src/libqhull/userprintf_rbox.c
    src/qhull/src/qvoronoi/qvoronoi.c
    src/qhull/src/qconvex/qconvex.c
)

# Fpocket sources (exclude main entry points)
set(FPOCKET_SOURCES
    src/aa.c
    src/asa.c
    src/atom.c
    src/calc.c
    src/cluster.c
    src/clusterlib.c
    src/descriptors.c
    src/dparams.c
    src/dpocket.c
    src/energy.c
    src/fparams.c
    src/fpocket.c
    src/fpout.c
    src/grid.c
    src/mdparams.c
    src/mdpbase.c
    src/mdpmem.c
    src/mdpocket.c
    src/mdpout.c
    src/memhandler.c
    src/neighbor.c
    src/pertable.c
    src/pocket.c
    src/pscoring.c
    src/psorting.c
    src/read_mmcif.c
    src/refine.c
    src/rpdb.c
    src/sort.c
    src/topology.c
    src/tparams.c
    src/tpocket.c
    src/utils.c
    src/voronoi.c
    src/voronoi_lst.c
    src/write_visu.c
    src/writepdb.c
    src/writepocket.c
)

# Build fpocket as a static library
add_library(fpocket_core STATIC ${FPOCKET_SOURCES} ${QHULL_SOURCES})
target_include_directories(fpocket_core PUBLIC
    ${CMAKE_SOURCE_DIR}/headers
    ${CMAKE_SOURCE_DIR}/src/qhull/src
    ${PLUGIN_INCLUDE}
    ${MOLFILE_DIR}
)
target_compile_definitions(fpocket_core PRIVATE
    ${FPOCKET_OS_FLAG}
    MNO_MEM_DEBUG
    MD_NOT_USE_GSL
)
target_compile_options(fpocket_core PRIVATE -O2 -std=gnu99 -Wno-incompatible-pointer-types -fPIC)

# Build the C wrapper as a shared library for CFFI
add_library(fpocket_wrapper SHARED python/fpocket_wrapper.c)
target_link_libraries(fpocket_wrapper PRIVATE
    fpocket_core
    ${MOLFILE_DIR}/libmolfile_plugin.a
)
target_include_directories(fpocket_wrapper PRIVATE
    ${CMAKE_SOURCE_DIR}/headers
    ${CMAKE_SOURCE_DIR}/src/qhull/src
    ${PLUGIN_INCLUDE}
    ${MOLFILE_DIR}
)
target_compile_definitions(fpocket_wrapper PRIVATE
    ${FPOCKET_OS_FLAG}
    MNO_MEM_DEBUG
    MD_NOT_USE_GSL
)
target_link_libraries(fpocket_wrapper PRIVATE m stdc++)

if(APPLE)
    set_target_properties(fpocket_wrapper PROPERTIES
        INSTALL_RPATH "@loader_path"
        MACOSX_RPATH TRUE
    )
else()
    set_target_properties(fpocket_wrapper PROPERTIES
        INSTALL_RPATH "$ORIGIN"
    )
endif()

install(TARGETS fpocket_wrapper LIBRARY DESTINATION .)
