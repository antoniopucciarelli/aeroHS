cmake_minimum_required(VERSION 3.16.0)

enable_language(Fortran)

project(aeroHS_PRJ)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

find_package(LAPACK REQUIRED)
list(APPEND linkedlibs ${LAPACK_LIBRARIES})

find_package(Gnuplot REQUIRED)

file(GLOB_RECURSE source_files ${PROJECT_SOURCE_DIR}/src/*.f90)
file(GLOB_RECURSE module_files ${PROJECT_SOURCE_DIR}/src/modules/*.f90)

add_executable(aeroHS ${source_files} ${module_files})
target_link_libraries(aeroHS ${linkedlibs})

# coping .plt files
file(COPY ${PROJECT_SOURCE_DIR}/src/plt/CP_plot.plt             DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(COPY ${PROJECT_SOURCE_DIR}/src/plt/VELfield.plt            DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(COPY ${PROJECT_SOURCE_DIR}/src/plt/velPLOT.py             DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(COPY ${PROJECT_SOURCE_DIR}/src/plt/CL_plot.plt             DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(COPY ${PROJECT_SOURCE_DIR}/src/plt/vel_plot.plt            DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(COPY ${PROJECT_SOURCE_DIR}/src/plt/pressure_plot.plt       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(COPY ${PROJECT_SOURCE_DIR}/src/plt/AIRFOILgnuplot.plt      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(COPY ${PROJECT_SOURCE_DIR}/src/plt/AIRFOILelement_plot.plt DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
