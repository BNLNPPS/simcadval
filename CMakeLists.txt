cmake_minimum_required(VERSION 3.24)
project(simcadval CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(ROOT REQUIRED COMPONENTS Geom)
find_package(OpenCASCADE REQUIRED)
find_package(Geant4 REQUIRED)
find_package(G4OCCT REQUIRED)

# DD4hep compact -> STEP + materials.xml. ROOT and OCCT stay in separate translation
# units (geo_desc.hh is the std-only interface between them).
add_executable(dd4hep2step src/dd4hep2step.cc src/tgeo_reader.cc src/occt_writer.cc)
target_include_directories(dd4hep2step PRIVATE src ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(dd4hep2step PRIVATE DD4hep::DDCore ROOT::Geom
  TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo TKPrim TKBO
  TKBool TKOffset TKShHealing TKLCAF TKCAF TKCDF TKXCAF TKXSBase TKDE TKDESTEP)

# Geantino ray scan + point sampling of one Geant4 world (DD4hep or STEP via G4OCCT).
add_executable(simcadval_scan src/simcadval_scan.cc src/step_world.cc)
target_include_directories(simcadval_scan PRIVATE src)
target_link_libraries(simcadval_scan PRIVATE G4OCCT::G4OCCT DD4hep::DDCore DD4hep::DDG4 ROOT::Geom
  ${Geant4_LIBRARIES})

install(TARGETS dd4hep2step simcadval_scan RUNTIME DESTINATION bin)
