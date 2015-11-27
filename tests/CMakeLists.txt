foreach(PROG spherical_harmonics spherical_harmonics_iter equilibrium_resistor logtest
             quantity_transfer
             ushape_2d mos1d_dg_n mos1d_dg_p mos1d_potential_kink
             random_numbers markov_chains simple_impurity_scattering
             hde_1d )
   add_executable(${PROG}-test src/${PROG}.cpp )
   target_link_libraries(${PROG}-test shesolvers ${OPENCL_LIBRARIES} viennagrid)
   add_test(${PROG} ${PROG}-test)
endforeach(PROG)

include_directories(${PROJECT_SOURCE_DIR}/external)

add_executable(external_linkage
               src/external_1.cpp
               src/external_2.cpp)

target_link_libraries(external_linkage shesolvers viennagrid)
if (ENABLE_OPENCL)
  target_link_libraries(external_linkage ${OPENCL_LIBRARIES})
endif (ENABLE_OPENCL)

