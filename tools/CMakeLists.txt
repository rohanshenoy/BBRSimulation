# Analysis tools: the bbrsim Python package and the plot scripts.
set(bbrsim_PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/python/bbrsim/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/bbrsim/hfss.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/bbrsim/io.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/bbrsim/paths.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/bbrsim/physics.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/bbrsim/select.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/bbrsim/sidecar.py
  )
set(bbrsim_PLOTS
  ${CMAKE_CURRENT_SOURCE_DIR}/plot_crack_angular.py
  ${CMAKE_CURRENT_SOURCE_DIR}/plot_crack_frequency.py
  ${CMAKE_CURRENT_SOURCE_DIR}/plot_cu_reflectance.py
  ${CMAKE_CURRENT_SOURCE_DIR}/plot_test_output.py
  )
install(FILES ${bbrsim_PYTHON} DESTINATION share/BBRsim/python/bbrsim)
install(PROGRAMS ${bbrsim_PLOTS} DESTINATION bin)
