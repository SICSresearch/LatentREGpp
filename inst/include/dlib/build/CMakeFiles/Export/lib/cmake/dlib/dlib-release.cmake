#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "dlib::dlib" for configuration "Release"
set_property(TARGET dlib::dlib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(dlib::dlib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/usr/lib64/libpthread.so;/usr/lib64/libnsl.so;/usr/lib64/libX11.so;/usr/lib64/libpng.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdlib.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS dlib::dlib )
list(APPEND _IMPORT_CHECK_FILES_FOR_dlib::dlib "${_IMPORT_PREFIX}/lib/libdlib.a" )

# Import target "dlib::dlib-shared" for configuration "Release"
set_property(TARGET dlib::dlib-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(dlib::dlib-shared PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/usr/lib64/libpthread.so;/usr/lib64/libnsl.so;/usr/lib64/libX11.so;/usr/lib64/libpng.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdlib.so.18.18.0"
  IMPORTED_SONAME_RELEASE "libdlib.so.18.18.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS dlib::dlib-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_dlib::dlib-shared "${_IMPORT_PREFIX}/lib/libdlib.so.18.18.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
