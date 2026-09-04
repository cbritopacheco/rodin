function(rodin_check_lfs_resources resources_dir)
  if(NOT EXISTS "${resources_dir}")
    return()
  endif()

  file(GLOB_RECURSE _rodin_resource_files
    "${resources_dir}/*.mesh"
    "${resources_dir}/*.msh"
    "${resources_dir}/*.sol"
    "${resources_dir}/*.stp"
    "${resources_dir}/*.step"
    "${resources_dir}/*.gf"
    "${resources_dir}/*.geo"
  )

  foreach(_rodin_resource_file IN LISTS _rodin_resource_files)
    file(READ "${_rodin_resource_file}" _rodin_resource_header LIMIT 64)
    if(_rodin_resource_header MATCHES "^version https://git-lfs.github.com/spec/v1")
      message(FATAL_ERROR
        "Unhydrated Git LFS resource detected: ${_rodin_resource_file}\n"
        "Run `git lfs pull`, or configure with -DRODIN_INSTALL_RESOURCES=OFF "
        "for a library-only install.")
    endif()
  endforeach()
endfunction()
