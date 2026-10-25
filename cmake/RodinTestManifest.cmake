# Export configured targets rather than inferring executables from CMake source
# text. In particular, this includes optional backends and newly added tests.
function(rodin_collect_test_targets directory)
  get_property(_targets DIRECTORY "${directory}" PROPERTY BUILDSYSTEM_TARGETS)
  foreach(_target IN LISTS _targets)
    get_target_property(_type "${_target}" TYPE)
    get_target_property(_source "${_target}" SOURCE_DIR)
    file(RELATIVE_PATH _owner "${PROJECT_SOURCE_DIR}" "${_source}")
    if(_owner MATCHES "^tests/" AND NOT _owner MATCHES "^tests/benchmarks(/|$)"
        AND _type STREQUAL "EXECUTABLE")
      get_property(_profiles TARGET "${_target}" PROPERTY RODIN_TEST_PROFILES)
      set_property(GLOBAL APPEND_STRING PROPERTY RODIN_TEST_MANIFEST
        "${_target}\t${_owner}\t$<TARGET_FILE:${_target}>\ttest\t${_profiles}\n")
    elseif(_owner MATCHES "^src/Rodin($|/)" AND
        _type MATCHES "^(STATIC|SHARED|OBJECT)_LIBRARY$")
      # These are the production compilation targets in the old default build.
      # The baseline builds them even if no selected test links a given module.
      set_property(GLOBAL APPEND_STRING PROPERTY RODIN_TEST_MANIFEST
        "${_target}\t${_owner}\t\tlibrary\t\n")
    endif()
  endforeach()
  get_property(_children DIRECTORY "${directory}" PROPERTY SUBDIRECTORIES)
  foreach(_child IN LISTS _children)
    rodin_collect_test_targets("${_child}")
  endforeach()
endfunction()

rodin_collect_test_targets("${PROJECT_SOURCE_DIR}")
get_property(_manifest GLOBAL PROPERTY RODIN_TEST_MANIFEST)
file(GENERATE OUTPUT "${PROJECT_BINARY_DIR}/RodinTestTargets-$<CONFIG>.tsv"
  CONTENT "target\towner\texecutable\tkind\tprofiles\n${_manifest}")

# Let the build tool refresh CMake's configuration before the runner reads the
# manifest. This target compiles nothing and prevents stale target selection.
add_custom_target(RodinTestManifest)
