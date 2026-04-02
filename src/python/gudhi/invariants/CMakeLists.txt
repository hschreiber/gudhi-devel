nanobind_add_module_with_targets(_persistence_objects_ext _persistence_objects.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_pers_cohomology_ext _helpers.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_build_helpers_ext _helpers.cc TARGETS TBB::tbb)
