
blt_add_library(
    NAME simdops
    HEADERS
       include/simdops/null.hpp
       include/simdops/power_qpx.hpp
       include/simdops/resetArch.hpp
       include/simdops/simdops.hpp
       include/simdops/x86_avx2.hpp
       include/simdops/x86_avx512f.hpp
       include/simdops/default_math.hpp
)
target_include_directories(simdops INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

