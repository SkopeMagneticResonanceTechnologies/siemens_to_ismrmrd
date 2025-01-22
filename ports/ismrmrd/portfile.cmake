if (VCPKG_TARGET_ARCHITECTURE MATCHES "x86")
    set(WIN32_INCLUDE_STDDEF_PATCH "x86-windows-include-stddef.patch")
endif()

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO SkopeMagneticResonanceTechnologies/ismrmrd
    REF 1a50114e524deeda14b71a1c4e9ec00ecd767d96
    SHA512 07094f79ff9b0d6aa5e91596c113c56096b6685b2332384890c4906f2b4ade82d40aa296a66b58546740a9ba2dae19fbfd0d005d6fca198977f70aaa57d6c518
    HEAD_REF static_linkage
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DUSE_SYSTEM_PUGIXML=ON
        -DUSE_HDF5_DATASET_SUPPORT=ON
        -DVCPKG_TARGET_TRIPLET=ON
        -DBUILD_TESTS=OFF
        -DBUILD_EXAMPLES=OFF
        -DBUILD_UTILITIES=OFF
)

vcpkg_cmake_install()
vcpkg_copy_pdbs()
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/ISMRMRD/)

if(EXISTS "${CURRENT_PACKAGES_DIR}/lib/ismrmrd.dll")
    file(COPY "${CURRENT_PACKAGES_DIR}/lib/ismrmrd.dll" DESTINATION "${CURRENT_PACKAGES_DIR}/bin")
    file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/lib/ismrmrd.dll")
endif()

if(EXISTS "${CURRENT_PACKAGES_DIR}/debug/lib/ismrmrd.dll")
    file(COPY "${CURRENT_PACKAGES_DIR}/debug/lib/ismrmrd.dll" DESTINATION "${CURRENT_PACKAGES_DIR}/debug/bin")
    file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/lib/ismrmrd.dll")
endif()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/share/ismrmrd/cmake")

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

if(VCPKG_LIBRARY_LINKAGE STREQUAL static)
    file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/bin/")
    file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/bin/")
endif()

file(INSTALL "${SOURCE_PATH}/LICENSE" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)
