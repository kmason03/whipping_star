include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${EIGEN3_INCLUDE_DIRS})
#include_directories(/cvmfs/larsoft.opensciencegrid.org/products/root/v6_22_06a/Linux64bit+3.10-2.17-e19-p383b-prof/lib/libMinuit.so)
message("eigen3: ${EIGEN3_INCLUDE_DIRS}")
add_executable (DLAppearenceSens DLAppearenceSens.cc)
add_executable (DLAppearenceSens_OneBinSaveTest DLAppearenceSens_OneBinSaveTest.cc)
add_executable (DL3plus1Sens DL3plus1Sens.cc)
add_executable (DL3plus1_FC DL3plus1_FC.cc)
add_executable (DL3plus1_FCTests DL3plus1_FCTests.cc)
add_executable (DL3plus1_massgen DL3plus1_massgen.cc)


target_link_libraries (DLAppearenceSens SBNfitlib tinyxmllib EventWeight)
target_link_libraries (DLAppearenceSens_OneBinSaveTest SBNfitlib tinyxmllib EventWeight)
target_link_libraries (DL3plus1Sens SBNfitlib tinyxmllib EventWeight)
#extra library added to get TMinuit
target_link_libraries (DL3plus1_FC SBNfitlib tinyxmllib EventWeight /usr/local/root/root-6.22.06/lib/libMinuit.so)
target_link_libraries (DL3plus1_FCTests SBNfitlib tinyxmllib EventWeight /usr/local/root/root-6.22.06/lib/libMinuit.so)
target_link_libraries (DL3plus1_massgen SBNfitlib tinyxmllib EventWeight)
