
# RESINXX DIRS =========================================================
set(RESINXX resinxx) 
set(RESINXX_GRID      ${RESINXX}/rixx_grid)
set(RESINXX_APP_FWK   ${RESINXX}/rixx_app_fwk)
set(RESINXX_CORE_GEOM ${RESINXX}/rixx_core_geom)
set(RESINXX_RES_MOD   ${RESINXX}/rixx_res_mod)

# MAIN RESINXX FILES ===================================================
set(RIXX_CPP_FILES
${RESINXX}/well_path.cpp
${RESINXX}/geometry_tools.cpp
)

# GRID FILES ===========================================================
set(RIXX_GRID_CPP_FILES
${RESINXX_GRID}/riextractor.cpp
${RESINXX_GRID}/rifaultncc.cpp
${RESINXX_GRID}/ricasedata.cpp
${RESINXX_GRID}/rigrid.cpp
${RESINXX_GRID}/ricell.cpp
)

# APP FWK FILES ========================================================
set(RIXX_APP_FWK_CPP_FILES
${RESINXX_APP_FWK}/cvfStructGrid.cpp
${RESINXX_APP_FWK}/cvfCellRange.cpp
${RESINXX_APP_FWK}/cafHexGridIntersectionTools.cpp
)

# CORE GEOM FILES ======================================================
set(RIXX_CORE_GEOM_CPP_FILES
${RESINXX_CORE_GEOM}/cvfAssert.cpp
${RESINXX_CORE_GEOM}/cvfAtomicCounter.cpp
${RESINXX_CORE_GEOM}/cvfBoundingBox.cpp
${RESINXX_CORE_GEOM}/cvfBoundingBoxTree.cpp
${RESINXX_CORE_GEOM}/cvfCharArray.cpp
${RESINXX_CORE_GEOM}/cvfMath.cpp
${RESINXX_CORE_GEOM}/cvfPlane.cpp
${RESINXX_CORE_GEOM}/cvfObject.cpp
${RESINXX_CORE_GEOM}/cvfRay.cpp
${RESINXX_CORE_GEOM}/cvfString.cpp
${RESINXX_CORE_GEOM}/cvfSystem.cpp
${RESINXX_CORE_GEOM}/cvfVector2.cpp
${RESINXX_CORE_GEOM}/cvfVector3.cpp
${RESINXX_CORE_GEOM}/cvfVector4.cpp
)

# RES MOD FILES ========================================================
set(RIXX_RES_MOD_CPP_FILES
${RESINXX_RES_MOD}/cvfGeometryTools.cpp
${RESINXX_RES_MOD}/RigCellGeometryTools.cpp
)

message(".............................................................")
message("RIXX_CPP_FILES: ${RIXX_CPP_FILES}")
message(".............................................................")
message("RIXX_GRID_CPP_FILES: ${RIXX_GRID_CPP_FILES}")
message(".............................................................")
message("RIXX_APP_FWK_CPP_FILES: ${RIXX_APP_FWK_CPP_FILES}")
message(".............................................................")
message("RIXX_CORE_GEOM_CPP_FILES: ${RIXX_CORE_GEOM_CPP_FILES}")
message(".............................................................")
message("RIXX_RES_MOD_CPP_FILES: ${RIXX_RES_MOD_CPP_FILES}")

# ALL ==================================================================
set(RIXX_ALL_CPP_FILES
${RIXX_CPP_FILES}
${RIXX_GRID_CPP_FILES}
${RIXX_APP_FWK_CPP_FILES}
${RIXX_CORE_GEOM_CPP_FILES}
${RIXX_RES_MOD_CPP_FILES}
)