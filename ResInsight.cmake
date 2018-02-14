
# ----------------------------------------------------------------------------------------
set(RESINSIGHT_FWK_VIZFWK_LIBRENDER "resinsight/Fwk/VizFwk/LibRender/")

set(LIBRENDER_SOURCE
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfGeometryBuilderDrawableGeo.cpp
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfPrimitiveSetIndexedUInt.cpp
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfScalarMapper.cpp
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfOutlineEdgeExtractor.cpp
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfDrawableGeo.cpp
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfDrawable.cpp
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfHitDetail.cpp
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfPrimitiveSet.cpp
)

set(LIBRENDER_HEADERS
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfGeometryBuilderDrawableGeo.h
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfPrimitiveSetIndexedUInt.h
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfScalarMapper.h
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfOutlineEdgeExtractor.h
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfDrawableGeo.h
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfDrawable.h
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfHitDetail.h
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfPrimitiveSet.h
	${RESINSIGHT_FWK_VIZFWK_LIBRENDER}cvfOpenGLTypes.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_FWK_VIZFWK_LIBCORE "resinsight/Fwk/VizFwk/LibCore/")

set(LIBCORE_SOURCE
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector4.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfString.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfSystem.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector2.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector3.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfCharArray.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfColor3.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfObject.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfAssert.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfDebugTimer.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfAtomicCounter.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfTimer.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfBase64.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfBase64.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfLogger.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfCodeLocation.cpp
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfLogEvent.cpp
)

set(LIBCORE_HEADERS
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector3.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector3.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector4.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector4.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVersion.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfObject.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfPlane.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfString.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfSystem.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfValueArray.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector2.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfVector2.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfBase.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfCharArray.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfCollection.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfCollection.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfColor3.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfColor4.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfConfigCore.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfMath.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfMath.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfMatrix3.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfMatrix3.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfMatrix4.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfMatrix4.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfObject.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfArray.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfArray.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfArrayWrapperConst.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfAssert.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfAtomicCounter.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfTimer.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfLibCore.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfBase64.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfFlags.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfFlags.inl
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfFunctorRange.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfLogger.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfCodeLocation.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfLogDestination.h
${RESINSIGHT_FWK_VIZFWK_LIBCORE}cvfLogEvent.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY "resinsight/Fwk/VizFwk/LibGeometry/")

set(LIBGEOMETRY_SOURCE
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfBoundingBox.cpp
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfGeometryBuilder.cpp
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfOutlineEdgeExtractor.cpp
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfEdgeKey.cpp
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfGeometryUtils.cpp
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfBoundingBoxTree.cpp
)

set(LIBGEOMETRY_HEADERS
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfBoundingBox.h
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfGeometryBuilder.h
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfOutlineEdgeExtractor.h
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfEdgeKey.h
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfGeometryUtils.h
${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}cvfBoundingBoxTree.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL "resinsight/Fwk/AppFwk/cafProjectDataModel/")

set(CAFPROJECTDATAMODEL_HEADERS
${RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL}cafFixedArray.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL_CAFPDMCORE 
"resinsight/Fwk/AppFwk/cafProjectDataModel/cafPdmCore/")

set(CAFPDMCORE_SOURCE
${RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL_CAFPDMCORE}cvfCharArray.cpp
)

set(CAFPDMCORE_HEADERS
${RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL_CAFPDMCORE}cafAppEnum.h
${RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL_CAFPDMCORE}cafAssert.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_FWK_APPFWK_COMMONCODE "resinsight/Fwk/AppFwk/CommonCode/")

set(COMMONCODE_SOURCE
${RESINSIGHT_FWK_APPFWK_COMMONCODE}cvfStructGridGeometryGenerator.cpp
${RESINSIGHT_FWK_APPFWK_COMMONCODE}cvfCellRange.cpp
${RESINSIGHT_FWK_APPFWK_COMMONCODE}cvfStructGrid.cpp
)

set(COMMONCODE_HEADERS	
${RESINSIGHT_FWK_APPFWK_COMMONCODE}cvfStructGridGeometryGenerator.h
${RESINSIGHT_FWK_APPFWK_COMMONCODE}cvfCellRange.h
${RESINSIGHT_FWK_APPFWK_COMMONCODE}cvfStructGrid.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL 
	"resinsight/ApplicationCode/ReservoirDataModel/")

set(RESERVOIRDATAMODEL_SOURCE
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigLocalGrid.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}cvfGeometryTools.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigMainGrid.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigWellPath.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigNNCData.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigWellPathIntersectionTools.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigCell.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigHexIntersectionTools.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigFault.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigGridBase.cpp
${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}RigActiveCellInfo.cpp
)

set(RESERVOIRDATAMODEL_HEADERS
${RESERVOIRDATAMODEL_HEADERS}cvfGeometryTools.h
${RESERVOIRDATAMODEL_HEADERS}cvfGeometryTools.inl
${RESERVOIRDATAMODEL_HEADERS}RigWellPath.h
${RESERVOIRDATAMODEL_HEADERS}RigWellPathIntersectionTools.h
${RESERVOIRDATAMODEL_HEADERS}RigGridBase.h
${RESERVOIRDATAMODEL_HEADERS}RigHexIntersectionTools.h
${RESERVOIRDATAMODEL_HEADERS}RigLocalGrid.h
${RESERVOIRDATAMODEL_HEADERS}RigMainGrid.h
${RESERVOIRDATAMODEL_HEADERS}RigNNCData.h
${RESERVOIRDATAMODEL_HEADERS}RigCell.h
${RESERVOIRDATAMODEL_HEADERS}RigFault.h
${RESERVOIRDATAMODEL_HEADERS}RigActiveCellInfo.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_APPLICATIONCODE_APPLICATION "resinsight/ApplicationCode/Application/")

set(APPLICATION_SOURCE
${RESINSIGHT_APPLICATIONCODE_APPLICATION}RiaDefines.cpp
)

set(APPLICATION_HEADERS
${RESINSIGHT_APPLICATIONCODE_APPLICATION}RiaDefines.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_APPLICATIONCODE_APPLICATION_TOOLS 
	"resinsight/ApplicationCode/Application/Tools/")

set(TOOLS_SOURCE
${RESINSIGHT_APPLICATIONCODE_APPLICATION_TOOLS}RiaLogging.cpp
)

set(TOOLS_HEADERS
${RESINSIGHT_APPLICATIONCODE_APPLICATION_TOOLS}RiaLogging.h
)

# ----------------------------------------------------------------------------------------
set(RESINSIGHT_APPLICATIONCODE_PROJECTDATAMODEL "resinsight/ApplicationCode/ProjectDataModel/")

set(PROJECTDATAMODEL_SOURCE
${RESINSIGHT_APPLICATIONCODE_PROJECTDATAMODEL}RimWellPath.cpp
)

set(PROJECTDATAMODEL_HEADERS
${RESINSIGHT_APPLICATIONCODE_PROJECTDATAMODEL}RimWellPath.h
)