make_road_dtm = function(dtm, xc)
{
  row <- terra::rowFromY(dtm, xc)
  roadz <- dtm[terra::cellFromRowCol(dtm, row=row)]
  m = matrix(roadz, nrow(dtm), ncol(dtm), byrow = TRUE)
  road_dtm = dtm
  road_dtm[] = m
  return(road_dtm)
}
