# load data
contours25ld <- shapefile(paste(dem.dir, "contour-lines25_ld.shp", sep = ""))
contours25ld$level <- as.numeric(contours25ld$level)
contours25md <- shapefile(paste(dem.dir, "contour-lines25_md.shp", sep = ""))
contours25md$level <- as.numeric(contours25md$level)
contours25hd <- shapefile(paste(dem.dir, "contour-lines25_hd.shp", sep = ""))
contours25hd$level <- as.numeric(contours25hd$level)

# transform to WGS84
contours25ld <- spTransform(contours25ld, wgs1984utm22s)
contours25md <- spTransform(contours25md, wgs1984utm22s)
contours25hd <- spTransform(contours25hd, wgs1984utm22s)

# apply affine transformation
contours25ld <- applyTransformation(affine.car25, contours25ld)
contours25md <- applyTransformation(affine.car25, contours25md)
contours25hd <- applyTransformation(affine.car25, contours25hd)

# GRASS GIS mapset = "carta25"
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, location = "dnos-sm-rs",
          mapset = "carta25", pid = Sys.getpid(), override = TRUE)
writeRAST6(grid.carta25, "grid.carta25", flags = "overwrite")
system("g.region rast=grid.carta25");gmeta6()

# export to GRASS GIS
writeVECT6(SDF = buffer25, "buffer25", v.in.ogr_flags = c("overwrite"))
system("v.info -c buffer25")
system("d.mon x0"); system("d.vect buffer25")
writeVECT6(SDF = contours25, "contours25", v.in.ogr_flags = c("overwrite"))
system("v.info -c contours25")
system(paste("d.mon x0")); system("d.vect contours25")
writeVECT6(SDF = contours25ld, "contours25ld", v.in.ogr_flags = c("overwrite"))
system("v.info -c contours25ld")
system(paste("d.mon x0")); system("d.vect contours25ld")
writeVECT6(SDF = contours25md, "contours25md", v.in.ogr_flags = c("overwrite"))
system("v.info -c contours25md")
system(paste("d.mon x0")); system("d.vect contours25md")
writeVECT6(SDF = contours25hd, "contours25hd", v.in.ogr_flags = c("overwrite"))
system("v.info -c contours25hd")
system(paste("d.mon x0")); system("d.vect contours25hd")
writeVECT6(SDF = points25, "points25", v.in.ogr_flags = c("overwrite"))
system("v.info -c points25")
system("d.mon x0"); system("d.vect points25")
writeVECT6(SDF = rivers25, "rivers25", v.in.ogr_flags = c("overwrite"))
system("v.info rivers25")
system("d.vect rivers25")
writeVECT6(SDF = lakes25, "lakes25", v.in.ogr_flags = c("overwrite"))
system("v.info lakes25")
system("d.mon x0");system("d.vect lakes25")

# contour lines to points
system("v.to.points -v --o type=line dmax=20 in=contours25 out=contours25pts")
system("v.info -c contours25pts")
system("v.to.points -v --o type=line dmax=20 in=contours25ld out=contours25ldpts")
system("v.to.points -v --o type=line dmax=20 in=contours25md out=contours25mdpts")
system("v.to.points -v --o type=line dmax=20 in=contours25hd out=contours25hdpts")

# merge point data
system("v.patch --o -e input=contours25pts,points25 output=level25pts")
system("v.info level25pts")
system("v.info -c level25pts")
system(paste("d.mon x0"));system("d.vect level25pts")

# transform 2D points to 3D
system("v.to.3d --o in=level25pts output=level25pts3d type=point layer=1 column=level")
system("v.info level25pts3d")
system(paste("d.mon x0")); system("d.vect level25pts3d -z zcolor=haxby")
system("v.to.3d --o in=contours25ldpts output=contours25ld3d type=point layer=1 column=level")
system("v.to.3d --o in=contours25mdpts output=contours25md3d type=point layer=1 column=level")
system("v.to.3d --o in=contours25hdpts output=contours25hd3d type=point layer=1 column=level")

# evaluate interpolation methods: v.surf.bspline
system("v.surf.bspline -c -e --v input=contours25ld3d method=bicubic")
system("v.surf.bspline -c sie=100, sin=100 --v input=contours25ld3d method=bicubic")
# Table of results: v.surf.bspline complete. Cross validation finished
# for sie = 100.000000 and sin = 100.000000
#    lambda |       mean |        rms |
#   0.01000 |    -0.0010 |     1.1011 |
#   0.05000 |    -0.0014 |     1.2166 |
#   0.10000 |    -0.0012 |     1.3432 |
#   0.20000 |    -0.0005 |     1.5159 |
#   0.30000 |     0.0001 |     1.6274 |
#   0.40000 |     0.0008 |     1.7082 |
system("g.region vect=contours25ld")
system("v.surf.bspline --o sie=100 sin=100 input=contours25ld3d raster=dem.bspline method=bicubic")
system("nviz elev=dem.bspline vector=contours25ld")
system("g.region vect=contours25md")
system("v.surf.bspline --o sie=100 sin=100 input=contours25md3d raster=dem.bspline method=bicubic")
system("nviz elev=dem.bspline vector=contours25md")
system("g.region vect=contours25hd")
system("v.surf.bspline --o sie=100 sin=100 input=contours25hd3d raster=dem.bspline method=bicubic")
system("nviz elev=dem.bspline vector=contours25hd")

# evaluate interpolation methods: v.surf.rst
system("g.region vect=contours25ld")
system("v.surf.rst --o input=contours25ld3d layer=0 elev=dem.rst smooth=0.01 dmin=0.01")
system("nviz elev=dem.rst vector=contours25ld")
system("g.region vect=contours25md")
system("v.surf.rst --o input=contours25md3d layer=0 elev=dem.rst smooth=0.01 dmin=0.01")
system("nviz elev=dem.rst vector=contours25md")
system("g.region vect=contours25hd")
system("v.surf.rst --o input=contours25hd3d layer=0 elev=dem.rst smooth=0.01 dmin=0.01")
system("nviz elev=dem.rst vector=contours25hd")
system("g.remove rast=dem.bspline,dem.rst")

# CONCLUSION: interpolate using v.surf.bspline.
# interpolate using v.surf.bspline
system("v.surf.bspline --o --v layer=0 lambda_i=0.01 sie=150 sin=150 input=level25pts3d raster=dem.carta25 method=bicubic")
system("r.info dem.carta25")
system("nviz elev=dem.carta25 vector=contours25")

# fill sinks
system("r.fill.dir --o input=dem.carta25 elevation=dem.carta25.fill direction=FDIR_10 type=answers")
system("r.info dem.carta25.fill")
system("r.mapcalc DIFF_10=dem.carta25.fill-dem.carta25");system("r.info DIFF_10")
system("d.mon x0");system("d.rast DIFF_10")
system("d.mon x0");system("d.rast FDIR_10")
system("nviz elev=dem.carta25.fill vector=contours25")

# burn streams
system("r.carve")
system("r.carve --o rast=dem.carta25.fill vect=rivers25 output=dem.carta25.fill.burn points=dem_carta25_burn_pts depth=10 width=10")

# difference
system("r.mapcalc 'tmp=dem.carta25.fill-dem.carta25.fill.burn'")
system("r.null tmp setnull=0")
system("r.info tmp")
system("r.colors map=tmp color=grey")
system("d.mon x0");system("d.rast.leg tmp")

# raster to vector
system("r.thin in=tmp out=tmp.thin")

system("r.to.vect in=tmp out=tmp feature=point")
system("v.info tmp")
system("d.mon x0");system("d.vect tmp")

# use stream points to create new DEM
system("v.patch --o input=level25pts3d,dem_carta25_burn_pts output=level25pts3d_burn")
system("v.info level25pts3d_burn")
system(paste("d.mon x0")); system("d.vect level25pts3d_burn")

# interpolate new DEm using v.surf.bspline
system("v.surf.bspline --o --v layer=0 lambda_i=0.01 sie=150 sin=150 input=level25pts3d_burn raster=dem.carta25.burn method=bicubic")
system("r.info dem.carta25.burn")
system("r.colors dem.carta25.burn color=haxby")
system("nviz elev=dem.carta25.burn")

# fill sinks of new DEM
system("r.fill.dir --o input=dem.carta25.burn elevation=dem.carta25.fill direction=FDIR_10 type=answers")
system("r.info dem.carta25.fill")
system("r.mapcalc DIFF_10=dem.carta25.fill-dem.carta25.burn");system("r.info DIFF_10")
system("d.mon x0");system("d.rast DIFF_10")
system("d.mon x0");system("d.rast FDIR_10")
system("nviz elev=dem.carta25.fill vector=contours25")

# import into GRASS GIS mapset "predicitions"
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, override = TRUE,
          location = "dnos-sm-rs", mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster");gmeta6()
system("g.copy --o rast=MASK@carta25,MASK")
system("g.copy --o rast=dem.carta25.egm@carta25,dem.carta25.egm")
system("r.mapcalc 'ELEV_90=dem.carta25.egm'")
system("r.info ELEV_90")