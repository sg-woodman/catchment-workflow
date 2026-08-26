///////////////////////////////////////////////////////////////////////////////////////
var yrange    = [1984,2025]; // year; range 1985 to present
var mrange    = [5,9]; // month; range 1 to 12
var maxcloud  = 20; // max threshold for cloud cover
var buffer    = 10; // units pixels
var scale     = 30; // units m
var crs       = 'EPSG:3161';
var roi       = aoi; // polygons from FLImports
var roi_bbox = roi.bounds()
///////////////////////////////////////////////////////////////////////////////////////

function maskL457sr(image) {
    // Bit 0 - Fill
    // Bit 1 - Dilated Cloud
    // Bit 2 - Unused
    // Bit 3 - Cloud
    // Bit 4 - Cloud Shadow
    var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('111111', 2)).eq(0);
    var saturationMask = image.select('QA_RADSAT').eq(0);

    // Apply the scaling factors to the appropriate bands.
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);

    // Replace the original bands with the scaled ones and apply the masks.
    return image.addBands(opticalBands, null, true)
        .addBands(thermalBand, null, true)
        .updateMask(qaMask)
        .updateMask(saturationMask);

}

function maskL89sr(image) {
    // Bit 0 - Fill
    // Bit 1 - Dilated Cloud
    // Bit 2 - Cirrus
    // Bit 3 - Cloud
    // Bit 4 - Cloud Shadow
    var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('111111', 2)).eq(0);
    var saturationMask = image.select('QA_RADSAT').eq(0);

    // Apply the scaling factors to the appropriate bands.
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);

    // Replace the original bands with the scaled ones and apply the masks.
    return image.addBands(opticalBands, null, true)
        .addBands(thermalBands, null, true)
        .updateMask(qaMask)
        .updateMask(saturationMask);

}

var get_indices = function(image){
/*
CALCULATES GREENNESS AND MOISTURE INDICES
GREENNESS:
EVI : ENHANCED VEGETATION INDEX
NDVI : NORMALIZED DIFFERENCE VEGETATION INDEX
NIRV : NEAR INFRARED REFLECTANCE OF VEGETATION
GCC : GREEN CHROMATIC COORDINATE (PERCENT GREENNESS)

MOISTURE/WATER:
NDMI : NORMALIZED DIFFERENCE MOISTURE INDEX
NDWI : NORMALIZED DIFFERENCE WATER INDEX (Gao 1996)
MNDWI : MODIFIED NORMALIZED DIFFERENCE WATER INDEX (Xu 2006)
NBR : NORMALIZED BURN RATIO (also sensitive to moisture)
LSWI : LAND SURFACE WATER INDEX

NOTE: Most indices are scaled by a factor of 10000.
*/

// Greenness indices
var EVI = image.expression('2.5 * (NIR_median-R_median) / (NIR_median + 6 * R_median - 7.5* B_median + 10000)',
    {NIR_median : image.select('NIR_median'),
    R_median   : image.select('R_median'),
    B_median   : image.select('B_median')}
    ),
NDVI = image.expression('(NIR_median-R_median) / (NIR_median+R_median)',
    {NIR_median : image.select('NIR_median'),
    R_median   : image.select('R_median')}
    ),
NIRV = image.expression('(NIR_median-R_median) / (NIR_median+R_median) * NIR_median',
    {NIR_median : image.select('NIR_median'),
    R_median   : image.select('R_median')}
    ),
GCC = image.expression('G_median / (B_median + G_median + R_median)',
    {R_median : image.select('R_median'),
    G_median : image.select('G_median'),
    B_median : image.select('B_median')}
    );

// Moisture/Water indices
var NDMI = image.expression('(NIR_median - SWIR1_median) / (NIR_median + SWIR1_median)',
    {NIR_median : image.select('NIR_median'),
    SWIR1_median : image.select('SWIR1_median')}
    ),
NDWI = image.expression('(NIR_median - SWIR2_median) / (NIR_median + SWIR2_median)',
    {NIR_median : image.select('NIR_median'),
    SWIR2_median : image.select('SWIR2_median')}
    ),
MNDWI = image.expression('(G_median - SWIR1_median) / (G_median + SWIR1_median)',
    {G_median : image.select('G_median'),
    SWIR1_median : image.select('SWIR1_median')}
    ),
NBR = image.expression('(NIR_median - SWIR2_median) / (NIR_median + SWIR2_median)',
    {NIR_median : image.select('NIR_median'),
    SWIR2_median : image.select('SWIR2_median')}
    ),
LSWI = image.expression('(NIR_median - SWIR1_median) / (NIR_median + SWIR1_median)',
    {NIR_median : image.select('NIR_median'),
    SWIR1_median : image.select('SWIR1_median')}
    );

var newbands = ee.Image.cat([EVI,NDVI,GCC,NDMI,NDWI,MNDWI,NBR,LSWI])
    .multiply(10000)
    .addBands(NIRV)
    .round()
    .short()
    .rename(['EVI','NDVI','GCC','NDMI','NDWI','MNDWI','NBR','LSWI','NIRV']);
return image.addBands(newbands);
};

var bandslsat89 = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL','QA_RADSAT'];
var bandslsat457 = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','QA_PIXEL','QA_RADSAT'];
var bandsout = ['B','G','R','NIR','SWIR1','SWIR2','pixel_qa','radsat_qa'];

var lsat4 = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2')
    .filter(ee.Filter.lt('CLOUD_COVER',maxcloud))
    .filterBounds(roi)
    .map(maskL457sr)
    .select(bandslsat457,bandsout);

var lsat5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filter(ee.Filter.lt('CLOUD_COVER',maxcloud))
    .filterBounds(roi)
    .map(maskL457sr)
    .select(bandslsat457,bandsout);

var lsat7_raw = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filter(ee.Filter.lt('CLOUD_COVER',maxcloud))
    .filterBounds(roi)
    .map(maskL457sr)
    .select(bandslsat457,bandsout);

print(lsat7_raw)

var fillGaps = function(image) {
  var filled = image.focal_mean(1.5, 'square', 'pixels', 10)

  var out = filled.blend(image)

  return out.set('system:time_start',image.get('system:time_start'))
}

var castl7 = function(image){
  // Define the bands you want to cast and their new types
    var l7bandTypes = {
      pixel_qa: 'uint16',
      radsat_qa: 'uint16'
    };

  var recast = image.cast(l7bandTypes)

  return recast
}

var lsat7 = lsat7_raw.map(fillGaps).map(castl7)


print(lsat7)

var lsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filter(ee.Filter.lt('CLOUD_COVER',maxcloud))
    .filterBounds(roi)
    .map(maskL89sr)
    .select(bandslsat89,bandsout);
print(lsat8)
var lsat9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filter(ee.Filter.lt('CLOUD_COVER',maxcloud))
    .filterBounds(roi)
    .map(maskL89sr)
    .select(bandslsat89,bandsout);

var col = lsat4.merge(lsat5).merge(lsat7).merge(lsat8).merge(lsat9)
  .filter(ee.Filter.calendarRange(yrange[0],yrange[1],'year'))
  .filter(ee.Filter.calendarRange(mrange[0],mrange[1],'month'))

print(col.limit(50), 'Full Collection')

var years = ee.List.sequence(yrange[0],yrange[1])
print (years)

// FIX: reduce(ee.Reducer.median()) on an EMPTY ImageCollection returns an image
// with ZERO bands (there's nothing for the reducer to infer band names from).
// For AOIs/years where no scene survives the cloud/bounds/date filters, that
// 0-band image then gets passed into get_indices(), which unconditionally
// calls .select('NIR_median') etc. -- and THAT is what throws the
// "Band pattern '...' was applied to an Image with no bands" error.
// The fix: build a fully-masked placeholder image with the exact band names
// the reducer would have produced, so every yearly composite -- even empty
// ones -- has an identical, valid 8-band structure before get_indices runs.
var medianBandNames = bandsout.map(function(b){ return b + '_median'; });
var emptyYearTemplate = ee.Image.constant(ee.List.repeat(0, medianBandNames.length))
    .rename(medianBandNames)
    .selfMask(); // fully masked -> "no valid data", but band names/count are intact

// Create yearly composites with proper band naming
var collectYear = ee.ImageCollection.fromImages(
    years.map(function(y){
        var start = ee.Date.fromYMD(y, 1, 1)
        var end = start.advance(12, 'month');
        var yearCol = col.filterDate(start, end);
        var nImages = yearCol.size(); // scenes actually contributing to this year

        // Use the real median composite when scenes exist; otherwise fall back
        // to the masked template so downstream band-dependent code never sees
        // a bandless image.
        var yearlyComposite = ee.Image(ee.Algorithms.If(
            nImages.gt(0),
            yearCol.reduce(ee.Reducer.median()),
            emptyYearTemplate
        ));

        return yearlyComposite
            .set('year', y)
            // nImages lets us filter out empty years cheaply later, without
            // needing to evaluate the whole band-dependent graph (bandNames()
            // triggers exactly the crash we're fixing, since it forces
            // evaluation of get_indices' selects).
            .set('nImages', nImages)
            .set('system:time_start', start.millis());
    }))
    .map(get_indices);

print(collectYear, 'Yearly Collection')

/**

 * Converts specified Landsat bands to int16 data type
 * @param {ee.Image} image - Input Landsat image
 * @param {Array<string>} bandsToConvert - Array of band names to convert to int16
 * @param {number} scaleFactor - Scale factor to apply before conversion (default: 10000)
 * @returns {ee.Image} - Image with specified bands converted to int16
   */
   function convertLandsatBandsToInt16(image, bandsToConvert, scaleFactor) {
     // Set default scale factor if not provided
     scaleFactor = scaleFactor || 10000;

  // Get all band names from the image
  var allBands = image.bandNames();

  // Filter to only include bands that exist in the image
  var validBands = ee.List(bandsToConvert).filter(ee.Filter.inList('item', allBands));

  // Convert specified bands to int16
  var convertedBands = validBands.map(function(bandName) {
    // Cast bandName to string to ensure proper type
    bandName = ee.String(bandName);
    return image.select([bandName])
                .multiply(scaleFactor)
                .int16()
                .rename([bandName]); // Pass as list
  });

  // Get bands that should remain unchanged
  var unchangedBandNames = allBands.removeAll(validBands);
  var unchangedBands = image.select(unchangedBandNames);

  // Combine converted and unchanged bands
  var result = ee.Image(ee.Algorithms.If(
    ee.List(convertedBands).size().gt(0),
    // Use iterate to combine all converted bands into a single image
    ee.Image(ee.List(convertedBands).iterate(function(img, acc) {
      return ee.Image(acc).addBands(ee.Image(img));
    }, ee.Image([]))).addBands(unchangedBands),
    image
  ));

  // Preserve original image properties
  return result.copyProperties(image, image.propertyNames());
}

// Example usage function for Landsat 4/5/7 surface reflectance
function convertLandsatBands(image) {
  var bandsToConvert = [
    'B_median', 'G_median', 'R_median', 'NIR_median',
    'SWIR1_median', 'SWIR2_median', 'pixel_qa_median', 'radsat_qa_median'
  ];

  return convertLandsatBandsToInt16(image, bandsToConvert, 10000);
}

// Filter out years with no contributing scenes (i.e., the masked template
// years from the empty-year fallback above). We filter on the 'nImages'
// property set earlier rather than on bandNames().length(): every composite
// now always has 17 bands (8 spectral/QA + 9 indices), even empty ones, so a
// band-count check can no longer distinguish real years from empty ones --
// and worse, computing bandNames() on a genuinely broken image is what
// crashed the original pipeline.
var nullimages = collectYear
    .filter(ee.Filter.gt('nImages', 0))
    .map(convertLandsatBands)
    .map(function(image){return image.clip(roi_bbox)});
print(nullimages, 'Complete images (years with at least 1 contributing scene)')

// Optional diagnostic: see exactly which years had zero contributing scenes
// for this AOI, so you know which years will be missing from the exported
// time series (rather than silently dropped).
var emptyYears = collectYear
    .filter(ee.Filter.eq('nImages', 0))
    .aggregate_array('year');
print(emptyYears, 'Years with zero contributing scenes for this AOI')

// ---------------------------------------------------------------------
// NA-FILLED VERSION FOR EXPORT: keeps ALL years, 1984-2025, in the output
// time series -- years with zero contributing scenes stay in as fully
// masked (NA) bands rather than being dropped, so every exported band
// stack has the same year-to-band-index mapping regardless of AOI.
// -9999 is used as an explicit sentinel "no data" value: after unmask(),
// GEE writes -9999 into the GeoTIFF wherever a pixel was masked (i.e. every
// pixel in an empty year, or any pixel outside your AOI's cloud-free
// coverage). Read into R (e.g. with terra::rast()), set the NAflag to
// -9999 and those pixels come back in as proper NA, e.g.:
//   r <- terra::rast("Landsat_NDVI_1984_2025.tif")
//   terra::NAflag(r) <- -9999
// ---------------------------------------------------------------------
var allYearsImages = collectYear
    .map(convertLandsatBands)
    .map(function(image){
        return image.clip(roi_bbox).unmask(-9999);
    });
print(allYearsImages, 'All years (missing years filled as NA = -9999)')

// Test visualization
var testImage = nullimages.filter(ee.Filter.eq('year', 2023)).first()
print(testImage, 'test image 2023')

var ndvi2023 = testImage.select("NDVI")
var ndviParams = {min: 0, max: 10000, palette: ['blue', 'white', 'green']};
Map.addLayer(ndvi2023.clip(roi_bbox), ndviParams, 'NDVI 2023');

var ndmi2023 = testImage.select("NDMI")
var ndmiParams = {min: -10000, max: 10000, palette: ['red', 'white', 'blue']};
Map.addLayer(ndmi2023.clip(roi_bbox), ndmiParams, 'NDMI 2023');

var swir2023 = testImage.select("SWIR1_median")
var ndmiParams = {min: -10000, max: 10000, palette: ['red', 'white', 'blue']};
Map.addLayer(swir2023.clip(roi_bbox), ndmiParams, 'Swir 2023');

///////////////////////////////////////////////////////////////////////////////////////
// NEW: EXPORT FUNCTIONS BY VARIABLE WITH YEARLY BANDS
///////////////////////////////////////////////////////////////////////////////////////

// Define all variables/bands to export
var variablesToExport = [
  // Spectral bands
  'B_median', 'G_median', 'R_median', 'NIR_median', 'SWIR1_median', 'SWIR2_median',
  // Vegetation indices
  'EVI', 'NDVI', 'NIRV', 'GCC',
  // Water/moisture indices
  'NDMI', 'NDWI', 'MNDWI', 'NBR', 'LSWI',
  // QA bands
  'pixel_qa_median', 'radsat_qa_median'
];

/**

 * Function to create a multi-band image for a specific variable across all years
 * @param {string} variableName - Name of the variable/band to extract
 * @param {ee.ImageCollection} imageCollection - Collection of yearly composites
 * @returns {ee.Image} - Multi-band image with bands labeled by year
   */
   function createVariableTimeSeries(variableName, imageCollection) {
     // Sort collection by year to ensure proper order
     var sortedCollection = imageCollection.sort('year');

  // Function to rename band with year suffix
  var addYearSuffix = function(image) {
    var year = ee.Number(image.get('year')).format('%d');
    var bandName = ee.String(variableName).cat('_').cat(year);
    return image.select([variableName]).rename([bandName]);
  };

  // Map the function to create renamed bands
  var renamedCollection = sortedCollection.map(addYearSuffix);

  // Convert to multi-band image using toBands()
  var timeSeriesImage = renamedCollection.toBands();

  return timeSeriesImage;
}

/**

 * Function to export a variable as a multi-band time series
 * @param {string} variableName - Name of the variable to export
   */
   function exportVariableTimeSeries(variableName) {
     // Uses allYearsImages so every export has one band per year in yrange,
     // with missing years present as an explicit -9999 NA layer rather than
     // silently absent -- keeps the year->band mapping identical across AOIs.
     var timeSeries = createVariableTimeSeries(variableName, allYearsImages);

  Export.image.toDrive({
    image: timeSeries.clip(roi_bbox),
    description: 'TimeSeries_' + variableName,
    folder: 'GEE',
    fileNamePrefix: 'Landsat_' + variableName + '_' + yrange[0] + '_' + yrange[1],
    region: roi,
    scale: scale,
    crs: crs,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    // Embeds an actual NoData tag (-9999) in the GeoTIFF header -- matches
    // the value used in allYearsImages' unmask(-9999) above. Without this,
    // the exported -9999s are just ordinary pixel values with no metadata
    // flag, and R/QGIS/etc. have no way to know they mean "missing".
    formatOptions: {noData: -9999}
  });

  print('Export task created for variable: ' + variableName);
}

/**

 * Function to export all variables as separate multi-band time series
   */
   function exportAllVariableTimeSeries() {
     variablesToExport.forEach(function(variable) {
   exportVariableTimeSeries(variable);
     });

  print('Export tasks created for all variables. Check the Tasks tab to run them.');
}

/**

 * Function to export specific variables (modify as needed)
   */
   function exportSpecificVariables() {
     // Example: Export only key vegetation and moisture indices
     var keyVariables = ['NDVI', 'EVI', 'NDMI', 'NDWI', 'B_median', 'G_median', 'R_median', 'NIR_median'];

  keyVariables.forEach(function(variable) {
    exportVariableTimeSeries(variable);
  });

  print('Export tasks created for key variables: ' + keyVariables.join(', '));
}

// NOTE: switched from `nullimages` to `allYearsImages` so the exported time
// series always has one band per calendar year in [yrange], with missing
// years present as an explicit NA (-9999) layer instead of being skipped.
// If you'd rather go back to silently dropping empty years, swap this back
// to `nullimages`.
var ndviTimeSeries = createVariableTimeSeries('NDVI', allYearsImages);
print(ndviTimeSeries, 'NDVI Time Series');

// Visualize NDVI for a specific year from the time series
var ndvi2020FromSeries = ndviTimeSeries.select('NDVI_2020');
Map.addLayer(ndvi2020FromSeries.clip(roi_bbox), ndviParams, 'NDVI 2020 from Time Series');

///////////////////////////////////////////////////////////////////////////////////////
// EXECUTION: Choose which export function to run
///////////////////////////////////////////////////////////////////////////////////////

// Option 1: Export all variables as separate time series (CAUTION: Creates many tasks!)
// exportAllVariableTimeSeries();

// Option 2: Export only specific key variables
//exportSpecificVariables();

// Option 3: Export individual variables (uncomment as needed)
exportVariableTimeSeries('NDVI');
// exportVariableTimeSeries('EVI');
// exportVariableTimeSeries('NDMI');
// exportVariableTimeSeries('B_median');

print('Check the Tasks tab to run the export tasks.');

///////////////////////////////////////////////////////////////////////////////////////
// ALTERNATIVE: Export to Assets for further processing
///////////////////////////////////////////////////////////////////////////////////////

/**

 * Function to export variable time series to Earth Engine Assets
   */
   function exportVariableToAsset(variableName) {
     var timeSeries = createVariableTimeSeries(variableName, allYearsImages);

  Export.image.toAsset({
    image: timeSeries.clip(roi_bbox),
    description: 'Asset_TimeSeries_' + variableName,
    assetId: 'users/yourusername/TimeSeries_' + variableName, // Change to your username
    region: roi,
    scale: scale,
    crs: crs,
    maxPixels: 1e13
  });
}
