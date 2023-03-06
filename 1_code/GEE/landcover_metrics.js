/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var table = ee.FeatureCollection("users/bgcasey/harvest_severity_birds/ss_xy");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// Run this script using the Earth Engine code editor at code.earthengine.google.com

//########################################################################################################
//##### User defined inputs ##### 
//########################################################################################################

// import ibutton xy locations
var ss_xy = ee.FeatureCollection("users/bgcasey/harvest_severity_birds/ss_xy");

// define a buffer size around point locations (for zonal stats)
var buf=600 

// Date range
var Date_Start = ee.Date('1984-01-01');
var Date_End = ee.Date('2020-12-31');
  
 
 //set number of months in time interval
var interval=12;

// import shapefile of study area
var study_area = ee.FeatureCollection("projects/ee-bgcasey-climate/assets/Alberta_boundary");
// var study_area = ee.FeatureCollection("users/bgcasey/studyarea_CL");

//########################################################################################################
//##### Setup ##### 
//########################################################################################################

// for zonal stats create buffer around points
var ss_xy_buff= ss_xy.map(function(pt){
    return pt.buffer(buf);
  });
print("ss_xy_buff", ss_xy_buff.limit(10))
  
var aoi = study_area.geometry();

// convert the geometry to a feature to get the batch.Download.ImageCollection.toDrive function to work
var aoi1=ee.FeatureCollection(aoi)

// Create list of dates for time series. It start at the first of each month in the date range and progress by num_months_in_interval
var n_months = Date_End.difference(Date_Start,'month').round();
var dates = ee.List.sequence(0, n_months, interval);
var make_datelist = function(n) {
  return Date_Start.advance(n,'month');
};
dates = dates.map(make_datelist);

print('list of dates for time series', dates)

//########################################################################################################
//##### Get spatial variables
//########################################################################################################


////////////////////////////////////////
// Canopy height
////////////////////////////////////////

var canopy_height = ee.Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1')
      .rename("canopy_height")
      .clip(aoi)

var canopy_standard_deviation = ee.Image('users/nlang/ETH_GlobalCanopyHeightSD_2020_10m_v1')
      .rename('canopy_standard_deviation')
      .clip(aoi)


//combine bands into single image
var canopy = canopy_height.addBands([canopy_standard_deviation])
// print("canopy", canopy)

var ev_fixed = canopy.reduceRegions({
  collection: ss_xy_buff,
  reducer: ee.Reducer.mean(),
  scale: 30
});

print("ev_fixed", ev_fixed.limit(10))
////////////////////////////////////////
// Landcover Indices
////////////////////////////////////////

var LC = ee.ImageCollection('projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2').
filterDate(Date_Start, Date_End);
print(LC)          


// choose reducers
var reducers = ee.Reducer.count().combine({
  reducer2: ee.Reducer.frequencyHistogram(),
  sharedInputs: true
});

var LC_1 = LC.map(function(img) {
  return img.reduceRegions({
    collection: ev_fixed,
    reducer: reducers, // set the names of output properties to the corresponding band names
    scale: 30,
    tileScale: 2
  }).map(function (feature) {
            var histogramResults = ee.Dictionary(feature.get('histogram'));
            var pixel_count= ee.Number(feature.get('count'))
      return feature.copyProperties(img, ['system:time_start']) //to get year properties from the stack
            .set(// get proportion of landcover from histogram 
                 // by dividing histogram pixel values by the total pixel_count.
                'Unclassified', ee.Number(histogramResults.get('0', 0)).divide(pixel_count),
                'Water', ee.Number(histogramResults.get('20', 0)).divide(pixel_count),
                'Snow_Ice', ee.Number(histogramResults.get('31', 0)).divide(pixel_count),
                'Rock_Rubble', ee.Number(histogramResults.get('32', 0)).divide(pixel_count),
                'Exposed_Barren_land', ee.Number(histogramResults.get('33', 0)).divide(pixel_count),
                'Bryoids', ee.Number(histogramResults.get('40', 0)).divide(pixel_count),
                'Shrubs', ee.Number(histogramResults.get('50', 0)).divide(pixel_count),
                'Wetland', ee.Number(histogramResults.get('80', 0)).divide(pixel_count),
                'Wetland-treed', ee.Number(histogramResults.get('81', 0)).divide(pixel_count),
                'Herbs', ee.Number(histogramResults.get('100', 0)).divide(pixel_count),
                'Coniferous', ee.Number(histogramResults.get('210', 0)).divide(pixel_count),
                'Broadleaf', ee.Number(histogramResults.get('220', 0)).divide(pixel_count),
                'Mixedwood', ee.Number(histogramResults.get('230', 0)).divide(pixel_count),
                'landcover_yr', img.date().format('YYYY'));
  })
}).flatten(); //  Flattens collections of collections into a feature collection of those collections

print("LC1", LC_1.limit(2))



// // // //########################################################################################################
// // // // // ### Save/export landcover data ###
// // // //########################################################################################################

// Export landcover data to a csv
Export.table.toDrive({
  folder: 'harvest_severity',
  collection: LC_1,
  description:'gee_metrics',
  fileFormat: 'csv'
    // selectors: [ // choose properties to include in export table
    //               'SS', 
    //               'srvy_yr',
    //               'landcover_yr',
    //               'count'
    //               ] 
});

//########################################################################################################
// // ### Map ###
//########################################################################################################

// Map.centerObject(aoi, 6) // center the map on the study area
 
// // add ibutton locations
// Map.addLayer(ss_xy,{color: 'bf1b29'}, "ss_xy")