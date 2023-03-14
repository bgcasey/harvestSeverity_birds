/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var table = ee.FeatureCollection("projects/ee-bgcasey-harvest-birds/assets/o18_HarvestAreas_HFI_2019_buff_neg30_surveyed_simp_2");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//###########################################################################################################
//## This script was adapted from the following:
//##
//##                                                                                                                
//## EXTRACTING & EXPORTING POST-HARVEST REGENERATION METRICS FOR FOREST HARVEST AREAS IN ALBERTA, USING LANDTRENDR 
//##                                                                                                               
//##         AUTHORS: Jennifer Hird (University of Calgary | Alberta Biodiversity Monitoring Institute)
//##               AND the authors of LandTrendr (see below)
//##         CONTACT: jennifer.hird@ucalgary.ca
//##                                                                                                                               #\\
//## LANDTRENDR GREATEST MAGNITUDE DISTURBANCE MAPPING                                            
//##
//##         date: 2018-04-19
//##         author: Zhiqiang Yang  | zhiqiang.yang@oregonstate.edu
//##                 Justin Braaten | jstnbraaten@gmail.com
//##                 Robert Kennedy | rkennedy@coas.oregonstate.edu
//##         website: https://github.com/eMapR/LT-GEE
///########################################################################################################


//########################################################################################################
//##### INPUTS ##### 
//########################################################################################################

///////////////////////
///// Point count locations
///////////////////////
var ss_xy= ee.FeatureCollection("projects/ee-bgcasey-harvest-birds/assets/ss_xy");

///////////////////////
/////Harvest polygons
///////////////////////

// // define interior buffer
// var buf=-30

var fullHA = ee.FeatureCollection("projects/ee-bgcasey-harvest-birds/assets/o18_HarvestAreas_HFI_2019_buff_neg30_surveyed_simp_2");
      // .filter(ee.Filter.bounds(ss_xy));

print('fullHA_size', fullHA.size());

// var HA_s = ee.FeatureCollection("projects/ee-bgcasey-harvest-birds/assets/HA_surveyed");
// print('HA_s', HA_s.limit(10));

// Export.table.toAsset(HA_surveyed, "HA_surveyed_export")

//// create an interior buffer in harvest polygons
// var fullHA_buf= fullHA.map(function(pt){
//     return pt.buffer(buf);
//   });

// // select a subset of polygons to try code  
// var subset = 5

// var filterSpecHAs = ee.Filter.lte('uniquID', subset);

// // //== OR ::
var filterSpecHAs = ee.Filter.and(ee.Filter.gt('uniquID',4800),ee.Filter.lte('uniquID',5000)); 


//== Apply filter to select portion of polygons 
var HA = fullHA.filter(filterSpecHAs);
print('HA_size', HA.size());


// print("HA", HA.size())
// // For final analysis use:
// var HA = fullHA;


// ///////////////////////
// /////Study area
// ///////////////////////

// Create geometry object, covering Alberta (for filtering Landsat data)
var geometry = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[-120.7177734375, 53.90433815627471],
          [-114.43359375, 48.574789910928864],
          [-109.3359375, 48.45835188280866],
          [-109.2041015625, 60.19615576604439],
          [-120.8935546875, 60.174306261926034]]]);

// define a geometry - there are lots of ways to do this, see the GEE User guide
var aoi = geometry; // should be a GEE geometry object - here we are getting it from an drawn polygon


///////////////////////
/////Parameters
///////////////////////

// visualization parameters
var imageVisParam = {"opacity":1,"bands":["yr1985","yr2001","yr2017"],
"min":473.9825265803124,"max":915.2623445519782,"gamma":1};

// define years and dates to include in landsat image collection
var startYear  = 1984;    // what year do you want to start the time series 
// var endYear    = 2017;    // what year do you want to end the time series ***** Set to 2017, year of HFI 2017 HA polygons
//** New for 2018 HAs (2020-06-17)
var endYear    = 2021;    // what year do you want to end the time series ***** Set to survey year
var startDay   = '06-01'; // what is the beginning of date filter | month-day
var endDay     = '09-30'; // what is the end of date filter | month-day


// define disturbance mapping filter parameters 
var treeLoss1  = 175;      // delta filter for 1 year duration disturbance, <= will not be included as disturbance - units are in units of segIndex defined in the following function definition
var treeLoss20 = 200;      // delta filter for 20 year duration disturbance, <= will not be included as disturbance - units are in units of segIndex defined in the following function definition
var preVal     = 400;      // pre-disturbance value threshold - values below the provided threshold will exclude disturbance for those pixels - units are in units of segIndex defined in the following function definition
var mmu        = 15;       // minimum mapping unit for disturbance patches - units of pixels


// define function to calculate a spectral index to segment with LT
var segIndex = function(img) {
    var index = img.normalizedDifference(['B4', 'B7'])                      // calculate normalized difference of band 4 and band 7 (B4-B7)/(B4+B7)
                  .multiply(1000)                                          // ...scale results by 1000 so we can convert to int and retain some precision
                  .select([0], ['NBR'])                                    // ...name the band
                  .set('system:time_start', img.get('system:time_start')); // ...set the output system:time_start metadata to the input image time_start otherwise it is null
    return index ;
};

var distDir = -1; // define the sign of spectral delta for vegetation loss for the segmentation index - 
                  // NBR delta is negetive for vegetation loss, so -1 for NBR, 1 for band 5, -1 for NDVI, etc


// define the segmentation parameters:
// reference: Kennedy, R. E., Yang, Z., & Cohen, W. B. (2010). Detecting trends in forest disturbance and recovery using yearly Landsat time series: 1. LandTrendr—Temporal segmentation algorithms. Remote Sensing of Environment, 114(12), 2897-2910.
//            https://github.com/eMapR/LT-GEE
var run_params = { 
  maxSegments:            6,
  spikeThreshold:         0.9,
  vertexCountOvershoot:   3,
  preventOneYearRecovery: true,
  recoveryThreshold:      0.5, //0.25, *** CHANGED BY JH (2018-12-07)
  pvalThreshold:          0.05,
  bestModelProportion:    0.75,
  minObservationsNeeded:  6
};



//########################################################################################################
//##### ANNUAL SR TIME SERIES COLLECTION BUILDING FUNCTIONS ##### 
//########################################################################################################

//----- MAKE A DUMMY COLLECTOIN FOR FILLTING MISSING YEARS -----
var dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))]); // make an image collection from an image with 6 bands all set to 0 and then make them masked values


//------ L8 to L7 HARMONIZATION FUNCTION -----
// slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, Remote Sensing of Environment, 185, 57-70.(http://dx.doi.org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients
var harmonizationRoy = function(oli) {
  var slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949]);        // create an image of slopes per band for L8 TO L7 regression line - David Roy
  var itcp = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029]);     // create an image of y-intercepts per band for L8 TO L7 regression line - David Roy
  var y = oli.select(['B2','B3','B4','B5','B6','B7'],['B1', 'B2', 'B3', 'B4', 'B5', 'B7']) // select OLI bands 2-7 and rename them to match L7 band names
            .resample('bicubic')                                                          // ...resample the L8 bands using bicubic
            .subtract(itcp.multiply(10000)).divide(slopes)                                // ...multiply the y-intercept bands by 10000 to match the scale of the L7 bands then apply the line equation - subtract the intercept and divide by the slope
            .set('system:time_start', oli.get('system:time_start'));                      // ...set the output system:time_start metadata to the input image time_start otherwise it is null
  return y.toShort();                                                                       // return the image as short to match the type of the other data
};


//------ RETRIEVE A SENSOR SR COLLECTION FUNCTION -----
var getSRcollection = function(year, startDay, endDay, sensor, aoi) {
  // ** ADDED IN BY JEN HIRD:::: create startDate & endDate ee.Date objects for use in functions
  // because get error about bad date/time otherwise
  // ---------------------------------------------
  var mo1 = ee.Number.parse(startDay.substr(0,2));
  var day1 = ee.Number.parse(startDay.substr(3,2));
  var mo2 = ee.Number.parse(endDay.substr(0,2));
  var day2 = ee.Number.parse(endDay.substr(3,2));
  var startDate = ee.Date.fromYMD(year, mo1, day1);
  var endDate = ee.Date.fromYMD(year, mo2, day2);
  // ---------------------------------------------
  
  // get a landsat collection for given year, day range, and sensor
  var srCollection = ee.ImageCollection('LANDSAT/'+ sensor + '/C01/T1_SR') // get surface reflectance images
                      .filterBounds(aoi)                                  // ...filter them by intersection with AOI
                      .filterDate(startDate, endDate); //*** ADDED BY JH (to deal with bad date error) ***
                      //.filterDate(year+'-'+startDay, year+'-'+endDay);    // ...filter them by year and day range
  
  // apply the harmonization function to LC08 (if LC08), subset bands, unmask, and resample           
  srCollection = srCollection.map(function(img) {
    var dat = ee.Image(
      ee.Algorithms.If(
        sensor == 'LC08',                                                  // condition - if image is OLI
        harmonizationRoy(img.unmask()),                                    // true - then apply the L8 TO L7 alignment function after unmasking pixels that were previosuly masked (why/when are pixels masked)
        img.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])                   // false - else select out the reflectance bands from the non-OLI image
          .unmask()                                                       // ...unmask any previously masked pixels 
          .resample('bicubic')                                            // ...resample by bicubic 
          .set('system:time_start', img.get('system:time_start'))         // ...set the output system:time_start metadata to the input image time_start otherwise it is null
      )
    );
    
    // make a cloud, cloud shadow, and snow mask from fmask band
    var qa = img.select('pixel_qa');                                       // select out the fmask band
    var mask = qa.bitwiseAnd(8).eq(0).and(                                 // include shadow
              qa.bitwiseAnd(16).eq(0)).and(                               // include snow
              qa.bitwiseAnd(32).eq(0));                                   // include clouds
    
    // apply the mask to the image and return it
    return dat.mask(mask); //apply the mask - 0's in mask will be excluded from computation and set to opacity=0 in display
  });

  return srCollection; // return the prepared collection
};


//------ FUNCTION TO COMBINE LT05, LE07, & LC08 COLLECTIONS -----
var getCombinedSRcollection = function(year, startDay, endDay, aoi) {
    var lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi);       // get TM collection for a given year, date range, and area
    var le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi);       // get ETM+ collection for a given year, date range, and area
    var lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi);       // get OLI collection for a given year, date range, and area
    var mergedCollection = ee.ImageCollection(lt5.merge(le7).merge(lc8)); // merge the individual sensor collections into one imageCollection object
    return mergedCollection;                                              // return the Imagecollection
};


//------ FUNCTION TO REDUCE COLLECTION TO SINGLE IMAGE PER YEAR BY MEDOID -----
/*
  LT expects only a single image per year in a time series, there are lost of ways to
  do best available pixel compositing - we have found that a mediod composite requires little logic
  is robust, and fast
  
  Medoids are representative objects of a data set or a cluster with a data set whose average 
  dissimilarity to all the objects in the cluster is minimal. Medoids are similar in concept to 
  means or centroids, but medoids are always members of the data set.
*/

// make a medoid composite with equal weight among indices
var medoidMosaic = function(inCollection, dummyCollection) {
  
  // fill in missing years with the dummy collection
  var imageCount = inCollection.toList(1).length();                                                            // get the number of images 
  var finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), inCollection, dummyCollection)); // if the number of images in this year is 0, then use the dummy collection, otherwise use the SR collection
  
  // calculate median across images in collection per band
  var median = finalCollection.median();                                                                       // calculate the median of the annual image collection - returns a single 6 band image - the collection median per band
  
  // calculate the different between the median and the observation per image per band
  var difFromMedian = finalCollection.map(function(img) {
    var diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));                                       // get the difference between each image/band and the corresponding band median and take to power of 2 to make negatives positive and make greater differences weight more
    return diff.reduce('sum').addBands(img);                                                                   // per image in collection, sum the powered difference across the bands - set this as the first band add the SR bands to it - now a 7 band image collection
  });
  
  // get the medoid by selecting the image pixel with the smallest difference between median and observation per band 
  return ee.ImageCollection(difFromMedian).reduce(ee.Reducer.min(7)).select([1,2,3,4,5,6], ['B1','B2','B3','B4','B5','B7']); // find the powered difference that is the least - what image object is the closest to the median of teh collection - and then subset the SR bands and name them - leave behind the powered difference band
};


//------ FUNCTION TO APPLY MEDOID COMPOSITING FUNCTION TO A COLLECTION -------------------------------------------
var buildMosaic = function(year, startDay, endDay, aoi, dummyCollection) {                                                                      // create a temp variable to hold the upcoming annual mosiac
  var collection = getCombinedSRcollection(year, startDay, endDay, aoi);  // get the SR collection
  var img = medoidMosaic(collection, dummyCollection)                     // apply the medoidMosaic function to reduce the collection to single image per year by medoid 
              .set('system:time_start', (new Date(year,8,1)).valueOf());  // add the year to each medoid image - the data is hard-coded Aug 1st 
  return ee.Image(img);                                                   // return as image object
};


//------ FUNCTION TO BUILD ANNUAL MOSAIC COLLECTION ------------------------------
var buildMosaicCollection = function(startYear, endYear, startDay, endDay, aoi, dummyCollection) {
  var imgs = [];                                                                    // create empty array to fill
  for (var i = startYear; i <= endYear; i++) {                                      // for each year from hard defined start to end build medoid composite and then add to empty img array
    var tmp = buildMosaic(i, startDay, endDay, aoi, dummyCollection);               // build the medoid mosaic for a given year
    imgs = imgs.concat(tmp.set('system:time_start', (new Date(i,8,1)).valueOf()));  // concatenate the annual image medoid to the collection (img) and set the date of the image - hard coded to the year that is being worked on for Aug 1st
  }
  return ee.ImageCollection(imgs);                                                  // return the array img array as an image collection
};



//********************************************************************************************************
//****************************** *JH FUNCTIONS*: BUILD MEDIAN COMPOSITE TS *******************************
//********************************************************************************************************

//== Function for computing median per summer season from an image collection
var calcSeasonalMedian = function(imgCollection, yr){
  var sysTime = ee.Image(imgCollection.first()).get('system:time_start');
  var seasMed = imgCollection.median();
  var year = ee.Image(ee.Number(yr)).rename("YEAR");
  return seasMed.set('year', yr).set('system:time_start',sysTime).addBands(year).toFloat();
};

//== Create sequential list of years 
var ltYrs = ee.List.sequence(startYear, endYear);
// print("ltYrs", ltYrs)

//== For each year in sequence, create a combined Landsat surface reflectance image collection,
//== and apply seasonal median function to create a single composite
var buildAnnualMedian = ltYrs.map(function(yr){
  var comboCollec = getCombinedSRcollection(yr, startDay, endDay, aoi);
  var medComp = calcSeasonalMedian(comboCollec, yr);
  return medComp;
});

//== Define function to calculate a spectral index to segment with LT (SOURCE: LandTrendr segIndex function)
var NBRindex = function(img) {
    var index = img.normalizedDifference(['B4', 'B7'])                      // calculate normalized difference of band 4 and band 7 (B4-B7)/(B4+B7)
                  .multiply(1000).multiply(distDir)                                          // ...scale results by 1000 so we can convert to int and retain some precision
                  .select([0], ['NBR'])                                    // ...name the band
                  .set('system:time_start', img.get('system:time_start')) // ...set the output system:time_start metadata to the input image time_start otherwise it is null
                  .set('year', img.get('year'));
    return index ;
};

//********************************************************************************************************
//********************************************************************************************************


//########################################################################################################
//##### UNPACKING LT-GEE OUTPUT STRUCTURE FUNCTIONS ##### 
//########################################################################################################

// ----- FUNCTION TO EXTRACT VERTICES FROM LT RESULTS AND STACK BANDS -----
var getLTvertStack = function(LTresult) {
  var emptyArray = [];                              // make empty array to hold another array whose length will vary depending on maxSegments parameter    
  var vertLabels = [];                              // make empty array to hold band names whose length will vary depending on maxSegments parameter 
  var iString;                                      // initialize variable to hold vertex number
  for(var i=1;i<=run_params.maxSegments+1;i++){     // loop through the maximum number of vertices in segmentation and fill empty arrays
    iString = i.toString();                         // define vertex number as string 
    vertLabels.push("vert_"+iString);               // make a band name for given vertex
    emptyArray.push(0);                             // fill in emptyArray
  }
  
  var zeros = ee.Image(ee.Array([emptyArray,        // make an image to fill holes in result 'LandTrendr' array where vertices found is not equal to maxSegments parameter plus 1
                                emptyArray,
                                emptyArray]));
  
  var lbls = [['yrs_','src_','fit_'], vertLabels,]; // labels for 2 dimensions of the array that will be cast to each other in the final step of creating the vertice output 

  var vmask = LTresult.arraySlice(0,3,4);           // slices out the 4th row of a 4 row x N col (N = number of years in annual stack) matrix, which identifies vertices - contains only 0s and 1s, where 1 is a vertex (referring to spectral-temporal segmentation) year and 0 is not
  
  var ltVertStack = LTresult.arrayMask(vmask)       // uses the sliced out isVert row as a mask to only include vertice in this data - after this a pixel will only contain as many "bands" are there are vertices for that pixel - min of 2 to max of 7. 
                      .arraySlice(0, 0, 3)          // ...from the vertOnly data subset slice out the vert year row, raw spectral row, and fitted spectral row
                      .addBands(zeros)              // ...adds the 3 row x 7 col 'zeros' matrix as a band to the vertOnly array - this is an intermediate step to the goal of filling in the vertOnly data so that there are 7 vertice slots represented in the data - right now there is a mix of lengths from 2 to 7
                      .toArray(1)                   // ...concatenates the 3 row x 7 col 'zeros' matrix band to the vertOnly data so that there are at least 7 vertice slots represented - in most cases there are now > 7 slots filled but those will be truncated in the next step
                      .arraySlice(1, 0, run_params.maxSegments+1) // ...before this line runs the array has 3 rows and between 9 and 14 cols depending on how many vertices were found during segmentation for a given pixel. this step truncates the cols at 7 (the max verts allowed) so we are left with a 3 row X 7 col array
                      .arrayFlatten(lbls, '');      // ...this takes the 2-d array and makes it 1-d by stacking the unique sets of rows and cols into bands. there will be 7 bands (vertices) for vertYear, followed by 7 bands (vertices) for rawVert, followed by 7 bands (vertices) for fittedVert, according to the 'lbls' list

  return ltVertStack;                               // return the stack
};





//########################################################################################################
//##### GREATEST DISTURBANCE EXTRACTION FUNCTIONS #####
//########################################################################################################

// ----- function to extract greatest disturbance based on spectral delta between vertices 
var extractDisturbance = function(lt, distDir, params, mmu) {
  // select only the vertices that represents a change
  var vertexMask = lt.arraySlice(0, 3, 4); // get the vertex - yes(1)/no(0) dimension
  var vertices = lt.arrayMask(vertexMask); // convert the 0's to masked
  
  // construct segment start and end point years and index values
  var left = vertices.arraySlice(1, 0, -1);    // slice out the vertices as the start of segments
  var right = vertices.arraySlice(1, 1, null); // slice out the vertices as the end of segments
  var startYear = left.arraySlice(0, 0, 1);    // get year dimension of LT data from the segment start vertices
  var startVal = left.arraySlice(0, 2, 3);     // get spectral index dimension of LT data from the segment start vertices
  var endYear = right.arraySlice(0, 0, 1);     // get year dimension of LT data from the segment end vertices 
  var endVal = right.arraySlice(0, 2, 3);      // get spectral index dimension of LT data from the segment end vertices
  
  var dur = endYear.subtract(startYear);       // subtract the segment start year from the segment end year to calculate the duration of segments 
  var mag = endVal.subtract(startVal);         // substract the segment start index value from the segment end index value to calculate the delta of segments 
  
  // concatenate segment start year, delta, duration, and starting spectral index value to an array 
  var distImg = ee.Image.cat([startYear.add(1), mag, dur, startVal.multiply(distDir)]).toArray(0); // make an image of segment attributes - multiply by the distDir parameter to re-orient the spectral index if it was flipped for segmentation - do it here so that the subtraction to calculate segment delta in the above line is consistent - add 1 to the detection year, because the vertex year is not the first year that change is detected, it is the following year
  
  // sort the segments in the disturbance attribute image delta by spectral index change delta  
  var distImgSorted = distImg.arraySort(mag.multiply(-1));                                  // flip the delta around so that the greatest delta segment is first in order

  // slice out the first (greatest) delta
  var tempDistImg = distImgSorted.arraySlice(1, 0, 1).unmask(ee.Image(ee.Array([[0],[0],[0],[0]])));  // get the first segment in the sorted array
 
  // ### ADDED BY JH (2019-01-18): slice out second-greatest delta from array
  var tempSecMag = distImgSorted.arraySlice(1,1,2).unmask(ee.Image(ee.Array([[0],[0],[0],[0]]))); // get second segment in the sorted array
  
  // make an image from the array of attributes for the greatest disturbance
  var finalDistImg = ee.Image.cat(tempDistImg.arraySlice(0,0,1).arrayProject([1]).arrayFlatten([['yod']]),     // slice out year of disturbance detection and re-arrange to an image band 
                                  tempDistImg.arraySlice(0,1,2).arrayProject([1]).arrayFlatten([['mag']]),     // slice out the disturbance magnitude and re-arrange to an image band 
                                  tempDistImg.arraySlice(0,2,3).arrayProject([1]).arrayFlatten([['dur']]),     // slice out the disturbance duration and re-arrange to an image band
                                  tempDistImg.arraySlice(0,3,4).arrayProject([1]).arrayFlatten([['preval']]),  // slice out the pre-disturbance spectral value and re-arrange to an image band
                                  tempSecMag.arraySlice(0,0,1).arrayProject([1]).arrayFlatten([['yod2nd']]),   // ### ADDED BY JH (2019-01-18): slice out year of disturbance for second-greatest change
                                  tempSecMag.arraySlice(0,1,2).arrayProject([1]).arrayFlatten([['mag2nd']]));  // ### ADDED BY JH (2019-01-18): slice out magnitude for second-greatest change
 

  // filter out disturbances based on user settings
  var threshold = ee.Image(finalDistImg.select(['dur']))                        // get the disturbance band out to apply duration dynamic disturbance magnitude threshold 
                    .multiply((params.tree_loss20 - params.tree_loss1) / 19.0)  // ...
                    .add(params.tree_loss1)                                     //    ...interpolate the magnitude threshold over years between a 1-year mag thresh and a 20-year mag thresh
                    .lte(finalDistImg.select(['mag']))                          // ...is disturbance less then equal to the interpolated, duration dynamic disturbance magnitude threshold 
                    .and(finalDistImg.select(['mag']).gt(0))                    // and is greater than 0  
                    .and(finalDistImg.select(['preval']).gt(params.pre_val));   // and is greater than pre-disturbance spectral index value threshold
  
  // apply the filter mask
  finalDistImg = finalDistImg/*.mask(threshold)*/.int16(); // ### Commented out by JH; don't want masking at this point in process ###
  
  // ### Commented out by JH - don't want masking; want to see all pixels ###
  // // patchify the remaining disturbance pixels using a minimum mapping unit
  // if(mmu > 1){
  //   var mmuPatches = finalDistImg.select(['yod'])           // patchify based on disturbances having the same year of detection
  //                           .connectedPixelCount(mmu, true) // count the number of pixel in a candidate patch
  //                           .gte(mmu);                      // are the the number of pixels per candidate patch greater than user-defined minimum mapping unit?
  //   finalDistImg = finalDistImg.updateMask(mmuPatches);     // mask the pixels/patches that are less than minimum mapping unit
  // } 
  
  return finalDistImg; // return the filtered greatest disturbance attribute image
};

//########################################################################################################
//##### BUILD COLLECTION AND RUN LANDTRENDR #####
//########################################################################################################

//----- BUILD LT COLLECTION -----
// build annual surface reflection collection
var annualSRcollection = buildMosaicCollection(startYear, endYear, startDay, endDay, aoi, dummyCollection); // put together the cloud-free medoid surface reflectance annual time series collection

// // apply the function to calculate the segmentation index and adjust the values by the distDir parameter - flip index so that a vegetation loss is associated with a postive delta in spectral value
// var annualIndexCollection = annualSRcollection.map(segIndex)                                             // map the function over every image in the collection - returns a 1-band annual image collection of the spectral index
//                                               .map(function(img) {return img.multiply(distDir)           // ...multiply the segmentation index by the distDir to ensure that vegetation loss is associated with a positive spectral delta
//                                               .set('system:time_start', img.get('system:time_start'))}); // ...set the output system:time_start metadata to the input image time_start otherwise it is null

// // add the segmentation index as a second band to the annual collection so that the LT fuction will return a "fitted" copy of the annual time series for this index   
// var indexNameFTV = ee.Image(annualIndexCollection.first()).bandNames().getInfo()[0]+'_FTV'; // create a name for the to-be-fitted band - get the name of the segmentation index and append "_FTV"
// var ltCollection = annualIndexCollection.map(function(img) {                                // start anonymous function to add the band
//   return img.addBands(img.select([0],[indexNameFTV])                                        // duplicate the segmentation index as a second band using the name that was just created
//                         .multiply(distDir))                                                // ...flip the values around so that it is back to its original orientation
//                         .set('system:time_start', img.get('system:time_start'));           // ...set the output system:time_start metadata to the input image time_start otherwise it is null                                 
// });


//********************************************************************************************************
//************************ JH: Use funcitons to create median composite time series **********************
//********************************************************************************************************

//== Build annual seasonal median composite
var medianSeries = ee.ImageCollection(buildAnnualMedian);
// print(medianSeries, "medianSeries")

//== Apply NBRindex function to calculate an NBR composite time series
var NBRMedseries = medianSeries.map(NBRindex);

//== add the segmentation index as a second band to the annual collection so that the LT fuction will return a "fitted" copy of the annual time series for this index   
var indexNameFTV = ee.Image(NBRMedseries.first()).bandNames().getInfo()[0]+'_FTV'; // create a name for the to-be-fitted band - get the name of the segmentation index and append "_FTV"
var NBRMedseriesFTV = NBRMedseries.map(function(img) {                                // start anonymous function to add the band
  return img.addBands(img.select([0],[indexNameFTV]))                                        // duplicate the segmentation index as a second band using the name that was just created
                        //.multiply(distDir))                                                // ...flip the values around so that it is back to its original orientation
                        .set('system:time_start', img.get('system:time_start'));           // ...set the output system:time_start metadata to the input image time_start otherwise it is null                                 
});

// //== PRINT OUTPUTS
// print('medianSeries:', medianSeries);
// print('NBRMedseries', NBRMedseries);
// print('NBRMedseriesFTV', NBRMedseriesFTV);
// print('ltCollection:', ltCollection);

// Create multiband image from NBR series image collection
var NBRImgResult = ee.Image(0);
for (var i = startYear; i <= endYear /*<= NBRseries.size()*/; i++){
  var img = ee.Image(NBRMedseries.filterMetadata('year','equals',i).first());
  var bandName = "yr"+i;
  var newBand = img.select([0]).rename(bandName).toFloat().multiply(distDir);
  var NBRImgResult = NBRImgResult.addBands(newBand);
}
NBRImgResult = NBRImgResult.select(NBRImgResult.bandNames().removeAll(["constant"]));


//********************************************************************************************************
//********************************************************************************************************


//----- RUN LANDTRENDR -----
// run_params.timeSeries = ltCollection;               // add LT collection to the segmentation run parameter object
run_params.timeSeries = NBRMedseriesFTV; // *** ADDED BY JH
var lt = ee.Algorithms.TemporalSegmentation.LandTrendr(run_params); // run LandTrendr spectral temporal segmentation algorithm

//== ### ADDED BY JH ###
// //== PRINT OUTPUTS
// print('lt',lt);

//== ### ADDED BY JH ###
//== Extract LandTrendr band
var ltResult = lt.select(["LandTrendr"]);
// //== PRINT OUTPUTS
// print('ltResult', ltResult);

// extract the segmentation-fitted index stack as 
var years = [];                                                           // make an empty array to hold year band names
// for (var i = startYear; i <= endYear; ++i) years.push('yr'+i.toString()); // fill the array with years from the startYear to the endYear and convert them to string
for (var i = startYear; i <= endYear; ++i) years.push(i.toString());

//== ### ADDED BY JH ###
// //== PRINT OUTPUTS
// print('years', years);

// Create multiband image of fitted index values
var ltFitStack = lt.select([1])                                           // select out the 2nd band data which is the segmentation-fitted spectral index 
                  .arrayFlatten([years]); 

//== ### ADDED BY JH ###
// //== PRINT OUTPUTS
// print('ltFitStack', ltFitStack);

//== ### ADDED BY JH ###
//== Use distDir to switch NBR back to positives (* -1)
var ltNBRFitted = ltFitStack.multiply(ee.Image.constant(distDir));
// print('ltNBRFitted',ltNBRFitted);


//########################################################################################################
//##### RUN THE GREATEST DISTURBANCE EXTRACT FUCTION #####
//########################################################################################################

// assemble the disturbance extraction parameters
var distParams = {
  tree_loss1: treeLoss1,
  tree_loss20: treeLoss20,  
  pre_val: preVal           
};

// run the dist extract function
var distImg = extractDisturbance(lt.select('LandTrendr'), distDir, distParams);

// //== ### ADDED BY JH: PRINT OUTPUTS
// print('distImg', distImg);


//########################################################################################################
//##################### CALCULATE POST-HARVEST/DISTURBANCE METRICS AND INFORMATION #######################
//########################################################################################################

//== Transform multi-band ltFitStack image into image collection (and convert back to positive values, for NBR)
var fitCollection = ee.ImageCollection(ltFitStack.bandNames().map(function(band) {
  var img = ltFitStack.select([band]).multiply(distDir).rename('NBR_fitted');
  var yr = ee.Number.parse(band);
  return img.set('year', yr);
}));

// // //== PRINT OUTPUTS
// print('fitCollection', fitCollection);


//-----------------------------------------
// ------- Pre-Disturbance Mean NBR -------    * WORKS! *

//== Create image from 'yod' band of distImg
// var yodImage = distImg.select(["yod"]);
var yodImage = distImg.select(["yod"]);

//== Function to mask out fitted NBR values that are before disturbance
var preDistColl = fitCollection.map(function(image) {
  var currYr = ee.Number.parse(image.get('year'));
  var yrImg = ee.Image.constant(currYr);
  var preDistDates = yrImg.lt(yodImage);
  var yearBand = ee.Image.constant(currYr).rename('Year').toFloat();
  return image.addBands(yearBand).updateMask(preDistDates);
});

//========= Using 5 years before harvest =========

// //== Create image from 'yod' band of distImg
// var yodImage = distImg.select(["yod"]);

//== Function to mask out fitted NBR values that are before disturbance
var preDist5yr = fitCollection.map(function(image) {
  var currYr = ee.Number.parse(image.get('year'));
  var yrImg = ee.Image.constant(currYr);
  var preDistDates = yrImg.lt(yodImage).and(yrImg.gte(yodImage.subtract(5)));
  var yearBand = ee.Image.constant(currYr).rename('Year').toFloat();
  return image.addBands(yearBand).updateMask(preDistDates);
});

//== Use reducer (mean) to calculate pre-disturbance mean NBR, per pixel, for 5 years before harvest
var preDistMean = preDist5yr.reduce(ee.Reducer.mean()).select(['NBR_fitted_mean']);

// //== PRINT OUTPUTS
// print('preDistMean', preDistMean);


//-------------------------------------------------
// ------- Post-Disturbance (YOD & EOD) NBR -------    * WORKS! *

//== Function to mask out fitted NBR values that are before the end of disturbance
var postDistColl = fitCollection.map(function(image) {
  var currYr = ee.Number.parse(image.get('year'));
  var yrImg = ee.Image.constant(currYr);
  var postDistDates = yrImg.gte(yodImage);
  var yearBand = ee.Image.constant(currYr).rename('Year').toFloat();
  return image.addBands(yearBand).updateMask(postDistDates);
});

// //== PRINT OUTPUTS
// print('postDistColl', postDistColl);

//== Create list of band names from ltFitStack multiband image
var listBands = ltFitStack.bandNames();

//== Map function over each band (which represents fitted NBR for that year),
//==  that masks out pixels for which that year that do not equal the value at 
//==  the biggest change (i.e., yod, calculated by LandTrendr functions)
var yodMaskedStack = listBands.map(function(band){
    var bandYrImg = ee.Image.constant(ee.Number.parse(band));
    var yodMask = bandYrImg.eq(yodImage);
    var relvVals = ltFitStack.select([band]).multiply(distDir).toFloat()
    .mask(yodMask).set('year', ee.Number.parse(band)).rename('yodVal');
    return relvVals;
  });

//== Convert outputted list to image collection
var yodValColl = ee.ImageCollection(yodMaskedStack);

// //== PRINT OUTPUTS
// print('yodMaskedStack', yodMaskedStack);
// print('yodValColl',yodValColl);

//== Combine image collection into one image (each pixel should be masked out in
//==  every image in the collection save one - the date of YOD)
var yodValImg = yodValColl.reduce(ee.Reducer.median()).rename('endDistVal').toFloat();

// //== PRINT OUTPUTS
// print('yodValImg',yodValImg);

//== Create empty start image filled with 0s, with correct (empty) bands and (empty) year property
var startImg = ee.List([ee.Image(0).select([0],['NBR_fitted'])
                                    .addBands(ee.Image(0).rename('Year'))
                                    .addBands(ee.Image(0).rename('increaseFlag'))
                                    .set('year',0)]);

//== Function to pass to Iterate(): create a new band where it is '1' when the current
//==  NBR value is higher than the previous NBR value
var flagIncrease = function(image, list) {
  var addImg = ee.Image.constant(1000);
  var currNBR = image.select(['NBR_fitted']).add(addImg);
  var prevImg = ee.Image(ee.List(list).get(-1));
  var prevNBR = prevImg.select(['NBR_fitted']).add(addImg);
  var incrFlag = ee.Image.constant(0).where(currNBR.gt(prevNBR), 1).rename('increaseFlag');
  var outImg = image.addBands(incrFlag).set('system:time_start', image.get('system:time_start'))
                    .set('year', image.get('year'));
  return ee.List(list).add(outImg);
};

//== Apply iteration to postDistColl
var incrFlagColl = ee.ImageCollection(ee.List(postDistColl.iterate(flagIncrease, startImg)));

// //== PRINT OUTPUTS
// print('incrFlagColl', incrFlagColl);

//== Remove first emtpy image from  collection (e.g. 'startImg')
var newFlagColl = incrFlagColl.filterMetadata('year','not_equals',0);

// //== PRINT OUTPUTS
// print('newFlagColl', newFlagColl);

//== Create multiband image from incrFlagColl image collection
var incrFlagImg = ee.Image(0);
for (var i = startYear; i <= endYear; i++){
  var img = ee.Image(newFlagColl.filterMetadata('year','equals',i).first());
  var bandName = ee.String(img.get('year'));
  var newBand = img.select(['increaseFlag']).rename(bandName);
  var incrFlagImg = incrFlagImg.addBands(newBand);
}
incrFlagImg = incrFlagImg.select(incrFlagImg.bandNames().removeAll(["constant"]));

//== Use 'increaseFlag' band to mask out 'Year' band in newFlagColl
var flagCollyrMasked = newFlagColl.map(function(image){
  var yrBand = image.select('Year');
  var incrFlag = image.select(['increaseFlag']);
  var flagMask = incrFlag.eq(1);
  var markedYr = ee.Image.constant(endYear + 1).where(incrFlag.eq(1), yrBand).rename('Year');
  return markedYr;
});

// //== PRINT OUTPUTS
// print('flagCollyrMasked', flagCollyrMasked);

//== Create multiband image from flagCollyrMasked image collection
var incrFlaggedYrImg = ee.Image(0);
for (var i = startYear; i <= endYear; i++){
  var img = ee.Image(flagCollyrMasked.filterMetadata('year','equals',i).first());
  var bandName = ee.String(img.get('year'));
  var newBand = img.select(['Year']).rename(bandName);
  var incrFlaggedYrImg = incrFlaggedYrImg.addBands(newBand);
}
incrFlaggedYrImg = incrFlaggedYrImg.select(incrFlaggedYrImg.bandNames().removeAll(["constant"]));

// //== PRINT OUTPUTS
// print('incrFlaggedYrImg', incrFlaggedYrImg);

//== Reduce collection, calculating minimum value for each band in each pixel,
//==  thereby extracting the earliest year where NBR values begin increasing from previous year
var extractEOC = flagCollyrMasked.select(['Year']).reduce(ee.Reducer.min());

// //== PRINT OUTPUTS
// print('extractEOC',extractEOC);

//== Extract the 'Year_min' band - the earliest year in which NBR increase from previous year
//==  is detected, and subtract 1 to get year of 'end of change'
var eocYrImage = extractEOC.select(['Year_min']).subtract(ee.Image.constant(1)).rename('eocYear');

//== Set any pixels where eoc was calculated as earlier than startYear, to the startYear value
var minYear = ee.Image.constant(startYear);
eocYrImage.where(eocYrImage.lt(startYear), minYear);

// //== PRINT OUTPUTS
// print('eocYrImage',eocYrImage);

//== Map function over each band (which represents fitted NBR for that year),
//==  that masks out pixels for which that year does not equal the value at 
//==  the end of all harvest change (i.e., eoc)
var eocMaskedStack = listBands.map(function(band){
    var bandYrImg = ee.Image.constant(ee.Number.parse(band));
    var eocMask = bandYrImg.eq(eocYrImage);
    var relvVals = ltFitStack.select([band]).multiply(distDir).toFloat()
    .mask(eocMask).set('year', ee.Number.parse(band)).rename('eocVal');
    return relvVals;
  });

//== Convert outputted list to image collection
var eocValColl = ee.ImageCollection(eocMaskedStack);

//== Combine image collection into one image (each pixel should be masked out in
//==  every image save one)
var eocValImg = eocValColl.reduce(ee.Reducer.median()).rename('endDistVal').toFloat();

// //== PRINT OUTPUTS
// print('eocValImg',eocValImg);


//-------------------------------------------------------------------
// ------- Length of Disturbance Period (Between YOD and EOC) -------    * WORKS! *

//== Calculate length of time between start of harvest and start of regeneration
var lenDistb = eocYrImage.subtract(yodImage).add(1);


//------------------------------------------------------------------
// ------- Length of Regeneration Period (endYear minus EOC) -------    * WORKS! *

//== Calculate length of regeneration period (EOC to endYear)
var lenRegen = ee.Image.constant(endYear).subtract(eocYrImage);


//--------------------------------------------------------------
// ------- Total Change (Relative to Pre-Distb Mean NBR) -------    * WORKS! *

//== Subtract value at end of disturbance from pre-disturbanc mean NBR
var totDistbVal = preDistMean.subtract(eocValImg);

// //== PRINT OUTPUTS
// print('totDistbVal',totDistbVal);


//----------------------------------------------------------------------
// ------- Relative Change (Relative to Full NBR Possible Range) -------    * WORKS! *

//== Calculate relative NBR change, based on Miller & Thode (2007; in RSE)
//== Formula: preDistbNBR - postDistbNBR / (SQRT(ABS(preDistbNBR/1000)))
var relDistbVal = totDistbVal.divide(preDistMean.divide(1000).abs().sqrt());


// //---------------------------------
// // ----- Rate of Regeneration -----    * WORKS! *

//== Function to mask out fitted NBR values that are before end of change
var postEocColl = fitCollection.map(function(image) {
  var currYr = ee.Number.parse(image.get('year'));
  var yrImg = ee.Image.constant(currYr);
  var postEoCDates = yrImg.gte(eocYrImage);
  var yearBand = ee.Image.constant(currYr).rename('Year').toFloat();
  return image.addBands(yearBand).updateMask(postEoCDates);
});

// //== PRINT OUTPUTS
// print('postEocColl',postEocColl);


//------------------------------------------
// ------- 5-Yr Percent Regeneration -------    * WORKING! *

//== Create image for year 5 years post-EOC
var yr5PostEocImg = eocYrImage.add(ee.Image.constant(5));

//== Map function over postEocColl to extract 5-year percent regeneration, based on EOC dates and values
var regen5yrColl = postEocColl.map(function(image) {
  var yr5flag = ee.Image(0).where(image.select('Year').eq(yr5PostEocImg), 1);
  var yr5vals = ee.Image(0).where(yr5flag.eq(1), image.select('NBR_fitted'));
  var yr5chgCalc = yr5vals.subtract(eocValImg);
  var yr5chg = yr5chgCalc.where(yr5flag.eq(0), 0);
  var regenPercCalc = yr5chg.divide(totDistbVal).multiply(100);
  regenPercCalc.where(yr5vals.eq(0), 0);
  return regenPercCalc.addBands(yr5vals)
                      .addBands(yr5flag)
                      .addBands(yr5chg)
                      .rename(['yr5regen','yr5val','yr5flag','yr5chg'])
                      .set('Year', image.get('year'));
  // return regenPercCalc.addBands(yr5vals).rename(['yr5regen','yr5val']);
});   

// //== PRINT OUTPUTS
// print('regen5yrColl', regen5yrColl);

//== Combine collection bands into single images
var regen5yrPerc = regen5yrColl.select('yr5regen').reduce(ee.Reducer.max()).toFloat();
var yr5postEocVal = regen5yrColl.select('yr5val').reduce(ee.Reducer.max()).toFloat();
var yr5postEocFlag = regen5yrColl.select('yr5flag').reduce(ee.Reducer.max()).toFloat();
var yr5postEocChg = regen5yrColl.select('yr5chg').reduce(ee.Reducer.max()).toFloat();

// //== PRINT OUTPUTS
// print('regen5yrPerc', regen5yrPerc);
// print('yr5postEocVal', yr5postEocVal)

//== Create a mask, to keep those where 5 year regen metrics are calculable, and mask out others;
//==  apply it
// var yr5RegenMask = yr5postEocVal.neq(0);
var yr5RegenMask = yr5postEocFlag.neq(0);
yr5RegenMask = yr5RegenMask.updateMask(yr5RegenMask);

regen5yrPerc = regen5yrPerc.updateMask(yr5RegenMask);
yr5postEocVal = yr5postEocVal.updateMask(yr5RegenMask);
yr5postEocFlag = yr5postEocFlag.updateMask(yr5RegenMask);
yr5postEocChg = yr5postEocChg.updateMask(yr5RegenMask);

//== Create a reverse mask, to flag those pixels where 5 year regen is not/cannot be recalculated,
//==  and add it to the maskFlagsImg image
// var noYr5RegenMask = yr5RegenMask.not(); 
var noYr5RegenMask = ee.Image.constant(1).where(yr5RegenMask.eq(1),0); //**WORKS
noYr5RegenMask = noYr5RegenMask.updateMask(noYr5RegenMask);
var maskFlagsImg = ee.Image.constant(0).where(noYr5RegenMask.eq(1),1).rename('No5yrRegen');


//-------------------------------------------
// ------- 10-Yr Percent Regeneration -------    * WORKING! *

//== Create image for year 10 years post-EOC
var yr10PostEocImg = eocYrImage.add(ee.Image.constant(10));

//== Map function over postEocColl to extract 10-year percent regeneration, based on EOC dates and values
var regen10yrColl = postEocColl.map(function(image) {
  var yr10flag = ee.Image(0).where(image.select('Year').eq(yr10PostEocImg), 1);
  var yr10vals = ee.Image(0).where(yr10flag.eq(1), image.select('NBR_fitted'));
  var yr10chgCalc = yr10vals.subtract(eocValImg);
  var yr10chg = yr10chgCalc.where(yr10flag.eq(0), 0);
  var regenPercCalc = yr10chg.divide(totDistbVal).multiply(100);
  regenPercCalc.where(yr10vals.eq(0), 0);
  return regenPercCalc.addBands(yr10vals)
                      .addBands(yr10flag)
                      .addBands(yr10chg)
                      .rename(['yr10regen','yr10val','yr10flag','yr10chg'])
                      .set('Year', image.get('year'));
  // return regenPercCalc.addBands(yr10vals).rename(['yr10regen','yr10val']);
});   //== ** New version, where inappropirate results are given 0, not masked out  *NOT WORKING :(

// //== PRINT OUTPUTS
// print('regen10yrColl', regen10yrColl);

//== Combine collection bands into single images
var regen10yrPerc = regen10yrColl.select('yr10regen').reduce(ee.Reducer.max()).toFloat();
var yr10postEocVal = regen10yrColl.select('yr10val').reduce(ee.Reducer.max()).toFloat();
var yr10postEocFlag = regen10yrColl.select('yr10flag').reduce(ee.Reducer.max()).toFloat();
var yr10postEocChg = regen10yrColl.select('yr10chg').reduce(ee.Reducer.max()).toFloat();

// //== PRINT OUTPUTS
// print('regen10yrPerc', regen10yrPerc);

//== Create a mask, to keep those where 10 year regen metrics are calculable, and mask out others;
//==  apply it
// var yr10RegenMask = yr10postEocVal.neq(0);
var yr10RegenMask = yr10postEocFlag.neq(0);
yr10RegenMask = yr10RegenMask.updateMask(yr10RegenMask); 

regen10yrPerc = regen10yrPerc.updateMask(yr10RegenMask);
yr10postEocVal = yr10postEocVal.updateMask(yr10RegenMask);
yr10postEocFlag = yr10postEocFlag.updateMask(yr10RegenMask);
yr10postEocChg = yr10postEocChg.updateMask(yr10RegenMask);

//== Create a reverse mask, to flag those pixels where 10 year regen is not/cannot be recalculated,
//==  and add it to the maskFlagsImg image
// var noYr10RegenMask = yr10RegenMask.not(); 
var noYr10RegenMask = ee.Image.constant(1).where(yr10RegenMask.eq(1),0); //**WORKS
noYr10RegenMask = noYr10RegenMask.updateMask(noYr10RegenMask);
var maskFlagsImg = maskFlagsImg.addBands(ee.Image.constant(0).where(noYr10RegenMask.eq(1),1).rename('No10yrRegen'));


//-------------------------------------------
// ------- 15-Yr Percent Regeneration -------    * WORKING! *

//== Create image for year 15 years post-EOC
var yr15PostEocImg = eocYrImage.add(ee.Image.constant(15));

//== Map function over postEocColl to extract 15-year percent regeneration, based on EOC dates and values
var regen15yrColl = postEocColl.map(function(image) {
  var yr15flag = ee.Image(0).where(image.select('Year').eq(yr15PostEocImg), 1);
  var yr15vals = ee.Image(0).where(yr15flag.eq(1), image.select('NBR_fitted'));
  var yr15chgCalc = yr15vals.subtract(eocValImg);
  var yr15chg = yr15chgCalc.where(yr15flag.eq(0), 0);
  var regenPercCalc = yr15chg.divide(totDistbVal).multiply(100);
  regenPercCalc.where(yr15vals.eq(0), 0);
  return regenPercCalc.addBands(yr15vals)
                      .addBands(yr15flag)
                      .addBands(yr15chg)
                      .rename(['yr15regen','yr15val','yr15flag','yr15chg'])
                      .set('Year', image.get('year'));
  // return regenPercCalc.addBands(yr15vals).rename(['yr15regen','yr15val']);
});   //== ** New version, where inappropirate results are given 0, not masked out  *NOT WORKING :(

// //== PRINT OUTPUTS
// print('regen15yrColl', regen15yrColl);

//== Combine collection bands into single images
var regen15yrPerc = regen15yrColl.select('yr15regen').reduce(ee.Reducer.max()).toFloat();
var yr15postEocVal = regen15yrColl.select('yr15val').reduce(ee.Reducer.max()).toFloat();
var yr15postEocFlag = regen15yrColl.select('yr15flag').reduce(ee.Reducer.max()).toFloat();
var yr15postEocChg = regen15yrColl.select('yr15chg').reduce(ee.Reducer.max()).toFloat();

// //== PRINT OUTPUTS
// print('regen15yrPerc', regen15yrPerc);

//== Create a mask, to keep those where 15 year regen metrics are calculable, and mask out others;
//==  apply it
// var yr15RegenMask = yr15postEocVal.neq(0);
var yr15RegenMask = yr15postEocFlag.neq(0);
yr15RegenMask = yr15RegenMask.updateMask(yr15RegenMask);

regen15yrPerc = regen15yrPerc.updateMask(yr15RegenMask);
yr15postEocVal = yr15postEocVal.updateMask(yr15RegenMask);
yr15postEocFlag = yr15postEocFlag.updateMask(yr15RegenMask);
yr15postEocChg = yr15postEocChg.updateMask(yr15RegenMask);

//== Create a reverse mask, to flag those pixels where 15 year regen is not/cannot be recalculated,
//==  and add it to the maskFlagsImg image
// var noYr15RegenMask = yr15RegenMask.not(); 
var noYr15RegenMask = ee.Image.constant(1).where(yr15RegenMask.eq(1),0); //**WORKS
noYr15RegenMask = noYr15RegenMask.updateMask(noYr15RegenMask);
var maskFlagsImg = maskFlagsImg.addBands(ee.Image.constant(0).where(noYr15RegenMask.eq(1),1).rename('No15yrRegen'));


//-------------------------------------------
// ------- 20-Yr Percent Regeneration -------    * WORKING! *

//== Create image for year 20 years post-EOC
var yr20PostEocImg = eocYrImage.add(ee.Image.constant(20));

//== Map function over postEocColl to extract 20-year percent regeneration, based on EOC dates and values
var regen20yrColl = postEocColl.map(function(image) {
  var yr20flag = ee.Image(0).where(image.select('Year').eq(yr20PostEocImg), 1);
  var yr20vals = ee.Image(0).where(yr20flag.eq(1), image.select('NBR_fitted'));
  var yr20chgCalc = yr20vals.subtract(eocValImg);
  var yr20chg = yr20chgCalc.where(yr20flag.eq(0), 0);
  var regenPercCalc = yr20chg.divide(totDistbVal).multiply(100);
  regenPercCalc.where(yr20vals.eq(0), 0);
  return regenPercCalc.addBands(yr20vals)
                      .addBands(yr20flag)
                      .addBands(yr20chg)
                      .rename(['yr20regen','yr20val','yr20flag','yr20chg'])
                      .set('Year', image.get('year'));
  // return regenPercCalc.addBands(yr20vals).rename(['yr20regen','yr20val']);
});   //== ** New version, where inappropirate results are given 0, not masked out  *NOT WORKING :(

// //== PRINT OUTPUTS
// print('regen20yrColl', regen20yrColl);

//== Combine collection bands into single images
var regen20yrPerc = regen20yrColl.select('yr20regen').reduce(ee.Reducer.max()).toFloat();
var yr20postEocVal = regen20yrColl.select('yr20val').reduce(ee.Reducer.max()).toFloat();
var yr20postEocFlag = regen20yrColl.select('yr20flag').reduce(ee.Reducer.max()).toFloat();
var yr20postEocChg = regen20yrColl.select('yr20chg').reduce(ee.Reducer.max()).toFloat();

// //== PRINT OUTPUTS
// print('regen20yrPerc', regen20yrPerc);

//== Create a mask, to keep those where 20 year regen metrics are calculable, and mask out others;
//==  apply it
// var yr20RegenMask = yr20postEocVal.neq(0);
var yr20RegenMask = yr20postEocFlag.neq(0);
yr20RegenMask = yr20RegenMask.updateMask(yr20RegenMask);

regen20yrPerc = regen20yrPerc.updateMask(yr20RegenMask);
yr20postEocVal = yr20postEocVal.updateMask(yr20RegenMask);
yr20postEocFlag = yr20postEocFlag.updateMask(yr20RegenMask);
yr20postEocChg = yr20postEocChg.updateMask(yr20RegenMask);

//== Create a reverse mask, to flag those pixels where 20 year regen is not/cannot be recalculated,
//==  and add it to the maskFlagsImg image
// var noYr20RegenMask = yr20RegenMask.not();
var noYr20RegenMask = ee.Image.constant(1).where(yr20RegenMask.eq(1),0); //**WORKS
noYr20RegenMask = noYr20RegenMask.updateMask(noYr20RegenMask);
var maskFlagsImg = maskFlagsImg.addBands(ee.Image.constant(0).where(noYr20RegenMask.eq(1),1).rename('No20yrRegen'));



//-----------------------------------
// ------- Total Regeneration -------    * WORKING! *

//== Extract index of last band in ltNBRbands
var ltNBRbands = ltNBRFitted.bandNames();
var lastBandIdx = ltNBRbands.length().subtract(1);
var lastBand = ltNBRbands.get(lastBandIdx);

//== Extract final band in NBR fitted time series
var endTSval = ltNBRFitted.select([lastBand]);

// //== PRINT OUTPUT
// print('endTSval',endTSval);

//== Calculate total percent regeneration (regen / disturbance)
var totPercRegen = endTSval.subtract(eocValImg).divide(totDistbVal)
                          .multiply(100).rename('percTotRegen');

// //== PRINT OUTPUT
// print('totPercRegen',totPercRegen);



//---------------------------------------------------
// ------- Years to 30% Regeneration (Y2R30%) -------  ** NEW **

//== Create image of NBR values at 30% recovered
var img30perc = totDistbVal.multiply(0.3).add(eocValImg);

//== Map function over postEocColl to extract years to 30% spectrally regenerated
var y2r30Coll = postEocColl.map(function(image){
  var currYr = image.select('Year');
  var at30regen = ee.Image(0)
    .where(image.select('NBR_fitted').gte(img30perc), currYr)
    .rename('y2r30');
  return at30regen.set('year', image.get('year')).updateMask(at30regen.gt(0));
});
// print('y2r30Coll', y2r30Coll);

//== Extract the minimum value per pixel, from the collection
var yrAt30img = y2r30Coll.min();

//== Calulate years it took to reach 30% regenerated
var y2r30img = yrAt30img.subtract(eocYrImage).rename('y2r30');
y2r30img.updateMask(y2r30img.gt(0));
// print('y2r30img',y2r30img);


//---------------------------------------------------
// ------- Years to 50% Regeneration (Y2R50%) -------  ** NEW **

//== Create image of NBR values at 50% recovered
var img50perc = totDistbVal.multiply(0.5).add(eocValImg);

//== Map function over postEocColl to extract years to 50% spectrally regenerated
var y2r50Coll = postEocColl.map(function(image){
  var currYr = image.select('Year');
  var at50regen = ee.Image(0)
    .where(image.select('NBR_fitted').gte(img50perc), currYr)
    .rename('y2r50');
  return at50regen.set('year', image.get('year')).updateMask(at50regen.gt(0));
});
// print('y2r50Coll', y2r50Coll);

//== Extract the minimum value per pixel, from the collection
var yrAt50img = y2r50Coll.min();

//== Calulate years it took to reach 50% regenerated
var y2r50img = yrAt50img.subtract(eocYrImage).rename('y2r50');
y2r50img.updateMask(y2r50img.gt(0));
// print('y2r50img',y2r50img);


//---------------------------------------------------
// ------- Years to 80% Regeneration (Y2R80%) -------  ** NEW **

//== Create image of NBR values at 80% recovered
var img80perc = totDistbVal.multiply(0.8).add(eocValImg);

//== Map function over postEocColl to extract years to 80% spectrally regenerated
var y2r80Coll = postEocColl.map(function(image){
  var currYr = image.select('Year');
  var at80regen = ee.Image(0)
    .where(image.select('NBR_fitted').gte(img80perc), currYr)
    .rename('y2r80');
  return at80regen.set('year', image.get('year')).updateMask(at80regen.gt(0));
});
// print('y2r80Coll', y2r80Coll);

//== Extract the minimum value per pixel, from the collection
var yrAt80img = y2r80Coll.min();

//== Calulate years it took to reach 80% regenerated
var y2r80img = yrAt80img.subtract(eocYrImage).rename('y2r80');
y2r80img.updateMask(y2r80img.gt(0));
// print('y2r80img',y2r80img);


//-----------------------------------------------------
// ------- Years to 100% Regeneration (Y2R100%) -------  ** NEW **

//== Create image of NBR values at 100% recovered
var img100perc = totDistbVal.add(eocValImg);

//== Map function over postEocColl to extract years to 100% spectrally regenerated
var y2r100Coll = postEocColl.map(function(image){
  var currYr = image.select('Year');
  var at100regen = ee.Image(0)
    .where(image.select('NBR_fitted').gte(img100perc), currYr)
    .rename('y2r100');
  return at100regen.set('year', image.get('year')).updateMask(at100regen.gt(0));
});
// print('y2r100Coll', y2r100Coll);

//== Extract the minimum value per pixel, from the collection
var yrAt100img = y2r100Coll.min();

//== Calulate years it took to reach 100% regenerated
var y2r100img = yrAt100img.subtract(eocYrImage).rename('y2r100');
y2r100img.updateMask(y2r100img.gt(0));
// print('y2r100img',y2r100img);



//#############################################################################
//######## MASKING (OF IRRELEVANT OR INAPPROPRIATE PIXELSFOR ANALYSIS) ########
//#############################################################################
 

//--------------------------------------------
// ------- Minimum Disturbance Masking -------    * WORKS! *

//== Create mask for masking out those pixels where the total disturbance is below a set 
//==  threshold (indicating that there was not reliably identifiable harvest happening in that pixel)
// var harvThreshMask = totDistbVal.gte(200); //** in NBR units, based on personal observations // *WORKS
var harvThreshMask = totDistbVal.gte(175); //** in NBR units, based on personal observations // *WORKS

//== Add to maskFlagsImg an image with flags where the harvest threshold mask was used (no harvest)
var reverseHarvMask = harvThreshMask.not(); // *WORKS
var maskFlagsImg = maskFlagsImg.addBands(ee.Image.constant(0).where(reverseHarvMask.eq(1),1).rename('NoHarvest'));


//-----------------------------------------
// ------- Out-f-Date-Range Masking -------    * WORKS! *  

//== Create mask to remove those pixels where the greatest change (yod) is withihn 5 years of the start
//== of the time series, and eoc is <5 years before the end of the time series
var minCutoff = ee.Image.constant(startYear + 5);
var maxCutoff = ee.Image.constant(endYear - 5);
// var relvDatesMask = yodImage.gte(minCutoff).and(yodImage.lte(maxCutoff)); // *WORKS 
var relvDatesMask = yodImage.gte(minCutoff).and(eocYrImage.lte(maxCutoff));

//== Create image with flags for pixels masked because of yod date range, and add to maskFlagsImg
var reverseRelMask = relvDatesMask.not(); // *WORKS
var maskFlagsImg = maskFlagsImg.addBands(ee.Image.constant(0).where(reverseRelMask.eq(1),1).rename('YODOutOfRange'));


//----------------------------------------
// ------- No Regeneration Masking -------    * WORKS! *

//== Create mask to ID those pixels where EOC equals endYear (2017); EOC was
//==  set to endYear in those cases where no post-YOD increases were found
var regenMask = eocYrImage.neq(endYear);

//== Create image band where pixels with no regeneration are flagged
var noRegenMask = regenMask.not(); // *WORKS
var maskFlagsImg = maskFlagsImg.addBands(ee.Image.constant(0).where(noRegenMask.eq(1),1).rename('NoRegenDetected'));


//----------------------------------------------
// ------- Mulitiple Disturbance Masking -------     * WORKS! *

//== Extract 'mag' and 'yod' for largest and second-largest drops in NBR in the time series,
//== to identify where there may be multiple disturbance events in a time series
var magImage = distImg.select(["mag"]);
var mag2ndImg = distImg.select(["mag2nd"]);
var yod2ndImg = distImg.select(["yod2nd"]);

//== Create masks for: where 2nd largest disturbance is > 175 NBR, and for dates
//==  that are > 3 years before YOD and > 3 years after EOC
var mag2gt175 = mag2ndImg.gt(175);
var yodMinus3yrs = yodImage.subtract(ee.Image.constant(3));
var eocPlus3yrs = eocYrImage.add(ee.Image.constant(3));
var yod2ndInRange = yod2ndImg.lt(yodMinus3yrs).or(yod2ndImg.gt(eocPlus3yrs));

//== Combine above masks to flag pixels in whic the 2nd largest disturbance is > 175 NBR,
//==  and is not within 3 years of 1st disturbance (YOD to EOC)
var multiDistbPix = mag2gt175.eq(1).and(yod2ndInRange.eq(1));

//== Create reverse of above and add to flagging image
var noMultiDistb = multiDistbPix.not();
var maskFlagsImg = maskFlagsImg.addBands(ee.Image.constant(0).where(multiDistbPix.eq(1),1).rename('MultiDistb'));

// //== PRINT OUTPUTS
// print('maskFlagsImg',maskFlagsImg);

//-----------------------------------
// ------- Combine Flag Image -------     

var updtMaskFlagsImg = maskFlagsImg.updateMask(maskFlagsImg);

//== Reorder the bands in maskFlagsImg so that the most important flags are at the end
//==  Order: No20yrRegen, No15yrRegen, No10yrRegen, No5yrRegen, NoRegenDetected, YODOutOfRange, MultiDistb, NoHarvest
var maskFlagsImgReordered = maskFlagsImg.select(['No20yrRegen','No15yrRegen', 'No10yrRegen', 'No5yrRegen', 'NoRegenDetected', 
                                                'YODOutOfRange', 'MultiDistb', 'NoHarvest']);

//== Reorder and select key flag mask bands; exclude those flags where metric is not calculated because period of 
//==  regeneration is not long enough
var maskKeyFlagsImgReordered = maskFlagsImg.select(['NoRegenDetected', 
                                                'YODOutOfRange', 'MultiDistb', 'NoHarvest']);
var keyFlagExists = maskKeyFlagsImgReordered.reduce(ee.Reducer.anyNonZero());
var zeroKeyFlagExists = ee.Image.constant(1).updateMask(keyFlagExists.not());


//########################################################################################################
//########## CALCULATE HARVEST AREA POLYGON ZONAL STATISTICS FOR POST-DISTURBANCE METRICS/INFO ###########
//########################################################################################################

//-------------------------------------------------------------------
// ------- Import Feature Collection of Harvest Area Polygons -------
 
//// print(HA)

// ***NEW for 2018 HAs (2020-06-17) * Seems to Work *
//== Function to add "OBJECTID" field (which doesn't exist in this fusion table, yet)
HA = HA.map(function(feat){
  return feat.set('OBJECTID', feat.get('uniquID'));
});

//== PRINT OUTPUTS
// print('HA polygon (first)', HA.first());
// print('HA polygon', HA)

//== Example HA polygon object, for visualization
// var featID = 206906;
var featID = ee.Number(HA.first().get('OBJECTID'));
// print('featID =', featID);

//== In new fusion table, OBJECTID is replaced with FIRST_OBJE
// var haObject = ee.Feature(HA.filterMetadata("FIRST_OBJE", "equals", featID).first());
var haObject = ee.Feature(HA.filterMetadata('OBJECTID', 'equals', featID).first());



//== PRINT OUTPUTS
// print('haObject', haObject);

// get a landsat collection for given year, day range, and sensor
var tempColl = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR') // get L5 surface reflectance images
                      .filterBounds(aoi)                                  // ...filter them by intersection with AOI
                      .filterDate(startYear+'-'+startDay, startYear+'-'+endDay);

// extract projection
var origProj = ee.Image(tempColl.first()).projection();


//---------------------------------------------------------------------------------------------------------
// ------- Extract Numbers of Pixels (Incl. Partials) per HA Polygon, for Unmasked & Masked Outputs -------

//== Create blank image filled with 1s, and then create masked version
var constImg = ee.Image.constant(1).rename('allNoPix');
var constImgMskd = constImg.updateMask(keyFlagExists.neq(1)).rename('unMskPix');

//== Use reducer to extract sum of all pixels per HA polygon
var sumAllPix = constImg.reduceRegions({
  collection: HA,
  reducer: ee.Reducer.sum(),
  crs: origProj,
  scale: 30,
  tileScale: 8
});
// print('sumAllPix', sumAllPix);


//== Use reducer to extract sum of unmasked (kept) pixels per HA polygon
var sumUnMskPix = constImgMskd.reduceRegions({
  collection: HA,
  reducer: ee.Reducer.sum(),
  crs: origProj,
  scale: 30,
  tileScale: 8
});
// print('sumUnMskPix', sumUnMskPix);

var sumZeroFlags = zeroKeyFlagExists.reduceRegions({
  collection: HA,
  reducer: ee.Reducer.sum(),
  crs: origProj,
  scale: 30,
  tileScale: 8
});
// print('sumZeroFlags', sumZeroFlags);

//-------------------------------------------------------------
// ------- Extract Flags (Causes of Masking) Information -------

//== Use reducer to calculate per-HA set-bucket histogram of masking flags
var AltCustmHistFlags = maskFlagsImgReordered.reduceRegions({
  collection: HA,
  reducer: ee.Reducer.fixedHistogram({
    min:0, max:9, steps:9 // bingo...
  }),
  crs: origProj,
  scale: 30,
  tileScale: 8
});
// print('AltCustmHistFlags',AltCustmHistFlags);


//-------------------------------------------------------------------------------
// ------- Calculate & Extract Number of Contiguous Pixels per HA Feature -------

//== Clip layer identifying unmasked pixels (not flagged at all) to polygon feature collection
var clipConstImgMskd = constImgMskd.clip(HA);

//== Create binary clipped image of pixels used in calculations (where weights are all 1)
var usedPixels = ee.Image(1).updateMask(clipConstImgMskd.gt(0));

//== Calculate per-pixel 'connected pixel count' that includes only those pixels attached on sides (not diagonally)
var contigPixImg = usedPixels.connectedPixelCount(11,false);

//== Use max() reducer to extract maximum "number of connected pixels" per each polygon
var maxContigPix = contigPixImg.reduceRegions({
  collection: HA,
  reducer: ee.Reducer.max(),
  crs: origProj,
  scale: 30,
  tileScale: 8
});
// print('maxContigPix', maxContigPix);
// print('maxContigPix.limit(4)', maxContigPix.limit(4));


//------------------------------------------------------------
// ------- Join Together All Metric Outputs for Export -------

//== Combine reducers (mean & std dev)
var twoReducers = ee.Reducer.mean().combine({
  reducer2: ee.Reducer.stdDev(),
  sharedInputs: true
});
// print('twoReducers',twoReducers); 

//== Combine output metric images into a collection
var outputColl = ee.ImageCollection.fromImages([
                                                preDistMean.rename('preDistMean').set('metricLayer','preDistMean'), 
                                                yodImage.rename('yod').set('metricLayer','yod'), 
                                                eocYrImage.rename('eocYr').set('metricLayer','eocYr'),
                                                lenDistb.rename('lengthDistb').set('metricLayer','lengthDistb'),
                                                lenRegen.rename('lengthRegen').set('metricLayer','lengthRegen'),
                                                totDistbVal.rename('totDistb').set('metricLayer','totDistb'), 
                                                relDistbVal.rename('relDistbVal').set('metricLayer','relDistbVal'),
                                                regen5yrPerc.rename('regen5yr').set('metricLayer','regen5yr'), 
                                                regen10yrPerc.rename('regen10yr').set('metricLayer','regen10yr'), 
                                                regen15yrPerc.rename('regen15yr').set('metricLayer','regen15yr'), 
                                                regen20yrPerc.rename('regen20yr').set('metricLayer','regen20yr'), 
                                                totPercRegen.rename('totRegen').set('metricLayer','totRegen'), 
                                                y2r30img.rename('yrTo30').set('metricLayer','yrTo30'),
                                                y2r50img.rename('yrTo50').set('metricLayer','yrTo50'),
                                                y2r80img.rename('yrTo80').set('metricLayer','yrTo80'),
                                                y2r100img.rename('yrTo100').set('metricLayer','yrTo100')
                                              ]); 
//// print('outputColl', outputColl);

//== Apply masks to keep only those pixels that are relevant/appropriate for extracting
//==  post-harvest regeneration info (harvest/disturbance occurs, occurs only once, is within appropriate
//==  date range, and regeneration is detected)
var outputCollmakd = outputColl.map(function(img){
  return img.updateMask(regenMask).updateMask(relvDatesMask)
            .updateMask(harvThreshMask).updateMask(noMultiDistb);
});

//== PRINT OUTPUTS
// print('outputCollmakd',outputCollmakd);


//#############################################################################
//#############################################################################


//== Function to extract zonal stats for metric images in collection
var extractZonalStats = outputCollmakd.map(function(img){  
  //== Extract name of img, create custom output names
  var imgName = ee.String(img.bandNames().get(0));
  var meanName = ee.String(imgName.cat("_mean"));
  var stdName = ee.String(imgName.cat("_stdev"));
  
  //== Combine reducers (mean & std dev)
  var twoReducers = ee.Reducer.mean().combine({
    reducer2: ee.Reducer.stdDev(),
    sharedInputs: true
  });
  
  //== Apply reducers to img
  var haStats = img.reduceRegions({
    collection: HA,
    reducer: twoReducers,
    crs: origProj,
    scale: 30,
    tileScale: 8
  });
  
  return haStats.map(function(poly){
    return poly.set('metric',imgName)
    .select({ 
              //== TO USE WITH FULL SET OF ALREADY CLIPPED POLYS:
              propertySelectors:['OBJECTID','mean','stdDev','metric'],
              newProperties:['OBJECTID','METRIC_MEAN','METRIC_STDEV','METRIC']});
  });

});
// print('extractZonalStats',extractZonalStats);


//== Flatten collection of collections
var outputMetrics = extractZonalStats.flatten();
print('outputMetrics.limit(4)',outputMetrics.limit(4));

//== PRINT OUTPUTS
// print('outputMetrics',outputMetrics);

//== Map function over outputMetric feature collection to add number of total pixels
//==  per polygon, and number of unmasked pixels
var metsForExport = outputMetrics.map(function(poly){

  var ID = ee.Number.parse(poly.get('OBJECTID'));
  var allPixPoly = sumAllPix.filterMetadata('OBJECTID','equals',ID).first();
  var mskPixPoly = sumUnMskPix.filterMetadata('OBJECTID','equals',ID).first();
  var zeroFlagPoly = sumZeroFlags.filterMetadata('OBJECTID','equals',ID).first();
  var flagPoly = AltCustmHistFlags.filterMetadata('OBJECTID','equals',ID).first();
  var maxCnnctdPixPoly = maxContigPix.filterMetadata('OBJECTID','equals',ID).first();
  

  //var casID = (HA.filterMetadata('OBJECTID','equals',ID).first());
  
  var sumAllVal = ee.Number(allPixPoly.get('sum'));
  var sumMskdVal = ee.Number(mskPixPoly.get('sum'));
  var sum0Flags = ee.Number(zeroFlagPoly.get('sum'));
  var maxCnnctdPix = ee.Number(maxCnnctdPixPoly.get('max'));
  
  var sum5yrFlags = ee.Array(flagPoly.get('No5yrRegen')).get([1,1]);
  var sum10yrFlags = ee.Array(flagPoly.get('No10yrRegen')).get([1,1]);
  var sum15yrFlags = ee.Array(flagPoly.get('No15yrRegen')).get([1,1]);
  var sum20yrFlags = ee.Array(flagPoly.get('No20yrRegen')).get([1,1]);
  var sumNoRegFlags = ee.Array(flagPoly.get('NoRegenDetected')).get([1,1]);
  var sumYOODFlags = ee.Array(flagPoly.get('YODOutOfRange')).get([1,1]);
  var sumMDFlags = ee.Array(flagPoly.get('MultiDistb')).get([1,1]);
  var sumNoHarvFlags = ee.Array(flagPoly.get('NoHarvest')).get([1,1]);
  
  

  //== Add properties to poly, using conditional if statements within function 
  return poly.set({
    TOT_PIXELS: sumAllVal, 
    UNMSKD_PIXELS: sumMskdVal, 
    PERC_UNMSKD: sumMskdVal.divide(sumAllVal).multiply(100),
    MAX_CONNCTDPIX: maxCnnctdPix,
    PERC_NOHARVEST: ee.Algorithms.If(sumAllVal.gt(0), sumNoHarvFlags.divide(sumAllVal).multiply(100),0),
    PERC_MULTIDISTB: ee.Algorithms.If(sumAllVal.gt(0), sumMDFlags.divide(sumAllVal).multiply(100),0),
    PERC_OUTOFRANGE: ee.Algorithms.If(sumAllVal.gt(0), sumYOODFlags.divide(sumAllVal).multiply(100),0),
    PERC_NOREGEN: ee.Algorithms.If(sumAllVal.gt(0), sumNoRegFlags.divide(sumAllVal).multiply(100),0),
    PERC_NO5REGEN: ee.Algorithms.If(sumAllVal.gt(0), sum5yrFlags.divide(sumAllVal).multiply(100),0),
    PERC_NO10REGEN: ee.Algorithms.If(sumAllVal.gt(0), sum10yrFlags.divide(sumAllVal).multiply(100),0),
    PERC_NO15REGEN: ee.Algorithms.If(sumAllVal.gt(0), sum15yrFlags.divide(sumAllVal).multiply(100),0),
    PERC_NO20REGEN: ee.Algorithms.If(sumAllVal.gt(0), sum20yrFlags.divide(sumAllVal).multiply(100),0)
    //cas_id: cas_id.select('cas_id').first()
  });
});

// -------------------------------------------------------

// Join NBR metrics with CASFRI data (CASID, station ID, survey year)

// Use an equals filter to specify how the collections match.
var casFilter = ee.Filter.equals({
  leftField: 'OBJECTID',
  rightField: 'OBJECTID'
});

// Define the join.
var innerJoin = ee.Join.inner();

// Apply the join.
var metsForExportfin = innerJoin.apply(metsForExport, HA, casFilter);


// What join does is it returns features where your join condition match and put the matching 
// features as properties named 'primary' and 'secondary' of a feature with null geometry. All that is left is to 
// reformat it to a proper feature collection with properties from both feature that I like to call cleaning the join. 

function cleanJoin(metsForExportfin){
  return ee.Feature(metsForExportfin.get('primary')).copyProperties(metsForExportfin.get('secondary'));
}

metsForExportfin  = metsForExportfin.map(cleanJoin);

// Print the result.
print('metsForExportfin.limit(4)',metsForExportfin.limit(4));

//########################################################################################################
//######################################### EXPORT OUTPUT TO DRIVE #######################################
//########################################################################################################

//== ----- Export feature collection results -----

//== Export to google drive 
Export.table.toDrive({
  collection: metsForExportfin,
  description: 'HA_Polys_ABMI'+ endYear,
  fileFormat: 'CSV',
  folder: 'HA_Polys',
  selectors: 'OBJECTID, HFI_ID, uniquID, YEAR, MAX_CONNCTDPIX,METRIC,METRIC_MEAN,METRIC_STDEV,PERC_MULTIDISTB,PERC_NO10REGEN,PERC_NO15REGEN, PERC_NO20REGEN,PERC_NO5REGEN,PERC_NOHARVEST,PERC_NOREGEN,PERC_OUTOFRANGE,PERC_UNMSKD,TOT_PIXELS,UNMSKD_PIXELS'
});


//----------------------------------------------------------------------------
// -------------- Add Harvest Area Polygons with Buffers to Map --------------
//----------------------------------------------------------------------------

// // == Add HA polygons to map
// Map.addLayer(HA, {color:'FFA0A0'}, 'Harvest Areas');

// // == Add single HA object to map
// Map.addLayer(haObject, {color:'AAAAAA'/*'FFFFFF'*/}, 'Harvest Area Polygon');


