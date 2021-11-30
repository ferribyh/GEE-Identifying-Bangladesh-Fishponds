//This is the GEE code for inland fishpond identification and classification
//in seven districts of Bangladesh (Bagerhat, Barisal, Bhola, Gopalganj, Jessore,
//Khulna, and Satkhira)

//Original code created by Zhiqi Yu denoted with a *Yu beforehand
//Adapted by Hannah Ferriby (ferribyh@msu.edu)


//INPUTS! 
//edited_region -> study region shapefiles
//gtMerge -> Final ground truthing shapefiles


//Functions
function toNDWI(image) {
  //for sentinel B8 and B3
  return image.normalizedDifference(['B3', 'B8']).rename('NDWI'); }
  
  
function toMNDWI1(image) {
  //for sentinel B11 and B3
  return image.normalizedDifference(['B3', 'B11']).rename('MNDWI'); }


function toAWEInsh(image) {
return image.expression(
    '4 * (GREEN - SWIR1) - (0.25 * NIR) - (2.75 * SWIR2)', {
    'GREEN': image.select('B3'),
    'SWIR1': image.select('B11'),
    'NIR': image.select('B8'),
    'SWIR2': image.select('B12')
  } ).rename('AWEInsh');
   }

//*Yu
function getIndexSeq(imgCollection, name, bands){
  var ndCollection = imgCollection.map(function(img){
    var nd = img.normalizedDifference(bands).rename(ee.String(img.get('system:index')).cat('_'+name));
    return nd;
  });
  return ndCollection.toBands();
}

//*Yu
function otsu(hist) {
  var counts = ee.Array(ee.Dictionary(hist).get('histogram'));
  var means = ee.Array(ee.Dictionary(hist).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });

  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
}

//*Yu
function otsu_mask(image, area, mask, bins, scale, is_gt) {
  var histogram = image.mask(mask).reduceRegion({
  reducer: ee.Reducer.histogram(bins, 0.001),
    geometry: area,
    scale: scale,
    maxPixels: 1e13
  });
  print("Histogram", histogram)
  
  histogram = histogram.map(function (key, value) {
    var hist = ee.Dictionary(value);
    if (hist.contains(ee.String('bucketMeans'))) {
      return value
    }
    return null
  })
  var bands = image.bandNames();
  print("Histogram filtered", histogram)
  print("Bands", bands)
  var bands_filtered = histogram.keys();
  
  print("Bands_filtered ", bands_filtered)
  var masks = bands_filtered.map(function(name) {
    var key = ee.String(name);
    var hist = ee.Dictionary(histogram.get(key));
  
    // compute threshold for the current band
    var thres = otsu(hist);
      
    // use threshold to segment the band
    var band = image.select(ee.String(name));
    if(is_gt){
       return band.gte(thres);
      }
    else{
      return band.lt(thres);
      }
    });
  print("Masks", masks)
  var water_mask = ee.ImageCollection(masks).toBands().reduce(ee.Reducer.allNonZero());
  return water_mask;
}

//*Yu
function add_spatial_attr(featureCollec, bound){
  // Add object-based spatial features to features
  var water_feature = featureCollec.map(function(f){
    var perimeter = f.perimeter(1);
    var area = f.area(1);
    var convex_hull = f.convexHull(1);
    var convex_area = convex_hull.area(1);
    var convex_perim = convex_hull.perimeter(1);
    var shaped = f.set('area', area)
            .set('perimeter', perimeter)
            .set('ipq', area.multiply(4*Math.PI).divide(perimeter.pow(2)))
            .set('soli', area.divide(convex_area))
            .set('pfd', perimeter.divide(4).log().multiply(2).divide(area.log()))
            .set('conv', convex_perim.divide(perimeter))
            .set('sqp', ee.Number(1).subtract(area.sqrt().multiply(4).divide(perimeter)))
    return shaped;
  });
  return water_feature;
}
//END OF FUNCTIONS



//SENTINEL-2 1C TOA
//Find images in research area in 2020 with less than 10% cloud pixels - 'CLOUDY_PIXEL_PERCENTAGE'
var getSentinel = ee.ImageCollection('COPERNICUS/S2').filterBounds(edited_region)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10));
//Define bands of water indices
var s2bands = {'green': 'B3', 'red': 'B4', 'nir': 'B8', 'swir1': 'B11', 'swir2': 'B12'};

//Pixel size
var scale = 10;  

//Individual districts of study region
var barisal = edited_region.filter(ee.Filter.eq('ADM2_EN', 'Barisal')).geometry();
var bhola = edited_region.filter(ee.Filter.eq('ADM2_EN', 'Bhola')).geometry();
var gopalganj = edited_region.filter(ee.Filter.eq('ADM2_EN', 'Gopalganj')).geometry();
var khulna = edited_region.filter(ee.Filter.eq('ADM2_EN', 'Khulna')).geometry();
var satkhira = edited_region.filter(ee.Filter.eq('ADM2_EN', 'Satkhira')).geometry();
var bagerhat = edited_region.filter(ee.Filter.eq('ADM2_EN', 'Bagerhat')).geometry();
var jessore = edited_region.filter(ee.Filter.eq('ADM2_EN', 'Jessore')).geometry();

//Each district indexed for loops
var districts = [barisal, bhola, gopalganj, khulna, satkhira, bagerhat, jessore];

//Single day for each district indexed for loops
var districtDayStart = ['2020-10-13', '2020-10-28', '2020-11-7', '2020-10-28', '2020-11-5', '2020-11-7', '2020-11-5'];
var districtDayEnd = ['2020-10-14', '2020-10-29', '2020-11-8', '2020-10-29', '2020-11-6', '2020-11-8', '2020-11-6'];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Creating MNDWI Mask
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
var mndwiMask = Array(7);

for (var i=0; i<7; i++) {
  var area_i = districts[i];
  var start = districtDayStart[i];
  var end = districtDayEnd[i];

  var dataIntermit = getSentinel.filterDate(start.toString(), end.toString()).filterBounds(area_i);

  var mndwi_i = getIndexSeq(dataIntermit, 'mndwi', [s2bands['green'], s2bands['swir1']]).clip(area_i);
  print('mndwi_i to Bands:', mndwi_i)
  var naive_mask = mndwi_i.gt(0).reduce(ee.Reducer.allNonZero());
  //Working with a binary image of 0s and 1s
  var vector_reduce = naive_mask.selfMask().reduceToVectors({
    geometry: area_i,
    scale: 10,
    eightConnected: true,
    maxPixels: 1e13
  });

  var add_buffer = vector_reduce.map(function (layer) {
    return layer.buffer(5*scale); 
  });

  var initial_mask = add_buffer.reduceToImage({
    properties: ['count'],
    reducer: ee.Reducer.anyNonZero()
  }).clip(area_i);

  mndwiMask[i]=initial_mask; 
}

print(mndwiMask)

///////////////////////////////////////////////////////////////////////////////////////////
//Otsu Segmentation & Edge Detection
///////////////////////////////////////////////////////////////////////////////////////////
var ndwi = Array(7);
var mndwi = Array(7);
var awei = Array(7);


for (var h=0; h<7; h++) {
  var area_h = districts[h];
  
  var dataIntermit_h = dataIntermit.filterBounds(area_h);
  
  var mndwi_h = getIndexSeq(getSentinel, 'mndwi', [s2bands['green'], s2bands['swir1']]).clip(area_h);
  var ndwi_h = getIndexSeq(getSentinel, 'ndwi', [s2bands['green'], s2bands['nir']]).clip(area_h);
  var awei_h = dataIntermit_h.map(function(img){
    return toAWEInsh(img).rename(ee.String(img.get('GENERATION_TIME')).cat('_awei'))
  }).toBands().clip(area_h);

  //If no access to ArcGIS Pro, the following is Gaussian smoothing followed by Laplacian 8x8
  //This replicated the Laplacian 5x5 in ArcGIS Pro
  
  var gaussianSmoothing = ee.Kernel.gaussian({radius: 3});
  var mndwi_smooth = mndwi_h.convolve(gaussianSmoothing);
  var ndwi_smooth = ndwi_h.convolve(gaussianSmoothing);
  var awei_smooth = awei_h.convolve(gaussianSmoothing);
  
  var laplacian = ee.Kernel.laplacian8({ normalize: false });
  var mndwi_lap = mndwi_smooth.convolve(laplacian);
  var ndwi_lap = ndwi_smooth.convolve(laplacian);
  var awei_lap = awei_smooth.convolve(laplacian);

  var ndwi_mask_h = otsu_mask(ndwi_lap, area_h, mndwiMask[h], 500, scale, true);
  var mndwi_mask_h = otsu_mask(mndwi_lap, area_h, mndwiMask[h], 500, scale, true);
  var awei_mask_h = otsu_mask(awei_lap, area_h, mndwiMask[h], 500, scale, true);

  ndwi[h] = ndwi_mask_h;
  mndwi[h] = mndwi_mask_h;
  awei[h] = awei_mask_h;

  print("NDWI mask", ndwi_mask_h)
  print("MNDWI mask", mndwi_mask_h)
  print("AWEI mask", awei_mask_h)
  Map.addLayer(ndwi_mask_h)
  Map.addLayer(mndwi_mask_h)
  Map.addLayer(awei_mask_h)
}


//IF combining more than one index together
// var barisal_mask = ee.ImageCollection.fromImages([awei0, mndwi0, ndwi0])
//             .toBands().reduce(ee.Reducer.mode());
            
// var bhola_mask = ee.ImageCollection.fromImages([awei1, mndwi1, ndwi1])
//             .toBands().reduce(ee.Reducer.mode());
            
// var gopalganj_mask = ee.ImageCollection.fromImages([awei2, mndwi2, ndwi2])
//             .toBands().reduce(ee.Reducer.mode());

// var khulna_mask = ee.ImageCollection.fromImages([awei3, mndwi3, ndwi3])
//             .toBands().reduce(ee.Reducer.mode());

// var satkhira_mask = ee.ImageCollection.fromImages([awei4, mndwi4, ndwi4])
//             .toBands().reduce(ee.Reducer.mode());

// var bagerhat_mask = ee.ImageCollection.fromImages([awei5, mndwi5, ndwi5])
//             .toBands().reduce(ee.Reducer.mode());
            
// var jessore_mask = ee.ImageCollection.fromImages([awei6, mndwi6, ndwi6])
//             .toBands().reduce(ee.Reducer.mode());
            

//Convert index rasters to polygons using best individual index per district
var waterVector0 = ndwi[0].selfMask().reduceToVectors({
  geometry:area_h,
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});
  
var waterVector1 = ndwi[1].selfMask().reduceToVectors({
  geometry: bhola,
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});

var waterVector2 = awei[2].selfMask().reduceToVectors({
  geometry: gopalganj,
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});

var waterVector3 = awei[3].selfMask().reduceToVectors({
  geometry: khulna,
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});

var waterVector4 = ndwi[4].selfMask().reduceToVectors({
  geometry: satkhira,
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});

var waterVector5 = awei[5].selfMask().reduceToVectors({
  geometry: bagerhat,
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});

var waterVector6 = awei[6].selfMask().reduceToVectors({
  geometry: jessore,
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});


// add spatial attributes to each feature
var waterShaped0 = add_spatial_attr(waterVector0, barisal);

var waterShaped1 = add_spatial_attr(waterVector1, bhola);

var waterShaped2 = add_spatial_attr(waterVector2, gopalganj);

var waterShaped3 = add_spatial_attr(waterVector3, khulna);

var waterShaped4 = add_spatial_attr(waterVector4, satkhira);

var waterShaped5 = add_spatial_attr(waterVector5, bagerhat);

var waterShaped6 = add_spatial_attr(waterVector6, jessore);

//Merge all district vectors together for final analysis
var waterShaped = waterShaped0.merge(waterShaped1).merge(waterShaped2).merge(waterShaped3)
  .merge(waterShaped4).merge(waterShaped5).merge(waterShaped6);


//Filter Ground Truthing data for only polygons
var gtFilter = gtMerge.map(function(f) {
  return f.set('geo_type', f.geometry().type())
});
var inter = gtFilter.filter(ee.Filter.eq('geo_type', 'Polygon'));

//Calculate spatial attributes for ground truthing data
var gtSpatial = add_spatial_attr(inter, edited_region);


//Creating training and validation data to cycle through with ground truthing data
//Divide randomly into 80% training and 20% validation
var dividers = [0, 0.2, 0.4, 0.6, 0.8, 1];
var groups = [];
var trainGroup = []

var sample = gtSpatial.randomColumn({
  seed: 3});
  
for(var i = 1; i < dividers.length; i++){
  var lower = dividers[i-1];
  var upper = dividers[i];
  var newGroup = sample.filter(ee.Filter.gte('random', lower)).filter(ee.Filter.lt('random', upper));
  groups.push(newGroup)
  var newGroupTrain = sample.filter(ee.Filter.lt('random', lower))
  var newGroupTrain2 = sample.filter(ee.Filter.gte('random', upper));
  var newGroupMerge = newGroupTrain.merge(newGroupTrain2)
trainGroup.push(newGroupMerge)
}

//Machine Learning Runs
//Support Vector Machine

for(var l=0; l<5; h++) {
  var trainingData = trainGroup[l]
  var validationData = groups[l]

  var bands = ['conv', 'ipq', 'pfd', 'soli', 'sqp'];
  
  //SVM - RBF kernel
  var classifier = ee.Classifier.libsvm({
  kernelType: 'RBF',
  gamma: 15,
  cost: 10
  });

  //Create training 
  var trained =classifier.train(trainingData, 'Fishpond', bands);
  
  // Classify the polygons
  var classified = waterShaped.classify(trained);
 
  //Calculate training accuracy
  var trainingClassify = trainingData.classify(trained);
  var trainAccuracy = trainingClassify.errorMatrix('Fishpond', 'classification');
  print('SVM Resubstitution error matrix: ', trainAccuracy);
  print('SVM Training overall accuracy: ', trainAccuracy.accuracy());
  
  // Classify the validation data.
  var validated = validationData.classify(trained);
  print("validated", validated)
  
  // Get a confusion matrix representing expected accuracy.
  var testAccuracy = validated.errorMatrix('Fishpond', 'classification');
  print('SVM Validation error matrix: ', testAccuracy);
  print('SVM Validation overall accuracy: ', testAccuracy.accuracy());

}


//CART
for(var j=0; j<5; j++) {
  var trainingData = trainGroup[j]
  var validationData = groups[j]

  //Define bands as spatial attributes
  var bands = ['conv', 'ipq', 'pfd', 'soli', 'sqp'];
  
  //Create training of smile Cart
  var trained = ee.Classifier.smileCart().train(trainingData, 'Fishpond', bands);
  
  // Classify the polygons
  var classified = waterShaped.classify(trained);
  
  //Calculate training accuracy
  var trainingClassify = trainingData.classify(trained);
  var trainAccuracy = trainingClassify.errorMatrix('Fishpond', 'classification');
  print('DT Resubstitution error matrix: ', trainAccuracy);
  print('DT Training overall accuracy: ', trainAccuracy.accuracy());

  // Classify the validation data.
  var validated = validationData.classify(trained);
  print("validated", validated)
  
  
  // Get a confusion matrix representing expected accuracy.
  var testAccuracy = validated.errorMatrix('Fishpond', 'classification');
  print('DT Validation error matrix: ', testAccuracy);
  print('DT Validation overall accuracy: ', testAccuracy.accuracy());

}



//Random Forest
for(var k=0; k<5; k++) {
  var trainingData = trainGroup[k]
  var validationData = groups[k]

  
  //Define bands as spatial attributes
  var bands = ['conv', 'ipq', 'pfd', 'soli', 'sqp'];
  
  var classifier = ee.Classifier.smileRandomForest({
  numberOfTrees: 1000,
  });
  
  //Create training of RF
  var trained = classifier.train(trainingData, 'Fishpond', bands);
  
  // Classify the polygons
  var classified = waterShaped.classify(trained);

  //Calculate training accuracy
  var trainingClassify = trainingData.classify(trained);
  var trainAccuracy = trainingClassify.errorMatrix('Fishpond', 'classification');
  print('RF Resubstitution error matrix: ', trainAccuracy);
  print('RF Training overall accuracy: ', trainAccuracy.accuracy());
  
  
  // Classify the validation data.
  var validated = validationData.classify(trained);
  print("validated", validated)
  
  // Get a confusion matrix representing expected accuracy.
  var testAccuracy = validated.errorMatrix('Fishpond', 'classification');
  print('RF Validation error matrix: ', testAccuracy);
  print('RF Validation overall accuracy: ', testAccuracy.accuracy());

}

