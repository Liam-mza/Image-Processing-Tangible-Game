class ImageProcessing extends PApplet {

color [] col;
boolean computeColor =true;

PImage img;
PVector tilt;

void settings() {
  size(960*6/10, 540*6/10, P3D);
}

void setup() {
  tilt = new PVector(0, 0, 0);
}

void draw() {
  if (cam.available() == true) {
        cam.read();
  }
  img = cam.get();

  img = thresholdHSB(img,92,137,21,255,58,247); //filtre HUE/BRIGHTNESS/SATURATION
  img = convolute(img); // etape "blurring"
  img = findConnectedComponents(img, true); //Detection des blob
  img = scharr(img); // etape edge detection
  img = thresholdBinary (img,185);  //elimination des pixels with low intensity
  List<PVector> lines = hough(img, nLines, regionSize); //recuperation des lignes caractéristiques des edges detected
  List<PVector> quad = detect_board(lines, img);
  tilt = getBoardTilt(quad, img);
  scale(0.6);
  //image(img, 0, 0);
  drawTheLines(cam.get(), lines, quad);
}

PImage thresholdBinary(PImage img, int threshold){
// create a new, initially transparent, 'result' image 
img.loadPixels();
PImage result = createImage(img.width, img.height, RGB); 
result.loadPixels();
for(int i = 0; i < img.width * img.height; i++) {
  float br =brightness(img.pixels[i]);
  if (br< threshold) {
  result.pixels[i]=color(0,0,0);
  }
  else {result.pixels[i]=color(255,255,255);}
}
return result;
}

PImage thresholdBinaryInv(PImage img, int threshold){
// create a new, initially transparent, 'result' image 
img.loadPixels();
PImage result = createImage(img.width, img.height, RGB); 
result.loadPixels();
for(int i = 0; i < img.width * img.height; i++) {
  float br =brightness(img.pixels[i]);
  if (br< threshold) {
  result.pixels[i]=color(255,255,255);
  }
  else {result.pixels[i]=color(0,0,0);}
}
return result;
}

PImage thresholdHSB(PImage img, int minH, int maxH, int minS, int maxS, int minB, int maxB) {
  PImage result = createImage(img.width, img.height, RGB); 
  img.loadPixels();
  result.loadPixels();
  for(int i = 0; i < img.width * img.height; i++) {
    float huePix =(hue(img.pixels[i]));
    float SatPix =(saturation(img.pixels[i]));
    float BrigPix =(brightness(img.pixels[i]));
    if (huePix>=minH && huePix<=maxH && SatPix >= minS && SatPix <= maxS && BrigPix >= minB && BrigPix <= maxB)
    result.pixels[i]= color (255,255,255);
    else result.pixels[i]=color(0,0,0);
  }
  return result;
}


boolean imagesEqual(PImage img1, PImage img2){
if(img1.width != img2.width || img1.height != img2.height)
return false;
for(int i = 0; i < img1.width*img1.height ; i++)
//assuming that all the three channels have the same value
  if(red(img1.pixels[i]) != red(img2.pixels[i])) {
    println(i);
    println(red(img1.pixels[i])+" et "+red(img2.pixels[i]));
    return false;}
return true; 
}

// methode de convolution avec flitre: Kernel3 permet de réaliser un blurring
PImage convolute(PImage img) {
float[][] kernel = { { 0, 0, 0 }, 
                     { 0, 2, 0 },
                     { 0, 0, 0 }};

float[][] kernel2 = { { 0, 1, 0 }, 
                      { 1, 0, 1 },
                      { 0, 1, 0 }};
                   
float[][] kernel3 = { { 9, 12, 9 }, 
                      { 12, 15, 12 },
                      { 9, 12, 9 }};                  
                     
float normFactor = 99.f;
// create a greyscale image (type: ALPHA) for output
img.loadPixels();
PImage result = createImage(img.width, img.height, ALPHA);

for(int x=0;x<img.width;x++){
  for(int y=0;y<img.height;y++){ 
    int sum=0; 
    for(int a=-1;a<2;a++){
      for(int b=-1;b<2;b++) {
        sum+=kernel3[a+1][b+1]*brightness(img.pixels[constrain(y-a,0,img.height-1) * img.width + constrain(x-b,0,img.width-1)]);
      } 
    }
    
    result.pixels[y * img.width + x]= color(sum/normFactor);
    if(y==0||x==0||x==img.width-1||y==img.height-1){
    result.pixels[y * img.width + x]=color(0);}
  }
}
return result; }


PImage scharr(PImage img) {
float[][] vKernel = { { 3,0,-3 },
                    { 10, 0, -10 },
                      { 3,0,-3 }};

float[][] hKernel = { { 3,10,3},
                      { 0,0,0},
                   { -3, -10, -3 } };

img.loadPixels();

PImage result = createImage(img.width, img.height, ALPHA);
        // clear the image
for (int i = 0; i < img.width * img.height; i++) {
          result.pixels[i] = color(0);
}
        
float max=0;
float normFactor = 1.f;
float[] buffer = new float[img.width * img.height];
// *************************************
// Implement here the double convolution

for(int x=0;x<img.width;x++){
  for(int y=0;y<img.height;y++){ 
    int sum=0; 
    int sum2=0; 
    for(int a=-1;a<2;a++){
      for(int b=-1;b<2;b++) {
        sum +=vKernel[a+1][b+1]*brightness(img.pixels[constrain(y-a,0,img.height-1) * img.width + constrain(x-b,0,img.width-1)]);
        sum2 +=hKernel[a+1][b+1]*brightness(img.pixels[constrain(y-a,0,img.height-1) * img.width + constrain(x-b,0,img.width-1)]);
      } 
    }
   float totVal = sqrt(pow(sum2, 2) + pow(sum, 2));
   buffer[y * img.width + x]= totVal;
   if (max<totVal){max=totVal;}
  }
}

// *************************************
for (int y = 1; y < img.height - 1; y++) { // Skip top and bottom edges 
  for (int x = 1; x < img.width - 1; x++) { // Skip left and right
int val=(int) ((buffer[y * img.width + x] / max)*255);
result.pixels[y * img.width + x]=color(val);
  }
}
return result; }
    
//Etape de détection des blobs
PImage findConnectedComponents(PImage input, boolean onlyBiggest) {
  
// First pass: label the pixels and store labels' equivalences
int [] labels = new int [input.width*input.height];
PImage result = createImage(input.width, input.height, ALPHA);
List<TreeSet<Integer>> labelsEquivalences = new ArrayList<TreeSet<Integer>>();
int currentLabel = 0; 


for (int i=0; i<input.width*input.height;++i){
  labels[i]=-1;
}
input.loadPixels();
// TODO!


for(int y=0;y<img.height;y++){
  for(int x=0;x<img.width;x++){ 
    if (input.pixels[y * img.width + x]!=color(0,0,0)){
    TreeSet<Integer> tmpEq = new TreeSet<Integer>();
    if (labels[y * img.width + constrain(x-1,0,img.width-1)]!=-1){
      tmpEq.add(labels[y * img.width + constrain(x-1,0,img.width-1)]);
    }
    for (int a=-1;a<2;++a){
      if (labels[constrain(y-1,0,img.height-1) * img.width + constrain(x+a,0,img.width-1)]!=-1){
      tmpEq.add(labels[constrain(y-1,0,img.height-1) * img.width + constrain(x+a,0,img.width-1)]);
      }
    }
    if (tmpEq.isEmpty()){
      labels[y * img.width + x]= currentLabel;
      tmpEq.add(currentLabel);
      labelsEquivalences.add(tmpEq);
      ++currentLabel;
    }
    else{
      int v= tmpEq.first();
      labels[y * img.width + x]= v;
      labelsEquivalences.get(v).addAll(tmpEq);
      for ( Integer e: tmpEq){
        labelsEquivalences.get(e).addAll(labelsEquivalences.get(v));
      }
    }
  }
  }
}
// Second pass: re-label the pixels by their equivalent class
// if onlyBiggest==true, count the number of pixels for each label
Integer [] count = new Integer[currentLabel];
for (int i=0; i<input.width*input.height;++i){
  if(labels[i]!=-1){
    labels[i]= labelsEquivalences.get(labels[i]).first();
    if (count[labels[i]]==null){count[labels[i]]=0;}
    else{++count[labels[i]];}
    
  }
}


// TODO!
// Finally:
// if onlyBiggest==false, output an image with each blob colored in one uniform color
// if onlyBiggest==true, output an image with the biggest blob in white and others in black

if (onlyBiggest==false){
  if(computeColor){
    col =new color[currentLabel];
    for(int i=0; i<currentLabel;++i){
       col[i]=color(random(0,255),random(0,255),random(0,255));
    }
    computeColor=false;
  }
  for (int i=0; i<input.width*input.height;++i){
    if (labels[i]!=-1){
    result.pixels[i]=col[labels[i]];
    
    }
    else {
      result.pixels[i]=color(0,0,0);
    }
  }
}
else {
  int max=0;
  int maxV=0;
   for (int i=0; i<count.length;++i){
    if (count[i]!=null &&count[i]>maxV){
    maxV =count[i];
    max=i;
  }
  }
  for (int i=0; i<input.width*input.height;++i){
    if (labels[i]== max){
      result.pixels[i]=color(255,255,255);
    }
    else {
      result.pixels[i]=color(0,0,0);
    }
  }
}

// TODO!
return result;
}



OpenCV opencv;

//PARAMETERS TO TUNE
boolean drawLines = true;
float discretizationStepsPhi = 0.07f; 
float discretizationStepsR = 1.7f;
int minVotes = 100;
int nLines = 10;
int regionSize = 100; //Definition of "local" for the local maxima
int max_quad_area = 2147483647;
int min_quad_area = 0;
int fps = 0; //If we manage to get the camera working

List<PVector> detect_board(List<PVector> lines, PImage img) {
  QuadGraph quadGraph = new QuadGraph();
  List<PVector> bestQuad = quadGraph.findBestQuad(lines, img.width, img.height, max_quad_area, min_quad_area, false);

  return bestQuad;
}

PVector getBoardTilt(List<PVector> bestQuad, PImage hough_test) {
  opencv = new OpenCV(this, 100, 100);
  TwoDThreeD twoToThree = new TwoDThreeD(hough_test.width, hough_test.height, fps);
  List<PVector> homogenQuad = bestQuad;
  for(PVector p : homogenQuad) {
    p = p.add(new PVector(0, 0, 1)); //Convert to homogeneous coordinates
  }
  PVector rotations = twoToThree.get3DRotations(homogenQuad); //Angles in radians
  //println(rotations.mult(180).div(PI)); //Degrees -- May have to add or subtract 180 to be close to 0 - +-180 is the same line
  return rotations;
}


//HOUGH_TRANSFORM ALGORITHM
List<PVector> hough(PImage edgeImg, int nLines, int regionSize) {
  
  // dimensions of the accumulator
  int phiDim = (int)(Math.PI / discretizationStepsPhi +1);
  //The max radius is the image diagonal, but it can be also negative 
  int rDim = (int)((sqrt(edgeImg.width*edgeImg.width + edgeImg.height*edgeImg.height) * 2) / discretizationStepsR +1);
  
  // our accumulator
  int[] accumulator = new int[phiDim * rDim];
  
  // pre-computed the sin and cos values
  float[] tabSin = new float[phiDim];
  float[] tabCos = new float[phiDim];

  float ang = 0;
  float inverseR = 1.f / discretizationStepsR;
  
  for (int accPhi = 0; accPhi < phiDim; ang += discretizationStepsPhi, accPhi++) {
    // we can also pre-multiply by (1/discretizationStepsR) since we need it in the Hough loop 
    tabSin[accPhi] = (float) (Math.sin(ang) * inverseR);
    tabCos[accPhi] = (float) (Math.cos(ang) * inverseR);
  }
  
  // Fill the accumulator: on edge points (ie, white pixels of the edge 
  // image), store all possible (r, phi) pairs describing lines going 
  // through the point.
  for (int y = 0; y < edgeImg.height; y++) {
    for (int x = 0; x < edgeImg.width; x++) {
      // Are we on an edge?
      if (brightness(edgeImg.pixels[y * edgeImg.width + x]) != 0) {
        // ...determine here all the lines (r, phi) passing through
        // pixel (x,y), convert (r,phi) to coordinates in the
        // accumulator, and increment accordingly the accumulator.
        // Be careful: r may be negative, so you may want to center onto
        // the accumulator: r += rDim / 2
        for(int phi=0; phi<phiDim; phi += 1) {
          int r = (int)(x*tabCos[phi] + y*tabSin[phi]);
          r += rDim / 2;
          accumulator[phi * rDim + r] += 1;
        }
      }
    }
  }
  
  if(!drawLines) {
    PImage houghImg = createImage(rDim, phiDim, ALPHA);
    houghImg.loadPixels();
    for (int i = 0; i < accumulator.length; i++) {
          houghImg.pixels[i] = color(min(255, accumulator[i]));
    }
    houghImg.updatePixels();
    houghImg.resize(400, 400);
    image(houghImg, 0, 0);
  }
  
  ArrayList<PVector> lines = new ArrayList<PVector>();
  ArrayList<Integer> bestCandidates = new ArrayList<Integer>();
  for (int idx = 0; idx < accumulator.length; idx++) {
    if (accumulator[idx] > minVotes) {
      int currLineCount = accumulator[idx];
      boolean maxima = true;
      for(int i=max(idx-regionSize, 0); i<min(idx+regionSize, accumulator.length); ++i) {
        if(accumulator[i] > currLineCount) {
          maxima = false;
        }
      }
      if(maxima) {
        bestCandidates.add(idx);
      }
    }
  }
  Collections.sort(bestCandidates, new HoughComparator(accumulator));
  for(int i=0; i<min(bestCandidates.size(), nLines); ++i) {
    int accPhi = (int) (bestCandidates.get(i) / (rDim));
    int accR = bestCandidates.get(i) - (accPhi) * (rDim);
    float r = (accR - (rDim) * 0.5f) * discretizationStepsR; 
    float phi = accPhi * discretizationStepsPhi; 
    lines.add(new PVector(r,phi));
  }
  return lines;
}




void drawTheLines(PImage hough_test, List<PVector> lines, List<PVector> bestQuad) {
    image(hough_test, 0, 0);
    println("UGDE");
    //DRAW ALL DETECTED LINES FROM THE HOUGH ALGORITHM
    for (int idx = 0; idx < lines.size(); idx++) {
      PVector line= lines.get(idx);
      float r = line.x; 
      float phi = line.y;
      // Cartesian equation of a line: y = ax + b
      // in polar, y = (-cos(phi)/sin(phi))x + (r/sin(phi))
      // => y = 0 : x = r / cos(phi)
      // => x = 0 : y = r / sin(phi)
      // compute the intersection of this line with the 4 borders of // the image
      int x0 = 0;
      int y0 = (int) (r / sin(phi));
      int x1 = (int) (r / cos(phi));
      int y1 = 0;
      int x2 = hough_test.width;
      int y2 = (int) (-cos(phi) / sin(phi) * x2 + r / sin(phi)); 
      int y3 = hough_test.width;
      int x3 = (int) (-(y3 - r / sin(phi)) * (sin(phi) / cos(phi)));
      
      // Finally, plot the lines
      stroke(204,102,0); 
      if (y0 > 0) {
        if (x1 > 0) {
          line(x0, y0, x1, y1);
        }
        else if (y2 > 0) {
          line(x0, y0, x2, y2);
        }
        else {
          line(x0, y0, x3, y3);
        }
      } 
      else {
        if (x1 > 0) {
          if (y2 > 0) {
            line(x1, y1, x2, y2); 
          }
          else {
            line(x1, y1, x3, y3);
          }
        } 
        else {
          line(x2, y2, x3, y3);
        }
      }
    }
    
    //DRAW THE QUAD
    for(PVector l : bestQuad) {
      fill(102, 203, 0);
      ellipse(l.x, l.y, 30, 30);
    }
}

class HoughComparator implements java.util.Comparator<Integer> { 
  
  int[] accumulator;
  
  public HoughComparator(int[] accumulator) {
    this.accumulator = accumulator; 
  }
  
  @Override
  public int compare(Integer l1, Integer l2) { 
    if (accumulator[l1] > accumulator[l2] || (accumulator[l1] == accumulator[l2] && l1 < l2)) {
      return -1; 
    }
    return 1;
  }
  
}

class QuadGraph {

  boolean verbose=false;
  
  List<int[]> cycles = new ArrayList<int[]>();
  int[][] graph;

  List<PVector> findBestQuad(List<PVector> lines, int width, int height, int max_quad_area, int min_quad_area, boolean verbose) {
    this.verbose=verbose;
    build(lines, width, height); //<>//
    findCycles(verbose);
    ArrayList<PVector> bestQuad=new ArrayList<PVector>();
    float bestQuadArea=0;
    for (int [] cy : cycles) {
      ArrayList<PVector> quad= new ArrayList<PVector>();
      PVector l1 = lines.get(cy[0]);
      PVector l2 = lines.get(cy[1]);
      PVector l3 = lines.get(cy[2]);
      PVector l4 = lines.get(cy[3]);


      quad.add(intersection(l1, l2));
      quad.add(intersection(l2, l3));
      quad.add(intersection(l3, l4));
      quad.add(intersection(l4, l1));
      quad=sortCorners(quad);

      PVector c1 = quad.get(0);
      PVector c2 = quad.get(1);
      PVector c3 = quad.get(2);
      PVector c4 = quad.get(3);

      if (isConvex(c1, c2, c3, c4) && 
        nonFlatQuad(c1, c2, c3, c4)) {
        float quadArea=validArea(c1, c2, c3, c4, max_quad_area, min_quad_area);
        if (quadArea>0 && quadArea>bestQuadArea) {
          bestQuadArea=quadArea;
          bestQuad=quad;
        }
      }
    }
    if (bestQuadArea>0)
      return bestQuad;
    else
      return new ArrayList<PVector>();
  }


  void build(List<PVector> lines, int width, int height) {

    int n = lines.size();

    // The maximum possible number of edges is n * (n - 1)/2
    graph = new int[n * (n - 1)/2][2];

    int idx =0;

    for (int i = 0; i < lines.size(); i++) {
      for (int j = i + 1; j < lines.size(); j++) {
        if (intersect(lines.get(i), lines.get(j), width, height)) {

          graph[idx][0]=i;
          graph[idx][1]=j;
          idx++;
        }
      }
    }
  }

  /** Returns true if polar lines 1 and 2 intersect 
   * inside an area of size (width, height)
   */
  boolean intersect(PVector line1, PVector line2, int width, int height) {

    double sin_t1 = Math.sin(line1.y);
    double sin_t2 = Math.sin(line2.y);
    double cos_t1 = Math.cos(line1.y);
    double cos_t2 = Math.cos(line2.y);
    float r1 = line1.x;
    float r2 = line2.x;

    double denom = cos_t2 * sin_t1 - cos_t1 * sin_t2;

    int x = (int) ((r2 * sin_t1 - r1 * sin_t2) / denom);
    int y = (int) ((-r2 * cos_t1 + r1 * cos_t2) / denom);

    if (0 <= x && 0 <= y && width >= x && height >= y)
      return true;
    else
      return false;
  }
  
  PVector intersection(PVector line1, PVector line2) {

    double sin_t1 = Math.sin(line1.y);
    double sin_t2 = Math.sin(line2.y);
    double cos_t1 = Math.cos(line1.y);
    double cos_t2 = Math.cos(line2.y);
    float r1 = line1.x;
    float r2 = line2.x;

    double denom = cos_t2 * sin_t1 - cos_t1 * sin_t2;

    int x = (int) ((r2 * sin_t1 - r1 * sin_t2) / denom);
    int y = (int) ((-r2 * cos_t1 + r1 * cos_t2) / denom);

    return new PVector(x,y);
  }

  void findCycles(boolean verbose) {
    cycles.clear();
    for (int i = 0; i < graph.length; i++) {
      for (int j = 0; j < graph[i].length; j++) {
        findNewCycles(new int[] {graph[i][j]});
      }
    }
    if (verbose) {
      for (int[] cy : cycles) {
        String s = "" + cy[0];
        for (int i = 1; i < cy.length; i++) {
          s += "," + cy[i];
        }
        System.out.println(s);
      }
    }
  }

  void findNewCycles(int[] path)
  {
    int n = path[0];
    int x;
    int[] sub = new int[path.length + 1];

    for (int i = 0; i < graph.length; i++)
      for (int y = 0; y <= 1; y++)
        if (graph[i][y] == n)
          //  edge refers to our current node
        {
          x = graph[i][(y + 1) % 2];
          if (!visited(x, path))
            //  neighbor node not on path yet
          {
            sub[0] = x;
            System.arraycopy(path, 0, sub, 1, path.length);
            //  explore extended path
            findNewCycles(sub);
          } else if ((path.length == 4) && (x == path[path.length - 1]))
            //  cycle found
          {
            int[] p = normalize(path);
            int[] inv = invert(p);
            if (isNew(p) && isNew(inv))
            {
              cycles.add(p);
            }
          }
        }
  }

  //  check of both arrays have same lengths and contents
  Boolean equals(int[] a, int[] b)
  {
    Boolean ret = (a[0] == b[0]) && (a.length == b.length);

    for (int i = 1; ret && (i < a.length); i++)
    {
      if (a[i] != b[i])
      {
        ret = false;
      }
    }

    return ret;
  }

  //  create a path array with reversed order
  int[] invert(int[] path)
  {
    int[] p = new int[path.length];

    for (int i = 0; i < path.length; i++)
    {
      p[i] = path[path.length - 1 - i];
    }

    return normalize(p);
  }

  //  rotate cycle path such that it begins with the smallest node
  int[] normalize(int[] path)
  {
    int[] p = new int[path.length];
    int x = smallest(path);
    int n;

    System.arraycopy(path, 0, p, 0, path.length);

    while (p[0] != x)
    {
      n = p[0];
      System.arraycopy(p, 1, p, 0, p.length - 1);
      p[p.length - 1] = n;
    }

    return p;
  }

  //  compare path against known cycles
  //  return true, iff path is not a known cycle
  Boolean isNew(int[] path)
  {
    Boolean ret = true;

    for (int[] p : cycles)
    {
      if (equals(p, path))
      {
        ret = false;
        break;
      }
    }

    return ret;
  }

  //  return the int of the array which is the smallest
  int smallest(int[] path)
  {
    int min = path[0];

    for (int p : path)
    {
      if (p < min)
      {
        min = p;
      }
    }

    return min;
  }

  //  check if vertex n is contained in path
  Boolean visited(int n, int[] path)
  {
    Boolean ret = false;

    for (int p : path)
    {
      if (p == n)
      {
        ret = true;
        break;
      }
    }

    return ret;
  }



  /** Check if a quad is convex or not.
   * 
   * Algo: take two adjacent edges and compute their cross-product. 
   * The sign of the z-component of all the cross-products is the 
   * same for a convex polygon.
   * 
   * See http://debian.fmi.uni-sofia.bg/~sergei/cgsr/docs/clockwise.htm
   * for justification.
   * 
   * @param c1
   */
  boolean isConvex(PVector c1, PVector c2, PVector c3, PVector c4) {

    PVector v21= PVector.sub(c1, c2);
    PVector v32= PVector.sub(c2, c3);
    PVector v43= PVector.sub(c3, c4);
    PVector v14= PVector.sub(c4, c1);

    float i1=v21.cross(v32).z;
    float i2=v32.cross(v43).z;
    float i3=v43.cross(v14).z;
    float i4=v14.cross(v21).z;

    if (   (i1>0 && i2>0 && i3>0 && i4>0) 
      || (i1<0 && i2<0 && i3<0 && i4<0))
      return true;
    else if(verbose)
      System.out.println("Eliminating non-convex quad");
    return false;
  }

  /** Compute the area of a quad, and check it lays within a specific range
   */
  float validArea(PVector c1, PVector c2, PVector c3, PVector c4, float max_area, float min_area) {

    float i1=c1.cross(c2).z;
    float i2=c2.cross(c3).z;
    float i3=c3.cross(c4).z;
    float i4=c4.cross(c1).z;

    float area = Math.abs(0.5f * (i1 + i2 + i3 + i4));

    
    if (area < max_area && area > min_area){
      return area;
    }
    return 0;
    
   }

  /** Compute the (cosine) of the four angles of the quad, and check they are all large enough
   * (the quad representing our board should be close to a rectangle)
   */
  boolean nonFlatQuad(PVector c1, PVector c2, PVector c3, PVector c4) {

    // cos(70deg) ~= 0.3
    float min_cos = 0.5f;

    PVector v21= PVector.sub(c1, c2);
    PVector v32= PVector.sub(c2, c3);
    PVector v43= PVector.sub(c3, c4);
    PVector v14= PVector.sub(c4, c1);

    float cos1=Math.abs(v21.dot(v32) / (v21.mag() * v32.mag()));
    float cos2=Math.abs(v32.dot(v43) / (v32.mag() * v43.mag()));
    float cos3=Math.abs(v43.dot(v14) / (v43.mag() * v14.mag()));
    float cos4=Math.abs(v14.dot(v21) / (v14.mag() * v21.mag()));

    if (cos1 < min_cos && cos2 < min_cos && cos3 < min_cos && cos4 < min_cos)
      return true;
    else {
      if(verbose)
        System.out.println("Flat quad");
      return false;
    }
  }


  ArrayList<PVector> sortCorners(ArrayList<PVector> quad) {

    // 1 - Sort corners so that they are ordered clockwise
    PVector a = quad.get(0);
    PVector b = quad.get(2);

    PVector center = new PVector((a.x+b.x)/2, (a.y+b.y)/2);

    Collections.sort(quad, new CWComparator(center));



    // 2 - Sort by upper left most corner
    PVector origin = new PVector(0, 0);
    float distToOrigin = 1000;

    for (PVector p : quad) {
      if (p.dist(origin) < distToOrigin) distToOrigin = p.dist(origin);
    }

    while (quad.get(0).dist(origin) != distToOrigin)
      Collections.rotate(quad, 1);

    return quad;
  }
}

class CWComparator implements Comparator<PVector> {

  PVector center;

  public CWComparator(PVector center) {
    this.center = center;
  }

  @Override
    public int compare(PVector b, PVector d) {
    if (Math.atan2(b.y-center.y, b.x-center.x)<Math.atan2(d.y-center.y, d.x-center.x))      
      return -1; 
    else return 1;
  }
}



class TwoDThreeD {

  // default focal length, well suited for most webcams
  float f = 700;

  // intrisic camera matrix
  float [][] K = {{f, 0, 0},
    {0, f, 0},
    {0, 0, 1}};
  float [][] invK;
  PVector invK_r1, invK_r2, invK_r3;
  Mat opencv_A, w, u, vt;
  double [][] V;

  // Real physical coordinates of the Lego board in mm
  //float boardSize = 380.f; // large Duplo board
  // float boardSize = 255.f; // smaller Lego board

  // the 3D coordinates of the physical board corners, clockwise
  float [][] physicalCorners = {
    {-128, 128, 0, 1},
    {128, 128, 0, 1},
    {128, -128, 0, 1},
    {-128, -128, 0, 1}
  };

  //Filtering variables: low-pass filter based on arFilterTrans from ARToolKit v5 */
  float[] q;
  float sampleRate;
  float cutOffFreq;
  float alpha;


  public TwoDThreeD(int width, int height, float sampleRate) {

    // set the offset to the center of the webcam image
    K[0][2] = 0.5f * width;
    K[1][2] = 0.5f * height;
    //compute inverse of K
    Mat opencv_K= new Mat(3, 3, CvType.CV_32F);
    opencv_K.put(0, 0, K[0][0]);
    opencv_K.put(0, 1, K[0][1]);
    opencv_K.put(0, 2, K[0][2]);
    opencv_K.put(1, 0, K[1][0]);
    opencv_K.put(1, 1, K[1][1]);
    opencv_K.put(1, 2, K[1][2]);
    opencv_K.put(2, 0, K[2][0]);
    opencv_K.put(2, 1, K[2][1]);
    opencv_K.put(2, 2, K[2][2]);
    Mat opencv_invK=opencv_K.inv();

    invK = new float[][]{
      { (float)opencv_invK.get(0, 0)[0], (float)opencv_invK.get(0, 1)[0], (float)opencv_invK.get(0, 2)[0] },
      { (float)opencv_invK.get(1, 0)[0], (float)opencv_invK.get(1, 1)[0], (float)opencv_invK.get(1, 2)[0] },
      { (float)opencv_invK.get(2, 0)[0], (float)opencv_invK.get(2, 1)[0], (float)opencv_invK.get(2, 2)[0] }};
    invK_r1=new PVector(invK[0][0], invK[0][1], invK[0][2]);
    invK_r2=new PVector(invK[1][0], invK[1][1], invK[1][2]);
    invK_r3=new PVector(invK[2][0], invK[2][1], invK[2][2]);

    opencv_A=new Mat(12, 9, CvType.CV_32F);
    w=new Mat();
    u=new Mat();
    vt=new Mat();
    V= new double[9][9];

    q=new float[4];
    q[3]=1;

    this.sampleRate=sampleRate;
    if (sampleRate>0) {
      cutOffFreq=sampleRate/2;
      alpha= (1/sampleRate)/(1/sampleRate + 1/cutOffFreq);
    }
  }

  PVector get3DRotations(List<PVector> points2D) {

    // 1- Solve the extrinsic matrix from the projected 2D points
    double[][] E = solveExtrinsicMatrix(points2D);


    // 2 - Re-build a proper 3x3 rotation matrix from the camera's
    //     extrinsic matrix E
    PVector firstColumn=new PVector((float)E[0][0], (float)E[1][0], (float)E[2][0]);
    PVector secondColumn=new PVector((float)E[0][1], (float)E[1][1], (float)E[2][1]);
    firstColumn.normalize();
    secondColumn.normalize();
    PVector thirdColumn=firstColumn.cross(secondColumn);
    float [][] rotationMatrix={{firstColumn.x, secondColumn.x, thirdColumn.x},
      {firstColumn.y, secondColumn.y, thirdColumn.y},
      {firstColumn.z, secondColumn.z, thirdColumn.z}};

    if (sampleRate>0)
      filter(rotationMatrix, false);

    // 3 - Computes and returns Euler angles (rx, ry, rz) from this matrix
    return rotationFromMatrix(rotationMatrix);
  }


  double[][] solveExtrinsicMatrix(List<PVector> points2D) {

    // p ~= K · [R|t] · P
    // with P the (3D) corners of the physical board, p the (2D)
    // projected points onto the webcam image, K the intrinsic
    // matrix and R and t the rotation and translation we want to
    // compute.
    //
    // => We want to solve: (K^(-1) · p) X ([R|t] · P) = 0

    float[][] projectedCorners = new float[4][3];

    if(points2D.size() >= 4)
    for (int i=0; i<4; i++) {
      // TODO:
      // store in projectedCorners the result of (K^(-1) · p), for each
      // corner p found in the webcam image.
      // You can use PVector dot function for computing dot product between K^(-1) lines and p.
      //Do not forget to normalize the result
      PVector point =points2D.get(i);
      projectedCorners[i][0]=point.dot(invK_r1)/point.dot(invK_r3);
      projectedCorners[i][1]=point.dot(invK_r2)/point.dot(invK_r3);
      projectedCorners[i][2]=1;
    }

    // 'A' contains the cross-product (K^(-1) · p) X P
    float[][] A= new float[12][9];

    for (int i=0; i<4; i++) {
      A[i*3][0]=0;
      A[i*3][1]=0;
      A[i*3][2]=0;

      // note that we take physicalCorners[0,1,*3*]: we drop the Z
      // coordinate and use the 2D homogenous coordinates of the physical
      // corners
      A[i*3][3]=-projectedCorners[i][2] * physicalCorners[i][0];
      A[i*3][4]=-projectedCorners[i][2] * physicalCorners[i][1];
      A[i*3][5]=-projectedCorners[i][2] * physicalCorners[i][3];

      A[i*3][6]= projectedCorners[i][1] * physicalCorners[i][0];
      A[i*3][7]= projectedCorners[i][1] * physicalCorners[i][1];
      A[i*3][8]= projectedCorners[i][1] * physicalCorners[i][3];

      A[i*3+1][0]= projectedCorners[i][2] * physicalCorners[i][0];
      A[i*3+1][1]= projectedCorners[i][2] * physicalCorners[i][1];
      A[i*3+1][2]= projectedCorners[i][2] * physicalCorners[i][3];

      A[i*3+1][3]=0;
      A[i*3+1][4]=0;
      A[i*3+1][5]=0;

      A[i*3+1][6]=-projectedCorners[i][0] * physicalCorners[i][0];
      A[i*3+1][7]=-projectedCorners[i][0] * physicalCorners[i][1];
      A[i*3+1][8]=-projectedCorners[i][0] * physicalCorners[i][3];

      A[i*3+2][0]=-projectedCorners[i][1] * physicalCorners[i][0];
      A[i*3+2][1]=-projectedCorners[i][1] * physicalCorners[i][1];
      A[i*3+2][2]=-projectedCorners[i][1] * physicalCorners[i][3];

      A[i*3+2][3]= projectedCorners[i][0] * physicalCorners[i][0];
      A[i*3+2][4]= projectedCorners[i][0] * physicalCorners[i][1];
      A[i*3+2][5]= projectedCorners[i][0] * physicalCorners[i][3];

      A[i*3+2][6]=0;
      A[i*3+2][7]=0;
      A[i*3+2][8]=0;
    }

    for (int i=0; i<12; i++)
      for (int j=0; j<9; j++)
        opencv_A.put(i, j, A[i][j]);

    Core.SVDecomp(opencv_A, w, u, vt);

    for (int i=0; i<9; i++)
      for (int j=0; j<9; j++)
        V[j][i]=vt.get(i, j)[0];

    double[][] E = new double[3][3];

    //E is the last column of V
    for (int i=0; i<9; i++) {
      E[i/3][i%3] = V[i][V.length-1] / V[8][V.length-1];
    }

    return E;
  }

  PVector rotationFromMatrix(float[][]  mat) {

    // Assuming rotation order is around x,y,z
    PVector rot = new PVector();

    if (mat[1][0] > 0.998) { // singularity at north pole
      rot.z = 0;
      float delta = (float) Math.atan2(mat[0][1], mat[0][2]);
      rot.y = -(float) Math.PI/2;
      rot.x = -rot.z + delta;
      return rot;
    }

    if (mat[1][0] < -0.998) { // singularity at south pole
      rot.z = 0;
      float delta = (float) Math.atan2(mat[0][1], mat[0][2]);
      rot.y = (float) Math.PI/2;
      rot.x = rot.z + delta;
      return rot;
    }

    rot.y =-(float)Math.asin(mat[2][0]);
    rot.x = (float)Math.atan2(mat[2][1]/Math.cos(rot.y), mat[2][2]/Math.cos(rot.y));
    rot.z = (float)Math.atan2(mat[1][0]/Math.cos(rot.y), mat[0][0]/Math.cos(rot.y));

    return rot;
  }

  int filter(float m[][], boolean reset) {

    float[] q= new float[4];
    float alpha, oneminusalpha, omega, cosomega, sinomega, s0, s1;

    mat2Quat(m, q);
    if (nomalizeQuaternion(q)<0) return -1;

    if (reset) {
      this.q[0] = q[0];
      this.q[1] = q[1];
      this.q[2] = q[2];
      this.q[3] = q[3];
    } else {
      alpha = this.alpha;

      oneminusalpha = 1.0 - alpha;

      // SLERP for orientation.
      cosomega = q[0]*this.q[0] + q[1]*this.q[1] + q[2]*this.q[2] + q[3]*this.q[3]; // cos of angle between vectors.
      if (cosomega < 0.0) {
        cosomega = -cosomega;
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
      }
      if (cosomega > 0.9995) {
        s0 = oneminusalpha;
        s1 = alpha;
      } else {
        omega = acos(cosomega);
        sinomega = sin(omega);
        s0 = sin(oneminusalpha * omega) / sinomega;
        s1 = sin(alpha * omega) / sinomega;
      }
      this.q[0] = q[0]*s1 + this.q[0]*s0;
      this.q[1] = q[1]*s1 + this.q[1]*s0;
      this.q[2] = q[2]*s1 + this.q[2]*s0;
      this.q[3] = q[3]*s1 + this.q[3]*s0;
      nomalizeQuaternion(this.q);
    }

    if (quat2Mat(this.q, m) < 0) return (-2);

    return (0);
  }


  int nomalizeQuaternion(float[] q) {// Normalise quaternion.
    float mag2 = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    if (mag2==0) return (-1);

    float mag = sqrt(mag2);

    q[0] /= mag;
    q[1] /= mag;
    q[2] /= mag;
    q[3] /= mag;

    return (0);
  }

  int mat2Quat(float m[][], float q[]) {
    float t, s;
    t = m[0][0] + m[1][1] + m[2][2] + 1.0;
    if (t > 0.0001) {
      s = sqrt(t) * 2.0;
      q[0] = (m[1][2] - m[2][1]) / s;
      q[1] = (m[2][0] - m[0][2]) / s;
      q[2] = (m[0][1] - m[1][0]) / s;
      q[3] = 0.25 * s;
    } else {
      if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {  // Column 0:
        s  = sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2.0;
        q[0] = 0.25 * s;
        q[1] = (m[0][1] + m[1][0] ) / s;
        q[2] = (m[2][0] + m[0][2] ) / s;
        q[3] = (m[1][2] - m[2][1] ) / s;
      } else if (m[1][1] > m[2][2]) {      // Column 1:
        s  = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2.0;
        q[0] = (m[0][1] + m[1][0] ) / s;
        q[1] = 0.25 * s;
        q[2] = (m[1][2] + m[2][1] ) / s;
        q[3] = (m[2][0] - m[0][2] ) / s;
      } else {            // Column 2:
        s  = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2.0;
        q[0] = (m[2][0] + m[0][2] ) / s;
        q[1] = (m[1][2] + m[2][1] ) / s;
        q[2] = 0.25 * s;
        q[3] = (m[0][1] - m[1][0] ) / s;
      }
    }
    return 0;
  }

  int quat2Mat( float q[], float m[][] )
  {
    float    x2, y2, z2;
    float    xx, xy, xz;
    float    yy, yz, zz;
    float    wx, wy, wz;

    x2 = q[0] * 2.0;
    y2 = q[1] * 2.0;
    z2 = q[2] * 2.0;

    xx = q[0] * x2;
    xy = q[0] * y2;
    xz = q[0] * z2;
    yy = q[1] * y2;
    yz = q[1] * z2;
    zz = q[2] * z2;
    wx = q[3] * x2;
    wy = q[3] * y2;
    wz = q[3] * z2;

    m[0][0] = 1.0 - (yy + zz);
    m[1][1] = 1.0 - (xx + zz);
    m[2][2] = 1.0 - (xx + yy);

    m[1][0] = xy - wz;
    m[0][1] = xy + wz;
    m[2][0] = xz + wy;
    m[0][2] = xz - wy;
    m[2][1] = yz - wx;
    m[1][2] = yz + wx;

    return 0;
  }
}

}
