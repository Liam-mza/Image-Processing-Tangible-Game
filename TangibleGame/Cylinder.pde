float cylinderBaseSize = 20;
float cylinderHeight = 40; 
int cylinderResolution = 40;
PShape openCylinder = new PShape();
PShape CoteCylinder = new PShape();


    
void CreateCyl() {
  float angle;
  float[] x = new float[cylinderResolution + 1]; 
  float[] y = new float[cylinderResolution + 1];
  
  //get the x and y position on a circle for all the sides
  for(int i = 0; i < x.length; i++) {
  angle = (TWO_PI / cylinderResolution) * i; x[i] = sin(angle) * cylinderBaseSize;
  y[i] = cos(angle) * cylinderBaseSize;
  }
  
  openCylinder = createShape();
  openCylinder.beginShape(TRIANGLE_STRIP);
  //draw the border of the cylinder
  for(int i = 0; i < x.length; i++) { 
    openCylinder.vertex(x[i], y[i] , 0);
    openCylinder.vertex(0,0,0);
  }
  
  for(int i = 0; i < x.length; i++) { 
    openCylinder.vertex(x[i], y[i], cylinderHeight);
    openCylinder.vertex(0,0,cylinderHeight);
  }
  
  openCylinder.endShape();
  
  CoteCylinder = createShape();
  CoteCylinder.beginShape(QUAD_STRIP);
  //draw the border of the cylinder
  for(int i = 0; i < x.length; i++) { 
    CoteCylinder.vertex(x[i], y[i] , 0);
    CoteCylinder.vertex(x[i], y[i], cylinderHeight);
  }
  CoteCylinder.endShape();
}

void displayCyl(){
  gameSurface.shape(openCylinder);
  gameSurface.shape(CoteCylinder);
}

void displayCylTopView(float x, float y, boolean b) {
  if(b) {
    topView.pushStyle();
    topView.fill(0, 0, 200);
    topView.ellipse((x-width/4)*180/400.0, (y-(height)/4)*180/400.0, 2*cylinderBaseSize*topViewLength/LongPlaque, 2*cylinderBaseSize*topViewLength/LongPlaque);
    topView.popStyle();
} else {
    topView.ellipse((x-width/4)*180/400.0, (y-(height)/4)*180/400.0, 2*cylinderBaseSize*topViewLength/LongPlaque, 2*cylinderBaseSize*topViewLength/LongPlaque);
  }
}
