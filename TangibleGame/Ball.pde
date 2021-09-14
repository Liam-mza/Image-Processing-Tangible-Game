class Ball {
  PVector location;
  PVector prevLocation;
  PVector velocity;
  PVector gravityForce;
  PVector dx;
  int radius;
  float gravityConstant=3.5;
  float normalForce = 32;
  float mu = 0.02;
  float frictionMagnitude = normalForce * mu;
  PShape sphere;
  float alphaX;
  float alphaZ;

  Ball(int r,int x,int y,int z,PImage img){
    radius=r;
    
    alphaX = 0;
    alphaZ = 0;
    
    location = new PVector(x,y,z);
    prevLocation = new PVector(location.x, location.y, location.z);
    dx = new PVector(location.x- prevLocation.x, location.z-prevLocation.z);
    
    gravityForce= new PVector(sin(Zangle) * gravityConstant,sin(Xangle) * gravityConstant,0);
    
    velocity = new PVector(0,0,0);
    
    sphere = createShape(SPHERE, r);
    sphere.setStroke(false);
    sphere.setTexture(img);
  }
  
  void update() {
        
        gravityForce= new PVector(sin(Zangle) * gravityConstant,0,sin(-Xangle) * gravityConstant);
        
        PVector friction = velocity.copy();
        friction.mult(-1);
        friction.normalize(); 
        friction.mult(frictionMagnitude);
        
        PVector regarde = new PVector(0, 0, -1);
        PVector cylBall = (new PVector((location.x-(ps.origin.x-width/2)), 0, (location.z-(ps.origin.y-height/2)))).normalize();
        thetaa = PVector.angleBetween(regarde, cylBall);
        if(location.x-(ps.origin.x-width/2) < 0) {
          thetaa = -thetaa; 
        }
        
        velocity.add(friction);
        velocity.add(gravityForce);
        location.add(velocity);
      }
      
  void display() {
        gameSurface.translate(location.x,location.y,location.z);
        sphere.rotate(alphaX, 0, 0, 1.0);
        sphere.rotateX(-alphaZ);
        gameSurface.shape(sphere);
  }
  
  void checkEdges() {
    if (location.x >= LongPlaque/2-radius) {
          velocity.x = -velocity.x;
          location.x= (LongPlaque/2)-radius;
        }
    else if (location.x <= -LongPlaque/2+radius) {
          velocity.x = -velocity.x;
          location.x=(-LongPlaque/2)+radius;
        }
    if (location.z >= LongPlaque/2-radius) {
          velocity.z = -(velocity.z);
          location.z= (LongPlaque)/2-radius;
        }
    else if (location.z <= -LongPlaque/2+radius) {
          velocity.z = -(velocity.z);
          location.z= (-LongPlaque)/2+radius;
        }
  } 
  
  
  void checkCylinderCollision(){
    for (PVector p: ps.particles){
      //shifting the coordonates to be as the ball coordinates ((0,0) at the center)
      PVector newp=  new PVector(-width/2+p.x,-40,-height/2+p.y);
      float dist= newp.dist(location);
      if (dist<=radius+cylinderBaseSize){ 
        //creating the normal vector
        PVector n = new PVector(location.x-newp.x,0,location.z-newp.z);
        n.normalize();
        
        //Avoiding the ball to enter the cylinder
        location.add(PVector.mult(n,radius+cylinderBaseSize-dist));
        
        //Calculating v2
        float dotV1 = velocity.dot(n);
        velocity=PVector.sub(velocity,n.mult(2*dotV1));
      }
    }
  }
  
  void rotateTexture() {
    dx.x = (location.x - prevLocation.x);
    dx.y = (location.z - prevLocation.z);
    alphaX = dx.x/(2*radius);
    alphaZ = dx.y/(2*radius);
    prevLocation.x = location.x;
    prevLocation.z = location.z;
  }
  
  
}
