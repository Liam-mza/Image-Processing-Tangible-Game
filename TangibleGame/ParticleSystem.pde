// A class to describe a group of Particles
class ParticleSystem {
ArrayList<PVector> particles; 
PVector origin;

ParticleSystem(PVector origin) {
this.origin = origin.copy();
particles = new ArrayList<PVector>();
particles.add(origin);
}

//A constructor to intialise an empty ParticleSystem that do nothing to wait the user to create one
ParticleSystem() {
this.origin = new PVector(0,0);
particles = new ArrayList<PVector>();
}

void addParticle() {
  PVector center;
  int numAttempts = 100;
  if (!particles.isEmpty()){
    score -= 6;
    for(int i=0; i<numAttempts; i++) {
      // Pick a cylinder and its center.
      int index = int(random(particles.size())); 
      center = particles.get(index).copy();
      // Try to add an adjacent cylinder.
      float angle = random(TWO_PI);
      center.x += sin(angle) * 2*cylinderBaseSize; 
      center.y += cos(angle) * 2*cylinderBaseSize; 
      if(checkPosition(center)) {
      particles.add(center);
      break; 
      }
    } 
  }
}

// Check if a position is available, i.e.
// - would not overlap with particles that are already created 
// (for each particle, call checkOverlap())
// - is inside the board boundaries
boolean checkPosition(PVector center) {
float maxdist=LongPlaque/2;
  
//check for overlap between two particules  
PVector rectOrigin= new PVector(width/2, height/2);
for (PVector p: particles){
  if (checkOverlap(center,p)){
  return false;
  }
}

//check that the particule is on the board
float distX= center.x-rectOrigin.x;
float distY= center.y-rectOrigin.y;
if (abs(distX)>maxdist-cylinderBaseSize||abs(distY)>maxdist-cylinderBaseSize){ 
return false;
}

//check overlap with the ball
PVector pball= new PVector(width/2+ball.location.x,height/2+ball.location.z);
if (center.dist(pball) < (cylinderBaseSize+ball.radius)){return false;}

return true;
}
// Check if a particle with center c1
// and another particle with center c2 overlap.
//return true s'il y un overlap
boolean checkOverlap(PVector c1, PVector c2) {
if( c1.dist(c2)<=2*cylinderBaseSize){
  return true;
}
else {
return false;}
}


// Iteratively update and display every particle,
// and remove them from the list if their lifetime is over. 
void run() {
  for (int i=0; i<particles.size();++i){
    PVector p = particles.get(i);
    if (eliminated(p)) {
      score += scoreFactor*round(ball.dx.mag());
      if(score > POINTS_TO_WIN) { gameFinished = true; }
      if (p.x==origin.x && p.y==origin.y) {
        lastScore = score;
        if(beginTheGame) {
          lastScore = score;
          beginTheGame = false;
        }
        score = 0;
        particles = new ArrayList<PVector>();
      }
      else {
        particles.remove(i);
      }
    } else {
      gameSurface.pushMatrix();
      gameSurface.rotateX(PI/2);
      gameSurface.translate(-width/2+p.x,-height/2+p.y, arrangeBaseCylWithPlaque);
      displayCyl();
      gameSurface.popMatrix();
      gameSurface.pushMatrix();
      gameSurface.translate(-width/2+origin.x, -epaisseur+cylinderHeight-73,-height/2+origin.y);
      gameSurface.rotateX(PI);
      gameSurface.rotateY(thetaa);
      gameSurface.scale(50);
      gameSurface.shape(robotnik);
      gameSurface.popMatrix();
    }
  }
}

void runTopView() {
  for (int i=0; i<particles.size();++i){
    PVector p = particles.get(i);
    if (!eliminated(p)) {
      if (p.x==origin.x && p.y==origin.y) {
        displayCylTopView(p.x, p.y, true);
      } else {
        displayCylTopView(p.x, p.y, false);
      }
    }
  }
}

boolean eliminated(PVector p){
  PVector newp=  new PVector(-width/2+p.x,-40,-height/2+p.y);
  float dist= newp.dist(ball.location);
  return dist<=ball.radius+cylinderBaseSize;
}
}
