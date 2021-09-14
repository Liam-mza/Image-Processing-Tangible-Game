PGraphics window;

class ConfettiCanon {
  
  int w;
  int h;
  int Xcoord;
  int Ycoord;
  ArrayList<Confetti> confettis;
  float timer;
  int count = 0;
  
  ConfettiCanon(int x, int y, int placeX, int placeY) {
    Xcoord = placeX;
    Ycoord = placeY;
    w = x;
    h = y;
    window = createGraphics(x, y, P2D);
    confettis = new ArrayList<Confetti>();
  }
  
  
  void drawConfettiCanon() {
    window.beginDraw();
    window.background(0);
    if(millis() - timer >= 10) {
      confettis.add(new Confetti(w/2, h/8, random(-8, 8), 0, random(0, 255), random(0, 255), random(0, 255)));
      timer = millis();
      //Number of confettis in ArrayList confettis converges to 599
      if(count >= 200) {
        count = 0;
        for(int i=0; i<confettis.size()/2; ++i) {
          confettis.remove(0);
        }
      }
      ++count;
    }
    for(Confetti p : confettis) {
      p.update();
      p.display();
    }
    window.endDraw();
    image(window, Xcoord, Ycoord);
  }
  
}



class Confetti { 
  PVector location;
  PVector velocity;
  PVector gravity;
  float red;
  float green;
  float blue;
  
  Confetti(float locationX, float locationY, float velocityX, float velocityY, float r, float g, float b) {
    location = new PVector(locationX, locationY);
    velocity = new PVector(velocityX, velocityY); 
    gravity = new PVector(0, 0.5);
    red = r;
    green = g;
    blue = b;
  }
  
  void update() {
    location.add(velocity);
    velocity.add(gravity);
  }
        
  void display() {
    window.noStroke();
    window.strokeWeight(2);
    window.fill(red, green, blue);
    window.ellipse(location.x, location.y, 24, 24);
  }
}
