import java.util.ArrayList;
import java.util.List; 
import java.util.TreeSet;
import java.util.Collections;
import java.util.LinkedList;
import processing.video.*;
import gab.opencv.*;
import java.util.Map;
import java.util.Comparator;
import processing.core.PVector;
import org.opencv.core.Mat;
import org.opencv.core.CvType;
import org.opencv.core.Core;

PGraphics gameSurface;
PGraphics scoreBoard;
PGraphics topView;
PGraphics barChart;

ConfettiCanon endingScene;

float Xangle=0.0;
float Zangle=0.0;
float Speed=1.0;
float change=0.03;
int mouseTreshold=1;
PImage texture;
float thetaa = 0;
int score = 0;
int lastScore = score;
int scoreFactor = 1;
boolean beginTheGame = true; //Just for the score (used in ps.run())
int POINTS_TO_WIN = 60;

int topViewLength = 180;
int bigOrSmall;

boolean bigger = true;
boolean gameFinished = false;

float arrangeBaseCylWithPlaque;

int LongPlaque= 400;
int epaisseur=15;
int ballRadius = 30;

boolean editMode=false;
Ball ball;
HScrollbarGame scrollBar;
Button button;
ArrayList<PVector> coorCyl = new ArrayList<PVector>();
LinkedList<int[]> dataChart = new LinkedList<int[]>();

ParticleSystem ps;

PShape robotnik;
PImage robotnikImg;

ImageProcessing imgproc;
/* Set boardInfluence to false to be able to play without the image processing influence */
boolean boardInfluence = true;
Movie cam;

void settings() {
  size(800, 800, P3D);
}

void setup() {
  gameSurface = createGraphics(width, height-200, P3D);
  scoreBoard = createGraphics(topViewLength, topViewLength, P2D);
  topView = createGraphics(topViewLength, topViewLength, P2D);
  barChart = createGraphics(chartSurfaceLength, topViewLength, P2D);
  endingScene = new ConfettiCanon(width, height, 0, 0);
  
  CreateCyl();
  ps = new ParticleSystem();
  scrollBar = new HScrollbarGame(410, 765, 250, 15);
  button = new Button(720, 762, 20, 20);
  robotnikImg = loadImage("robotnik.png");
  robotnik = loadShape("robotnik.obj");
  robotnik.setStroke(false);
  robotnik.setTexture(robotnikImg);
  texture = loadImage("space2.jpg"); //Choose between earth.jpg or space2.jpg
  ball = new Ball(ballRadius, 0, -37, 0, texture);
  arrangeBaseCylWithPlaque = 8;
  
  cam = new Movie(this, "testvideo.avi");
  cam.loop();
  imgproc = new ImageProcessing();
  String []args = {"Image processing window"};
  PApplet.runSketch(args, imgproc);
}

void draw() {
  if (gameFinished) {
    endingScene.drawConfettiCanon();
  } else {
    scrollBar.update();
    scrollBar.display();
    button.drawBut();
    drawGame();
    image(gameSurface, 0, 0);
    drawTop();
    image(topView, 10, width-200+10);
    drawScores();
    image(scoreBoard, 210, width-200+10);
    drawChart();
    image(barChart, 410, width-200+10);
  }
}

void drawTop() {
  topView.beginDraw();
  //topView.noStroke();
  topView.fill(0, 100, 130, 100);
  topView.rect(0, 0, topViewLength-1, topViewLength-1);

  topView.stroke(0);

  topView.pushStyle();
  topView.fill(150, 0, 0);
  topView.ellipse(topViewLength*ball.location.x/LongPlaque+topViewLength/2.0, topViewLength*ball.location.z/LongPlaque+topViewLength/2.0, 2*ballRadius*topViewLength/LongPlaque, 2*ballRadius*topViewLength/LongPlaque);
  topView.popStyle();

  topView.fill(255);
  ps.runTopView();
  topView.endDraw();
}

void drawScores() {
  scoreBoard.beginDraw();
  scoreBoard.fill(255);
  scoreBoard.rect(0, 0, topViewLength-1, topViewLength-1);

  scoreBoard.pushStyle();
  scoreBoard.fill(220);
  scoreBoard.noStroke();
  scoreBoard.rect(3, 3, topViewLength-6, topViewLength-6);
  scoreBoard.popStyle();

  scoreBoard.fill(0);
  scoreBoard.textSize(16);
  scoreBoard.text("Current Score: ", 10, 22);
  String totalScore = ""+score;
  scoreBoard.text(totalScore, 10, 39);
  scoreBoard.text("Velocity: ", 10, 85);

  String velocity = ""+round(ball.dx.mag());
  scoreBoard.text(velocity, 10, 102);
  scoreBoard.text("Last Score: ", 10, 80+73);
  String lstScore = ""+lastScore;
  scoreBoard.text(lstScore, 10, 80+73+17);
  scoreBoard.endDraw();
}

float clampAng(float i) {
  if(i > PI/2) {
    return i-PI;
  }
  if(i < -PI/2) {
    return i+PI;
  }
  return i;
}

void drawGame() {
  gameSurface.beginDraw();
  //CAMERA AND LIGHTS
  if (bigger) {
    gameSurface.camera(width/2, height/2-50, 700, 400, 400, 0, 0, 1, 0);
  } else {
    gameSurface.camera(width/2, height/2-50, 800, 250, 250, 0, 0, 1, 0);
  }
  gameSurface.directionalLight(100, 125, 150, 0, 1, -1); 
  gameSurface.ambientLight(102, 102, 102);
  gameSurface.background(255);

  if (!editMode) {

    //DISPLAY ANGLE AND SPEED
    displayAnglesAndSpeed();

    //TRANSFORMATION
    gameSurface.translate(width/2, height/2, 0);
    if(boardInfluence) {
      Xangle = constrain(clampAng(imgproc.tilt.x), -PI/3, PI/3);
      Zangle = constrain(clampAng(imgproc.tilt.y), -PI/3, PI/3);
    }
    gameSurface.rotateZ(Zangle);
    gameSurface.rotateX(Xangle);

    //cylinder
    float ratio=(frameRate)/2;
    if (ratio==0) { 
      ratio=1;
    }
    if (frameCount%round(ratio)==0) {
      ps.addParticle();
    }
    ps.run();

    //PLATE
    gameSurface.pushMatrix();
    gameSurface.strokeWeight(0);

    gameSurface.fill(255);
    gameSurface.box(LongPlaque, epaisseur, LongPlaque);
    gameSurface.popMatrix();

    //COORDONATES LINES
    //drawCoordinateLines();

    //DISPLAY THE BALL
    gameSurface.pushMatrix();
    ball.update();
    ball.checkCylinderCollision();
    ball.checkEdges();
    ball.rotateTexture();
    ball.display();
    gameSurface.popMatrix();
  } else {
    if (bigger) {
      gameSurface.camera(width/2, height/2, 570, 400, 400, 0, 0, 1, 0);
    } else {
      gameSurface.camera(width/2, height/2, 500, 250, 250, 0, 0, 1, 0);
    }
    gameSurface.translate(width/2, height/2, 0);
    gameSurface.rotateX(PI/2);
    gameSurface.pushMatrix();
    gameSurface.strokeWeight(0);
    gameSurface.fill(255);
    gameSurface.box(LongPlaque, epaisseur, LongPlaque);
    gameSurface.popMatrix();

    //Draw a shadow of the ball
    gameSurface.pushMatrix();
    gameSurface.noStroke();
    gameSurface.fill(0, 100);
    gameSurface.translate(ball.location.x, epaisseur, -ball.location.z);
    gameSurface.sphere(ball.radius);
    gameSurface.popMatrix();

    if (!ps.particles.isEmpty()) {
      gameSurface.pushMatrix();
      gameSurface.rotateX(PI/2);
      gameSurface.translate(-width/2+ps.origin.x, gameSurface.height/2-ps.origin.y+100, -50);
      displayCyl();
      gameSurface.popMatrix();
    }
  }
  gameSurface.endDraw();
}

void mouseDragged() {
  if (mouseY <= 600) {
    if (pmouseY-mouseY>mouseTreshold) {
      Xangle=constrain(Xangle+change*Speed, -PI/3, PI/3);
    } else if (pmouseY-mouseY<-mouseTreshold) {
      Xangle=constrain(Xangle-change*Speed, -PI/3, PI/3);
    }
    if (pmouseX-mouseX>mouseTreshold) {
      Zangle=constrain(Zangle-change*Speed, -PI/3, PI/3);
    } else if (pmouseX-mouseX<-mouseTreshold) {
      Zangle=constrain(Zangle+change*Speed, -PI/3, PI/3);
    }
  }
}


void mouseWheel(MouseEvent event) {
  float e = event.getCount();
  if (e<0) {
    Speed +=0.1;
  } else if (e>0) {
    Speed -=0.1;
  }
  if (Speed<0) {
    Speed=0;
  }
  if (Speed>3) {
    Speed=3;
  }
}

void keyPressed() { 
  if (key == CODED) {
    if (keyCode == SHIFT && !editMode) {
      editMode=true;
    }
  }
}

void keyReleased() {
  if (key == CODED) {
    if (keyCode == SHIFT && editMode) {
      editMode=false;
    }
  }
}

void mouseClicked() {
  if (!bigger) {
    if (editMode) {
      if (mouseX>=75 && mouseX<=425 && mouseY>=75 && mouseY<=425) {
        PVector p = new PVector(mouseX, mouseY);
        PVector pball = new PVector(width/2+ball.location.x, height/2+ball.location.z);
        if (p.dist(pball) > (cylinderBaseSize+ball.radius)) {
          ps=new ParticleSystem(new PVector(mouseX, mouseY));
          coorCyl.add(new PVector(mouseX, mouseY));
        }
      }
    }
  } else {
    if (editMode) {
      if (mouseX>=225 && mouseX<=575 && (mouseY+100)>=225 && (mouseY+100)<=575) {
        PVector p = new PVector(mouseX, mouseY+100);
        PVector pball = new PVector(width/2+ball.location.x, (height)/2+ball.location.z);
        if (p.dist(pball) > (cylinderBaseSize+ball.radius)) {
          ps=new ParticleSystem(p);
          coorCyl.add(p);
        }
      }
    }
    button.updateOnMouseClicked();
  }
}

void displayAnglesAndSpeed() {
  if (!bigger) {
    gameSurface.pushMatrix();
    gameSurface.textMode(MODEL);
    gameSurface.textSize(20);
    String sx= "RotationX: "+(Xangle*180/PI);
    gameSurface.text(sx, -175+35, -175, 0);
    String sz= "RotationZ: "+(Zangle*180/PI);
    gameSurface.text(sz, 75+75+35, -175, 0);
    String s= "Speed: "+Speed;
    gameSurface.text(s, 325+150+35, -175, 0);
    gameSurface.popMatrix();
  } else {
    gameSurface.pushMatrix();
    gameSurface.textMode(MODEL);
    gameSurface.textSize(17);

    gameSurface.fill(0);
    gameSurface.text("RotationX: ", 5+100, 25, 0);
    gameSurface.fill(0, 80, 110);
    String sx = ""+(Xangle*180/PI);
    gameSurface.text(sx, 5+190, 25, 0);

    gameSurface.fill(0);
    gameSurface.text("RotationZ: ", 250+100, 25, 0);
    gameSurface.fill(0, 80, 110);
    String sz = ""+(Zangle*180/PI);
    gameSurface.text(sz, 250+190, 25, 0);

    gameSurface.fill(0);
    gameSurface.text("Speed: ", 495+100, 25, 0);
    gameSurface.fill(0, 100, 130);
    String s = ""+Speed;
    gameSurface.text(s, 495+160, 25, 0);

    gameSurface.popMatrix();
  }
}

void drawCoordinateLines() {
  gameSurface.pushMatrix();
  gameSurface.strokeWeight(2);

  gameSurface.textMode(SHAPE);
  gameSurface.textSize(40);

  gameSurface.stroke(255, 0, 0);
  gameSurface.fill(255, 0, 0);
  gameSurface.text('X', 400, -5, 0);
  gameSurface.line(-400, 0, 0, 400, 0, 0);

  gameSurface.stroke(0, 255, 0);
  gameSurface.fill(0, 255, 0);
  gameSurface.text('Y', 8, 300, 0);
  gameSurface.line(0, -300, 0, 0, 300, 0);

  gameSurface.stroke(0, 0, 255);
  gameSurface.fill(0, 0, 255);
  gameSurface.textSize(30);
  gameSurface.text('Z', 4, 0, 400);
  gameSurface.line(0, 0, -400, 0, 0, 400);
  gameSurface.popMatrix();
}
