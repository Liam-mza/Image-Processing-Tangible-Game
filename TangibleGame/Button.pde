class Button {
  
  boolean butOver = false;
  boolean squares = true;
  int butX;
  int butY;
  int butWidth;
  int butHeight;
  
  Button(int x, int y, int w, int h) {
    butX = x;
    butY = y;
    butWidth = w;
    butHeight = h;
  }
  
  void drawBut() {
    updateBut();
    pushStyle();
    if(butOver) {
      fill(0, 50, 80);
    } else {
      fill(0, 100, 130);
    }
    rect(butX, butY, butWidth, butHeight);
    popStyle();
  }
  
  void updateOnMouseClicked() {
    if(butOver) {
      squares = !squares;
    }
  }
  
  void updateBut() {
    if(overBut(butX, butY, butWidth, butHeight)) {
      butOver = true;
    } else {
      butOver = false;
    }
  }

  boolean overBut(int x, int y, int wid, int heig)  {
    if (mouseX >= x && mouseX <= x+wid && 
        mouseY >= y && mouseY <= y+heig) {
      return true;
    } else {
      return false;
    }
  }
  
}
