int currTime = 0;
int timeScale = 8; //mid value
int currScale = timeScale;
float timer = millis();
int chartSurfaceLength = 2*topViewLength+20;
float translateAmount = 0;
int squareSize = 8;

int clamp(int value, int min, int max) {
  if(value < min) { return min; }
  if(value > max) { return max; }
  return value;
}

void drawChart() {
  barChart.beginDraw();
  barChart.stroke(0);
  barChart.fill(255);
  barChart.rect(0, 0, 2*topViewLength+20-1, topViewLength-40-1);
  barChart.translate(translateAmount, 0);
  if(millis() - timer >= 1000) {
    dataChart.add(new int[]{currTime, -clamp(score, -topViewLength/2+20+1, topViewLength/2-20)});
    drawChartBlocks();
    currTime = currScale*(dataChart.size()-1);
    currTime += currScale;
    if(currTime >= chartSurfaceLength) {
      translateAmount -= currScale;
    }
    timer = millis();
  } else {
    drawChartBlocks();
  }
  barChart.endDraw();
}

void drawChartBlocks() {
  barChart.pushStyle();
  barChart.fill(0, 100, 130);
  if(button.squares) {
    for(int i=0; i<dataChart.size(); ++i) {
      int tempScr = 0;
      int currScore = -dataChart.get(i)[1];
      if(currScore <= 0) {
        while(tempScr < -currScore-squareSize) {
          barChart.rect(dataChart.get(i)[0], topViewLength/2.0-20+tempScr, currScale, squareSize);
          tempScr += squareSize;
        }
        barChart.rect(dataChart.get(i)[0], topViewLength/2.0-20+tempScr, currScale, -currScore - tempScr);
      } else {
        while(tempScr > -currScore+squareSize) {
          barChart.rect(dataChart.get(i)[0], topViewLength/2.0-20+tempScr, currScale, -squareSize);
          tempScr -= squareSize;
        }
        barChart.rect(dataChart.get(i)[0], topViewLength/2.0-20+tempScr, currScale, -currScore - tempScr);
      }
      //COMMENT PREVIOUS LINES FROM "int tempScr = 0;" AND UNCOMMENT THIS LINE TO DRAW WITHOUT SQUARES AND GET SMOOTH RECTANGLES
      //barChart.rect(dataChart.get(i)[0], topViewLength/2.0-20, currScale, dataChart.get(i)[1]);
    }
  } else {
    for(int i=0; i<dataChart.size(); ++i) {
      barChart.rect(dataChart.get(i)[0], topViewLength/2.0-20, currScale, dataChart.get(i)[1]);
    }
  }
  barChart.popStyle();
}

void updateTimeScaleWithSlider(float sliderPos) {
  currScale = round(2*sliderPos*timeScale);
  int dataLength = dataChart.size();
  for(int i=0; i<dataLength; ++i) {
    //prevents from O(n^2) runtime
    int[] data = dataChart.remove();
    data[0] = currScale*i;
    dataChart.add(data);
  }
  translateAmount = min(chartSurfaceLength-currScale*(dataChart.size()), 0);
}
