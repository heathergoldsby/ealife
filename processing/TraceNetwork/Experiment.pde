int[] INFO_AREA= {
  10, 430, 400, 400
};
PFont FONT;
int FONT_SIZE=16;

class Experiment {
  int currentUpdate;
  int lineCount;
  boolean paused;
  Grid grid; // the activity-displaying grid.

  ArrayList<Datafile> datafiles;

  Experiment() {    
    FONT = createFont("SansSerif", FONT_SIZE);
    grid = new Grid(this);
    datafiles = new ArrayList<Datafile>();
  }

  void addDatafile(Datafile df) {
    datafiles.add(df);
  }

  boolean atBegin() {
    return (currentUpdate <= 0);
  }

  boolean atEnd() {
    boolean end=true;
    for(Datafile i : datafiles) {
      end = end && i.atEnd(currentUpdate);
    }
    return end;
  }
  
  void reset() {
    currentUpdate=0;
    paused = true;
    redraw();
  }

  void step(int direction) {
    // are we doing anything?
    if((direction > 0) && atEnd()) {
      pause();
      return;
    }
    if((direction < 0) && atBegin()) {
      pause();
      return;
    }
    currentUpdate += direction;
    redraw();
  }

  void draw() {
    grid.draw(currentUpdate);
    for(Datafile i : datafiles) {
      i.draw(currentUpdate);
    }

    stroke(0);
    fill(0);
    textFont(FONT);
    lineCount=0;
    infoLine("Update: " + nf(currentUpdate,0));
    infoLine("");
    infoLine("m/M=start/stop movie");
    infoLine("p=toggle pause");
    infoLine("q=quit");
    infoLine("r=reset");
    infoLine("s/S=step fwd/back");
  }

  void infoLine(String l) {
    text(l, INFO_AREA[0], INFO_AREA[1]+FONT_SIZE*lineCount++);
  }

  void togglePause() {
    if(paused) {
      unpause();
    } 
    else {
      pause();
    }
  }

  void pause() {
    paused = true;
    noLoop();
  }

  void unpause() {
    paused = false;
    loop();
  }

  void mouseClicked(int x, int y) {
    grid.mouseClicked(x,y);
    redraw();
  }
}

