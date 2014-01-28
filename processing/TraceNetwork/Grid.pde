int[] GRID_AREA = {
  10, 10, 400, 400
};
int[] GRID_CELLS = { 
  5, 5
};
int[] CELL_SIZE = {
  80, 80
};

int[] GENOME_AREA= {
  430, 15, 8
};

/*! 
 */
class Grid {
  Experiment expr;
  int tracex, tracey;

  Grid(Experiment expr) {
    this.expr = expr;
    tracex=-1;
    tracey=-1;
  }

  void gridTranslate() {
    pushMatrix();
    translate(GRID_AREA[0], GRID_AREA[1]);
  }

  void cellTranslate(float offset) {
    pushMatrix();
    translate(GRID_AREA[0]+offset*CELL_SIZE[0], GRID_AREA[1]+offset*CELL_SIZE[1]);
  }

  void draw(int update) {
    gridTranslate();
    stroke(150);
    for (int y=0; y<GRID_CELLS[1]; ++y) {
      for (int x=0; x<GRID_CELLS[0]; ++x) {
        pushMatrix();    
        translate(x*CELL_SIZE[0], y*CELL_SIZE[1]);
        fill(240);
        rect(0, 0, CELL_SIZE[0], CELL_SIZE[1]);
        popMatrix();
      }
    }
    popMatrix();
  }

  float[] cell2xy(int cell) {
    float[] xy = { 
      cell%GRID_CELLS[0]*CELL_SIZE[0], cell/GRID_CELLS[0]*CELL_SIZE[1]
    };  
    return xy;
  }

  int idx2cell(int x, int y) {
    return y*GRID_CELLS[0] + x;
  }

  int xy2cell(int x, int y) {
    x -= GRID_AREA[0];
    y -= GRID_AREA[1];
    x = x / (GRID_AREA[2] / GRID_CELLS[0]);
    y = y / (GRID_AREA[3] / GRID_CELLS[1]);
    return y * GRID_CELLS[0] + x;
  }

  boolean inGrid(int x, int y) {
    return (x > GRID_AREA[0]) && (x<(GRID_AREA[0]+GRID_AREA[2])) && (y > GRID_AREA[1]) && (y<(GRID_AREA[1]+GRID_AREA[3]));
  }

  void toggleTrace(int x, int y) {
    if ((tracex==x) && (tracey==y)) {
      tracex = -1; 
      tracey=-1;
    } 
    else {
      tracex = x; 
      tracey = y;
    }
  }
  void mouseClicked(int x, int y) {
    if (inGrid(x, y)) {
      x -= GRID_AREA[0];
      y -= GRID_AREA[1];
      x = x / (GRID_AREA[2] / GRID_CELLS[0]);
      y = y / (GRID_AREA[3] / GRID_CELLS[1]);
      toggleTrace(x, y);
    }
  }


  void cellTranslate(int x, int y, float xoffset, float yoffset) {
    pushMatrix();
    translate(GRID_AREA[0] + x*CELL_SIZE[0] + xoffset*CELL_SIZE[0], 
    GRID_AREA[1] + y*CELL_SIZE[1] + yoffset*CELL_SIZE[1]);
  }


  //  void link(int src, int dst, int[] clr) {
  //    cellTranslate(0.5);
  //    stroke(clr[0], clr[1], clr[2]);
  //    float[] s = cell2xy(src);
  //    float[] d = cell2xy(dst);
  //    strokeWeight(2);
  //    line(s[0], s[1], d[0], d[1]);
  //    strokeWeight(1);
  //    popMatrix();
  //  }

  void occupied(int x, int y, int[] clr) {
    // did we already occupy this cell?  if yes, don't recolor:
        
    cellTranslate(x, y, 0.0, 0.0);
    fill(clr[0], clr[1], clr[2]);
    if (tracing(x, y)) {
      strokeWeight(2);
      stroke(0);
      rect(1, 1, CELL_SIZE[0]-2, CELL_SIZE[1]-2);
    } 
    else {
      stroke(150);
      rect(0, 0, CELL_SIZE[0], CELL_SIZE[1]);
    }
    strokeWeight(1);
    popMatrix();
  }

  void node(int x, int y, float xoffset, float yoffset, int[] clr) {
    cellTranslate(x, y, xoffset, yoffset);
    stroke(0);
    fill(clr[0], clr[1], clr[2]);
    ellipse(0, 0, 20, 20);
    popMatrix();
  }

  void cellHeading(int x, int y, int heading, int[] clr) {
    cellTranslate(x, y, 0.5, 0.5);
    rotate(heading * TWO_PI/8.0);
    fill(clr[0], clr[1], clr[2]);
    rect(CELL_SIZE[0]/4, 0, CELL_SIZE[0]/4, CELL_SIZE[1]/20);
    popMatrix();
  }

  void genomeTranslate() {
    pushMatrix();
    translate(GENOME_AREA[0], GENOME_AREA[1]);
  }

  void genome(ArrayList<String> g) {
    genomeTranslate();
    translate(10, 0);
    stroke(0);
    fill(0);
    textFont(SMALL_FONT);
    for (int i=0; i<g.size(); ++i) {
      text(g.get(i), 0, GENOME_AREA[2]*i);
    }
    popMatrix();
  }

  boolean tracing(int x, int y) {
    return (tracex==x) && (tracey==y);
  }

  void traceDot(int x, int y, int ip) {
    if (tracing(x, y)) {
      genomeTranslate();
      translate(0, GENOME_AREA[2]*ip - GENOME_AREA[2]/2);
      fill(RED[0], RED[1], RED[2]);
      ellipse(0, 0, 6, 6);
      popMatrix();
    }
  }
}

