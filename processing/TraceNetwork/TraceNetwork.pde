String ROOT="/Users/dk/tmp/gls/gls-var";
String EXPR="/d_3";
//String EXPR="/d_9";

int SIZE_X=800;
int SIZE_Y=850;
boolean MAKE_MOVIE=false;

Experiment expr;

// Setup the grid, get our list of directories:
void setup() {
  size(SIZE_X, SIZE_Y);
  frameRate(5);
  smooth();
  textAlign(LEFT);
  expr = new Experiment();
  expr.pause();
  expr.addDatafile(new Trace(ROOT+EXPR+"/trace.dat", expr));
  expr.addDatafile(new TaskLocations(ROOT+EXPR+"/task_locations.dat", expr));
  expr.addDatafile(new Genome(ROOT+EXPR+"/lod_last_genome.dat", expr));
}

// Draw (and potentially advance) the current update:
void draw() {
  background(255);
  if (!expr.paused) {
    expr.step(1);
  }
  expr.draw();
  if (MAKE_MOVIE) {
    saveFrame("frames/#####.png");
  }
}

// UI -- on key press
void keyPressed() {
  switch(key) {
  case 'M':
     {
       MAKE_MOVIE = false;
       break;
     }
  case 'm':
    {
      MAKE_MOVIE = true;
      expr.reset();
      expr.togglePause();
      break;
    }
  case 'P':
  case 'p': 
    {
      expr.togglePause();
      break;
    }
  case 'Q':
  case 'q': 
    {
      exit();
    }
  case 'r': 
    {
      expr.reset();
      break;
    }
  case 'S': 
    {
      expr.step(-1);
      break;
    }
  case 's': 
    {
      expr.step(1);
      break;
    }
  }
}


// UI - on mouse click
void mouseClicked() {
  expr.mouseClicked(mouseX, mouseY);
}

