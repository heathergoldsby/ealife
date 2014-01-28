PFont SMALL_FONT;
int SMALL_FONT_SIZE=8;


class Genome extends Datafile {
  ArrayList<String> genome;

  Genome(String filename, Experiment expr) {
    super(filename, expr);
    SMALL_FONT = createFont("Monospaced", SMALL_FONT_SIZE);
  }

  boolean atEnd(int update) {
    return true;
  }

  void processLineOnLoad(String[] e) {
    if (genome == null) {
      genome = new ArrayList<String>();
    }  
    int pos = int(e[0]);
    int isa = int(e[1]);
    String inst = e[2];
    genome.add(inst);
  }

  void draw(int update) {
    expr.grid.genome(genome);
  }
}

