class TaskLocations extends Datafile {
  TaskLocations(String filename, Experiment expr) {
    super(filename, expr);
  }

  void processLine(String[] e) {
    int update = int(e[0]);
    int x = int(e[1]);
    int y = int(e[2]);
    String task = e[3];

    if (task.equals("not")) {
      expr.grid.node(x, y, 0.5, 0.5, GREEN);
    } 
    else {
      expr.grid.node(x, y, 0.5, 0.5, RED);
    }
  }
}

