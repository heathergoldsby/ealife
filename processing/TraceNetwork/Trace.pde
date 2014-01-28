class Trace extends Datafile {
  Trace(String filename, Experiment expr) {
    super(filename, expr);
  }      

  // draw line e
  void processLine(String[] e) {
    int update = int(e[0]);
    int x = int(e[1]);
    int y = int(e[2]);
    int heading = int(e[3]);
    int germ = int(e[4]);
    String inst = trim(e[5]);
    int ip = int(e[6]);

    // occupied cell:
//    expr.grid.occupied(x, y, LIGHT_GREY);

    // yellow dot == broadcast
    if (inst.equals("bc_msg")) {
      expr.grid.node(x, y, 0.75, 0.25, YELLOW);
    }

    // white dot == rx
    if(inst.equals("rx_msg")) {
      expr.grid.node(x, y, 0.25, 0.75, WHITE);
    }

    // blue dot == soma
    if (germ != 1) {
      expr.grid.node(x, y, 0.25, 0.25, BLUE);
    }

    if (inst.equals("tx_msg")) {
      expr.grid.cellHeading(x, y, heading, GREEN);
    } 
    else {
      expr.grid.cellHeading(x, y, heading, BLACK);
    }
    
    expr.grid.traceDot(x,y,ip);
  }
}

