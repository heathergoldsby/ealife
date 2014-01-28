class Datafile {
  HashMap<Integer,ArrayList<String[]>> data;
  int end;
  Experiment expr;

  Datafile(String filename, Experiment expr) {
    this.expr = expr;
    print("loading datafile: " + filename + "...");
    String[] file = loadStrings(filename);
    data = new HashMap<Integer,ArrayList<String[]>>();

    for(int i=1; i<file.length; ++i) {
      String l = trim(file[i]);
      if((l.length()==0) || (l.charAt(0)=='#')) {
        continue;
      }
      String[] e = splitTokens(l);
      processLineOnLoad(e);
    }
    println(" done.");
  }

  boolean atEnd(int update) {
    return update >= end;
  }
  
  void processLineOnLoad(String[] e) {
    int k = int(e[0]); // update
    end = max(end, k);
    if(!data.containsKey(k)) {
      data.put(k, new ArrayList<String[]>());
    }
    data.get(k).add(e);
  }

  void processLine(String[] e) {
  }

  void topHalf(int update) {
  }

  void bottomHalf(int update) {
  }

  void processAll(ArrayList<String[]> all) {
  }

  void draw(int update) {
    topHalf(update);
    if(data.containsKey(update)) {
      ArrayList<String[]> ar = data.get(update);
      processAll(ar);
      for(int i=0; i<ar.size(); ++i) {
        processLine(ar.get(i));
      }
    }
    bottomHalf(update);
  }
}

