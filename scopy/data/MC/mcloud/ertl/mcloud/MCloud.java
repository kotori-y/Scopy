
/*
 * Molecule Cloud, Copyright (c) 2012, Novartis Institutes for BioMedical Research Inc
 * written by Peter Ertl and Bernhard Rohde
 *
 * this software is released under the terms of FreeBSD license,
 * see license.txt in the distribution
 *
 * for usage details see readme.txt in the distribution
 *
 * cite as:
 * P. Ertl, B. Rohde: The Molecule Cloud - compact visualization of large collections of molecules
 * J. Cheminformatics 2012, 4:12
 * http://www.jcheminf.com/content/4/1/12/
 */

package ertl.mcloud;

import java.io.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.awt.Image;
import java.awt.image.*;
import java.awt.geom.*;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
class MCloud {
  
  JFrame molFrame = null;
  MolPanel molPanel;
  static int pwidth = 1000;
  static int pheight = 700;

  ArrayList boxes = new ArrayList(); 
  // options
  boolean doImage = false;
  boolean doOptimization = true;
  static String signature = null;
  int stop = 0; // debug stop
  int skip = 0;
  boolean noMolecules = false;
  boolean nologscale = false;
  boolean doGui = true;
  boolean debug = false;
  ArrayList layoutPoints = null;

  // --------------------------------------------------------------------------
  void process(String fileName, int upto) {
    
    readData(fileName,upto);
    doSizes();

    // sorting acording to the size (number of atoms)
    MolBox.sortType = MolBox.SIZE;
    Collections.sort(boxes);

    cloudLayout();

    // creating and showing MolPanel;
    if (doGui) {
      molFrame = new JFrame();
      molFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      molFrame.setTitle("Molecule Cloud");

      molPanel = new MolPanel(pwidth,pheight,this);

      molFrame.getContentPane().setLayout(new BorderLayout());
      molFrame.getContentPane().add("Center",molPanel);

      molFrame.pack();
      molFrame.setLocation(200,10);
      molFrame.setVisible(true); 

      molPanel.repaint();
    }

    if (doOptimization) optimizeLayout();

    if (!noMolecules) MolBox.showMolecules = true;
    if (doGui) molPanel.repaint();
    if (doImage) saveImage();
  }
  // --------------------------------------------------------------------------
  void optimizeLayout() {
    System.err.println("optimizing layout ...");
    double v;

    MolBox.sortType = MolBox.COUNT;
    Collections.sort(boxes); // with largest count on the top

    int nnotmoved = 0;
    int move = 1;
    for (int step=1;step<=2000;step++) {
      if (step % 100 == 0) System.err.println(step + "\r");

      boolean moved = false;
      for (int i=0;i<boxes.size();i++) {
        MolBox box = (MolBox)boxes.get(i);
  
        double vmin = score(i);
  
        double bestx=0.,besty=0.;
        boolean ok = false;
  
        box.x += move;
        // making sure not moving outside panel area
        if (box.x - box.width/2 > 0 && box.x + box.width/2 < pwidth) {
        v = score(i); 
        if (v < vmin) {
          vmin = v;
          bestx = box.x;
          besty = box.y;
          ok = true;
        }
        }
        box.x -= 2 * move;
        if (box.x - box.width/2 > 0 && box.x + box.width/2 < pwidth) {
        v = score(i); 
        if (v < vmin) {
          vmin = v;
          bestx = box.x;
          besty = box.y;
          ok = true;
        }
        }
        box.x += move;
  
        box.y += move;
        if (box.y - box.height/2 > 0 && box.y + box.height/2 < pheight) {
        v = score(i); 
        if (v < vmin) {
          vmin = v;
          bestx = box.x;
          besty = box.y;
          ok = true;
        }
        }
        box.y -= 2 * move;
  
        if (box.y - box.height/2 > 0 && box.y + box.height/2 < pheight) {
        v = score(i); 
        if (v < vmin) {
          vmin = v;
          bestx = box.x;
          besty = box.y;
          ok = true;
        }
        }
        box.y += move;
  
        if (ok) {
          box.x = bestx;
          box.y = besty;
          moved = true;
        }
  
      }

      if (doGui) molPanel.repaint();
      if (!moved) {
        if (++nnotmoved == 3) {
          if (debug) System.err.println("optimisation took " + step + " steps");
          break;
        }
      }
      else if (nnotmoved > 0) {
        System.err.println("reset nnm");
        nnotmoved = 0;
      }
    };

  }
  // --------------------------------------------------------------------------
  private double score(int n) {
    return score(n,false);
  }
  // --------------------------------------------------------------------------
  // calculates current layout score
  private double score(int n, boolean upto) {
    double v = 0.;

    // distances between boxes
    MolBox box = (MolBox)boxes.get(n);
    double sq1 = Math.sqrt(box.width/2*box.width/2 + box.height/2 * box.height/2);
    for (int i=0;i<boxes.size();i++) {
      if (i == n) {
        if (upto) break;
        else continue;
      }
      MolBox boxi = (MolBox)boxes.get(i);
      double sq2 = Math.sqrt(boxi.width/2*boxi.width/2 + boxi.height/2 * boxi.height/2);
      // distance, dx, dy (if larger > 0)
      double[] dd = distance(box,boxi);
      double limit = sq1+sq2;

      double ov = overlap(box,boxi);

      if (ov > 0) {
        v += ov * ov;
        double q = ov * ov;
      }
      else if (dd[0] < limit) {
        v += limit - dd[0];
      }

    }
    
    // borders have repulsive force
    double limit = MolBox.avrgWidth;
    double d,dd;
    d = box.x - box.w;
    if (d <= limit) {
      dd = (limit - d);
      v += dd*dd;
    }
    d = pwidth - box.x - box.w;
    if (d <= limit) {
      dd = (limit - d);
      v += dd*dd;
    }
    d = box.y - box.h;
    if (d <= limit) {
      dd = (limit - d);
      v += dd*dd;
    }
    d = pheight - box.y - box.h;
    if (d <= limit) {
      dd = (limit - d);
      v += dd*dd;
    }

    // additional push from corners
    limit = MolBox.avrgWidth * 3;
    double cs = 8.;
    d = distance(box,0,0);
    if (d <= limit) v += (limit-d) * (limit-d) * cs;
    d = distance(box,pwidth,0);
    if (d <= limit) v += (limit-d) * (limit-d) * cs;
    d = distance(box,0,pheight);
    if (d <= limit) v += (limit-d) * (limit-d) * cs;
    d = distance(box,pwidth,pheight);
    if (d <= limit) v += (limit-d) * (limit-d) * cs;

    return v;
  }
  // --------------------------------------------------------------------------
  private double[] distance(MolBox box1, MolBox box2) {
    double dx = box2.x - box1.x;
    double dy = box2.y - box1.y;
    double r12 = Math.sqrt(dx * dx + dy * dy); // distance between centers

    // dx and dy are distances betwen 2 boxes (borders)
    dx = 0.;
    if (box1.x + box1.w < box2.x - box2.w) 
      dx = box2.x - box2.w - (box1.x + box2.w);
    else if (box1.x - box1.w > box2.x + box2.w) 
      dx = box1.x - box2.w - (box1.x + box2.w);
    dy = 0.;
    if (box1.y + box1.h < box2.y - box2.h) 
      dy = box2.y - box2.h - (box1.y + box2.h);
    else if (box1.y - box1.h > box2.y + box2.h) 
      dy = box1.y - box2.h - (box1.y + box2.h);


    double[] d = new double[3];
    d[0] = r12;
    d[1] = dx;
    d[2] = dy; 
    return d;
  }
  // --------------------------------------------------------------------------
  private double distance(MolBox box, int x, int y) {
    double dx = box.x - x;
    double dy = box.y - y;
    double r = Math.sqrt(dx * dx + dy * dy);
    return r;
  }
  // --------------------------------------------------------------------------
  private double overlap(MolBox box1, MolBox box2) {
    if (box1.x + box1.w < box2.x - box2.w) return 0.; 
    if (box1.x - box1.w > box2.x + box2.w) return 0.; 
    if (box1.y + box1.h < box2.y - box2.h) return 0.;
    if (box1.y - box1.h > box2.y + box2.h) return 0.;

    // there is some overlap
    double xov = Math.min(box1.x+box1.w,box2.x+box2.w) - Math.max(box1.x-box1.w,box2.x-box2.w);
    double yov = Math.min(box1.y+box1.h,box2.y+box2.h) - Math.max(box1.y-box1.h,box2.y-box2.h);

    return xov * yov;
  }
  // --------------------------------------------------------------------------
  // layout with largest box in the center and smaller ones around
  void cloudLayout() {

    System.err.println("initial layout ...");
    // largest first
    MolBox.sortType = MolBox.SIZE;
    Collections.sort(boxes);

    int n = -1;
    for (Iterator i=boxes.iterator();i.hasNext();) {
      n++;
      MolBox box = (MolBox)i.next();
      if (n == 0) {
        // the largest box is located in the center, slightly left
        box.setCenter(pwidth/2-box.w,pheight/2);
      }
      else {
        // iterative placement of other boxes
        double vbest = Double.MAX_VALUE;
        int xbest=0,ybest=0;
        double step = 15.;
        for (int r=5;r<Math.sqrt(pwidth*pwidth + pheight*pheight)/2.;r+=step) {
          int nslices = (int)(2. * Math.PI * r / step);
          if (nslices < 1) nslices = 1;
          double sliceangle = 2. * Math.PI / nslices;
          for (int k=0;k<nslices;k++) {
            double angle = k * sliceangle;
            int x = (int)(pwidth/2. + r * Math.sin(angle));
            int y = (int)(pheight/2. + r * Math.cos(angle));
            if (stop > 0 && n == 1) layoutPoints.add(new Point(x,y)); // debug
            box.setCenter(x,y);
            if (isOutside(box,x,y)) continue;
            double v = score(n,true);
            if (v < vbest) {
              vbest = v;
              xbest = x;
              ybest = y;
            } 
          }
        }
        box.setCenter(xbest,ybest);
        if (stop > 0 && n >= stop) {
          box.setCenter(-1000,-1000);
        }
      }
    }
  }
  // --------------------------------------------------------------------------
  boolean isOutside(MolBox box, int x, int y) {
    if (box.x - box.w <= 0) return true;
    if (box.x + box.w >= pwidth) return true;
    if (box.y - box.h <= 0) return true;
    if (box.y + box.h >= pheight) return true;
    return false;
  }
  // --------------------------------------------------------------------------
  public static void main(String[] args) {

    System.err.println("MoleculeCloud v1.21");
    if (args.length == 0) {
      System.err.println("Usage: java -jar MCloud.jar -f data [-parameters]");
      System.err.println("data contains smiles tab separated from its frequency, one record per line");
      System.err.println("the most important parameters are:");
      System.err.println("-i = save image into mcloud.png");
      System.err.println("-x value -y value = image dimensions");
      System.err.println("-nogui = runs without gui (to generate large images)");
      System.err.println("-n = process only first n molecules from the data file");
      System.err.println("-skip n = skip first n structures (usually used 1 to skip large phenyl)");
      System.exit(0);
    }

    MCloud mcloud = new MCloud();
    int upto = Integer.MAX_VALUE;

    String fileName = null;
    // parameters
    if (args.length == 1)
      fileName = args[0];
    else {
      for (int i=0;i<args.length;i++) {
        if (args[i].startsWith("-f")) fileName = args[++i];
        else if (args[i].startsWith("-i")) mcloud.doImage = true;
        else if (args[i].equals("-noopt")) mcloud.doOptimization = false;
        else if (args[i].equals("-nologscale")) mcloud.nologscale = true;
        else if (args[i].equals("-n")) upto = Integer.parseInt(args[++i]); 
        else if (args[i].equals("-x")) pwidth = Integer.parseInt(args[++i]);
        else if (args[i].equals("-stop")) {
          mcloud.stop = Integer.parseInt(args[++i]);
          mcloud.layoutPoints = new ArrayList();
        }
        else if (args[i].equals("-y")) pheight = Integer.parseInt(args[++i]);
        else if (args[i].equals("-skip")) mcloud.skip = Integer.parseInt(args[++i]);
        else if (args[i].equals("-nogui")) mcloud.doGui = false;
      }
    }

    if (fileName == null) {
      System.err.println("no file name");
    }


    mcloud.process(fileName,upto);

  }
  // --------------------------------------------------------------------------
  void readData(String fileName, int upto) {

    double maxScale = 0.;
    double minScale = Double.MAX_VALUE;
    MolBox.avrgWidth = MolBox.avrgHeight = 0.;
    try {
      BufferedReader in = new BufferedReader(new FileReader(fileName));
      String line;
      int n = 0;
      while ((line = in.readLine()) != null) {
        n++;
        StringTokenizer st = new StringTokenizer(line,"\t");
        String smiles = st.nextToken();
        if (skip > 0 && n <= skip) {
          System.err.println("skipping " + smiles + " "+n);
          continue;
        }
        double scale = Double.parseDouble(st.nextToken()); // count
        if (!nologscale) scale = Math.log(scale);
        if (scale > maxScale) maxScale = scale;
        if (scale < minScale) minScale = scale;
        String data = null;
        if (st.hasMoreTokens()) data = st.nextToken(); // coloring etc

        MolBox box = new MolBox(smiles,scale,data);
        boxes.add(box);

        if (n == upto) break; // reading only upto items
      }
    }
    catch (IOException e) {
      System.err.println("Cannot read file " + fileName +  "!");
      e.printStackTrace();
      System.exit(-1);
    }

    // largest molecule will be maxfactor-times larger than the smallest one
    double maxfactor = 5.;
    for (Iterator i=boxes.iterator();i.hasNext();) {
      MolBox box = (MolBox)i.next();
      double scale = (box.scale - minScale) / (maxScale - minScale);
      scale = 1. + (maxfactor - 1.) * scale;
      box.rescale(scale*scale);
    }
  }
  // --------------------------------------------------------------------------
  void doSizes() { 

    System.err.println("calculating sizes ...");
    double idealdensity = 85.;

    // 1st pass, getting molecule dimensions in any standard scale
    double sumsize = 0.;
    for (Iterator i=boxes.iterator();i.hasNext();) {
      MolBox box = (MolBox)i.next();
      double size = box.width * box.height;
      sumsize += size;
    }
    double pv = sumsize * 100. / (pwidth * pheight); 
    double rescale = idealdensity / pv;
    if (debug) System.err.println("pv0="+pv+", rescale="+rescale);

    for (int step=1;step<=10;step++) { // scale iteration steps

      // 2nd pass, now with corrected scale
      MolBox.avrgWidth = 0.; MolBox.avrgHeight = 0.;
      sumsize = 0;
      for (Iterator i=boxes.iterator();i.hasNext();) {
        MolBox box = (MolBox)i.next();
        // new scale
        box.rescale(rescale);
        MolBox.avrgWidth += box.width;
        MolBox.avrgHeight += box.height;

        double size = box.width * box.height;
        sumsize += size;
      }

      pv = sumsize * 100. / (pwidth * pheight); 
      rescale = idealdensity / pv;
      //if (pv > idealdensity) rescale *= 0.9;
      //else rescale *= 1.1;
      if (debug) System.err.println("corrected pv="+pv+", rescale="+rescale);

      MolBox.avrgWidth /= boxes.size();
      MolBox.avrgHeight /= boxes.size();

      if (Math.abs(idealdensity-pv) < 1.) break;
    }

  }
  // --------------------------------------------------------------------------
  void draw(Graphics2D g) {

    // display of layout points for testing
    if (layoutPoints != null) {
      for (Iterator i=layoutPoints.iterator();i.hasNext();) {
        Point p = (Point)i.next();
        //g.drawLine(p.getX(),p.getY(),p.getX(),p.getY());
        g.setColor(Color.black);
        g.drawLine(p.x-2,p.y,p.x+2,p.y);
        g.drawLine(p.x,p.y-2,p.x,p.y+2);
      }
    }    

    for (Iterator i=boxes.iterator();i.hasNext();) {
      MolBox box = (MolBox)i.next();
      box.paint(g);
    }

    if (signature != null) {
      // 1000 + 700 > 12
      int fontsize = (int)Math.round((pwidth + pheight) * 12. / 1700.);
      g.setColor(Color.blue);
      //g.setFont(new Font("SansSerif",Font.PLAIN,fontsize));
      g.setFont(new Font("SansSerif",Font.BOLD,fontsize));
      int w = ((FontMetrics)g.getFontMetrics()).stringWidth(signature); 
      g.drawString(signature,pwidth-w-fontsize,pheight-fontsize/2);

    }

  }
  // --------------------------------------------------------------------------
  public void saveImage() {
    int width = pwidth;
    int height = pheight;
    String imageFile = "mcloud.png";
    try {
      File output = new File(imageFile);
    
      BufferedImage image = new BufferedImage(width,height,BufferedImage.TYPE_USHORT_555_RGB);
      Graphics2D g = image.createGraphics();

      // white bg
      g.setColor(Color.white);
      g.fillRect(0,0,width,height);
      draw(g);

      try { 
        javax.imageio.ImageIO.setUseCache(false);
        javax.imageio.ImageIO.write(image, "png", output);
        System.err.println("File " + imageFile  + " written");
      }
      catch (NoClassDefFoundError e) {
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
