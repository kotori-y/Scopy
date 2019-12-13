
/*
 *
 * Molecule Cloud, Copyright (c) 2012, Novartis Institutes for BioMedical Research Inc
 * written by Peter Ertl and Bernhard Rohde
 *
 * this software is released under the terms of FreeBSD license,
 * see license.txt in the distribution
 *
 */

package ertl.mcloud;

import java.awt.*;
import java.util.Random;
import java.awt.geom.*;

// "fake" MolHandler, draws just rectangles of random size instead of molecules
// to test the layout engine
// you should implement something like this using real molecule processing
// using Avalon Cheminformatics Toolkit, or other cheminformatics engine

class TestMolHandler implements MolHandler {
  //private Molecule molecule;
  private double width, height;

  public TestMolHandler(String smiles) {
    //molecule = new Molecule(smiles);

    // random width and height for testing
    double a = 150. + Math.random() * 300.;
    double b = 150. + Math.random() * 300.;
    a = 300.; b = 200.;
    width = Math.max(a,b);
    height = Math.min(a,b);

  }

  public void rescale(double scale) {
    // rectanle area should be scaled by "scale" 
    // therefore rectangle sides need to be scaled by sqrt(scale)
    double sq = Math.sqrt(scale);

    // molecule.rescale(scale);
    //width = molecule.getWidth();
    //height = molecule.getHeight();

    // for testing
    width = width * sq;
    height = height * sq;
  }

  public double getWidth() {
    return width;
  }

  public double getHeight() {
    return height;
  }

  public void draw(Graphics2D g, double x, double y, double w, double h) {
    // shade
    BasicStroke orgStroke = (BasicStroke)g.getStroke(); 
    BasicStroke borderStroke = new BasicStroke(2.f,BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND);
    g.setColor(Color.gray);
    g.setStroke(borderStroke);
    g.draw(new Line2D.Double(x-w+1,y-h+height,x-w+width,y-h+height));
    g.draw(new Line2D.Double(x-w+width,y-h+1,x-w+width,y-h+height));
    g.setStroke(orgStroke);
    // clip
    Rectangle orgClip = g.getClipBounds();
    g.setClip(new Rectangle2D.Double(x-w,y-h,width,height));

    // draw molecule 
    // add real molecule drawing code here
    // molecule.draw(g);

    // drawing just rectangle for testing 
    Random random = new Random();
    int red = random.nextInt(255);
    int green = random.nextInt(255);
    int blue = random.nextInt(255);
    Color color = new Color(red,green,blue);
    g.setColor(color);
    g.fill(new Rectangle2D.Double(x-w,y-h,width,height));
    g.setColor(Color.black);
    g.draw(new Rectangle2D.Double(x-w,y-h,width,height));


    g.setClip(orgClip);
  }

  static {
    System.err.println("\n");
    System.err.println("TestMolHandler does not display molecules, only color rectangles");
    System.err.println("To display molecules you have to connect AvalonMolHandler.");
    System.err.println("See readme.txt for instructions.");
    System.err.println("\n");
  }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
