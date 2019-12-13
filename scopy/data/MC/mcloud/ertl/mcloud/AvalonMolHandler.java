
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

import java.util.Random;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;

import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;

import java.io.BufferedReader;
import java.io.StringReader;

import avalon.jni.JNISmi2Mol;

import novartis.chemistry.molecule.Molecule;
import novartis.chemistry.molecule.StandardMOL;
import novartis.chemistry.molecule.MoleculeDepicter;

import novartis.utilities.Box;
import novartis.utilities.FortranInputStream;
import novartis.utilities.MetaDraw;

// MolHandler rendering the chemical structures through the Avalon Toolkit

class AvalonMolHandler implements MolHandler {
  final static double STD_BOND = 1.54;
  final static double STD_BOX  = 200.0;
  private double width, height;

  StandardMOL m = null;
  String[] commands = null;
  Box box = null;
  double scale = 1.0;

  public AvalonMolHandler(String smiles)
  {
    // random width and height for testing
    double a = 150. + Math.random() * 300.;
    double b = 150. + Math.random() * 300.;
    width = Math.max(a,b);
    height = Math.min(a,b);

    try
    {
      String molfile = JNISmi2Mol.getSmi2Mol().smiToMOL(smiles);
      BufferedReader br = new BufferedReader(new StringReader(molfile.toString()));
      FortranInputStream in = new FortranInputStream(br);
      m = new StandardMOL();
      m.readMOLFile(in);
      br.close();
      m.recolor();
      box = new Box(0, 0, 0, 0);
      double len;
      len = MoleculeDepicter.getWorldWindow(m, box);
      MoleculeDepicter d = new MoleculeDepicter();
      commands = d.computeDepiction(m, 0, 0, MoleculeDepicter.USE_COLORS, null, null);
      width = STD_BOX*len/STD_BOND*box.getWidth();
      height = STD_BOX*len/STD_BOND*box.getHeight();
    }
    catch (Exception e)
    {
    }
  }

  public void rescale(double scale)
  {
    // rectangle area should be scaled by "scale" 
    // therefore rectangle sides need to be scaled by sqrt(scale)
    this.scale *= Math.sqrt(scale);
  }

  public double getWidth()
  {
    return this.scale*width;
  }

  public double getHeight()
  {
    return this.scale*height;
  }

  public void draw(Graphics2D g, double x, double y, double w, double h)
  {
    // shade
    BasicStroke orgStroke = (BasicStroke)g.getStroke(); 
    BasicStroke borderStroke = new BasicStroke(2.f,BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND);
    g.setColor(Color.gray);
    g.setStroke(borderStroke);
    g.draw(new Line2D.Double(x-w+1,y-h+scale*height,x-w+scale*width,y-h+scale*height));
    g.draw(new Line2D.Double(x-w+scale*width,y-h+1,x-w+scale*width,y-h+scale*height));
    g.setStroke(orgStroke);
    // clip
    Rectangle orgClip = g.getClipBounds();
    g.setClip(new Rectangle2D.Double(x-w,y-h,scale*width,scale*height));

    // coloring the rectangle's background in random pastell color
    Random random = new Random();
    int red = (int)(0.9*255+0.1*random.nextInt(255));
    int green = (int)(0.9*255+0.1*random.nextInt(255));
    int blue = (int)(0.9*255+0.1*random.nextInt(255));
    Color color = new Color(red,green,blue);
    g.setColor(color);
    g.fill(new Rectangle2D.Double(x-w,y-h,scale*width,scale*height));
    g.setColor(Color.black);
    g.draw(new Rectangle2D.Double(x-w,y-h,scale*width,scale*height));

    g.setClip(orgClip);
    // really draw the molecule
    MetaDraw.depictCommands(g, commands, new Rectangle((int)(x-w),(int)(y-h),(int)(scale*width),(int)(scale*height)), true, true);
  }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
