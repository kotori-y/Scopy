
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

interface MolHandler {

  // construct molecule and depiction from SMILES, 
  // MolHandler(String smiles);

  // rescale the molecule
  void rescale(double rescale);

  // width of the molecule box
  double getWidth();

  // height of the molecule box
  double getHeight();

  // draws molecule into center x, y and within box w, h
  void draw(Graphics2D g, double x, double y, double w, double h);

}

