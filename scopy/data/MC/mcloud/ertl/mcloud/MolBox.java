package ertl.mcloud;

import java.awt.*;
import java.awt.geom.Rectangle2D;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class MolBox implements Comparable {
    private String smiles = null;
    private String data = null; // coloring etc
    private MolHandler molHandler = null;

    double x,y; // center
    double width,height;
    double w,h; // width/2 ....
    double scale;

    static boolean showMolecules = false;
    static double avrgWidth,avrgHeight;
    static int sortType = 0;
    static final int SIZE = 1, COUNT = 2;

    MolBox(String smiles, double scale, String data) {
        this.smiles = smiles;
        this.scale = scale;
        this.data = data;
        // test handler just to display rectangles (for testing)
        //molHandler = new TestMolHandler(smiles);

        // see readme txt for how to connect the Avalon Handler
        molHandler = new AvalonMolHandler(smiles);
    }

    void rescale(double rescale) {
        molHandler.rescale(rescale);
        width = molHandler.getWidth();
        height = molHandler.getHeight();
        w = width / 2.;
        h = height / 2.;
    }

    void setCenter(double x, double y) {
        this.x = (int)Math.round(x);
        this.y = (int)Math.round(y);
    }

    public String toString() {
        String s = "MolBox:";
        s += " " + smiles;
        return s;
    }

    void paint(Graphics2D g) {
        if (showMolecules) {
            molHandler.draw(g,x,y,w,h);
        }
        else {
            g.setColor(Color.red);
            g.draw(new Rectangle2D.Double(x-w,y-h,width,height));
            g.setColor(Color.black);
        }
    }

    // smaller on the top (so that they are visible)
    public int compareTo(Object o) {
        MolBox box = (MolBox)o;

        if (sortType == COUNT) {
            if (scale > box.scale) return 1;
            else if (scale < box.scale) return -1;
        }
        else if (sortType == SIZE) {
            if (width * height < box.width * box.height) return 1;
            else if (width * height > box.width * box.height) return -1;
        }
        else
            throw new RuntimeException("Bad sort type");
        return 0;
    }

    public void print() {
        System.out.println(this);
    }

}