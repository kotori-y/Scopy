package ertl.mcloud;

import javax.swing.*;
import java.awt.*;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class MolPanel extends JPanel {
    int xsize, ysize;
    MCloud mcloud; // defined draw method (used also in image saving)

    MolPanel(double xsize, double ysize, MCloud mcloud) {
        super();

        this.xsize = (int)xsize;
        this.ysize = (int)ysize;
        this.mcloud = mcloud;

        setBackground(Color.white);
        setPreferredSize(new Dimension((int)xsize,(int)ysize));
    }
    // --------------------------------------------------------------------------
    public void paintComponent(Graphics gg) {
        super.paintComponent(gg); // white bg
        Graphics2D g = (Graphics2D)gg;
        mcloud.draw(g);
    }
}