package gov.nih.ncgc.viz;

import java.awt.*;
import java.awt.geom.*;
import java.awt.font.TextLayout;
import java.awt.image.BufferedImage;
import javax.swing.UIManager;


public class LabelRenderer {
    private Font font = UIManager.getFont("Label.font")
        .deriveFont(11.f).deriveFont(Font.BOLD);
        
    private Color textBackgroundColor = new Color (0x10, 0x10, 0x10, 140);
    private RoundRectangle2D rect = new RoundRectangle2D.Double();
    // fill label to the width of the component
    private boolean fillLabel = true; 
    private boolean labelTop = false; // show label at bottom
    private int labelOffset = 0;
        
    public LabelRenderer () {
    }
        
    public void setLabelFilled (boolean fillLabel) { 
        this.fillLabel = fillLabel; 
    }
    public boolean getLabelFilled () { return fillLabel; }
    public void setLabelTop (boolean labelTop) { 
        this.labelTop = labelTop; 
    }
    public boolean getLabelTop () { return labelTop; }
    public void setLabelFont (Font font) { this.font = font; }
    public Font getLabelFont () { return font; }
    public void setLabelOffset (int offset) { this.labelOffset = offset; }
    public int getLabelOffset () { return labelOffset; }
        
    public void renderLabel (Graphics2D g2, int offset, 
                             int padx, int pady, 
                             int width, int height, 
                             String label) {
        if (label == null) {
            return;
        }
            
        BufferedImage img = g2.getDeviceConfiguration()
            .createCompatibleImage(width, height, Transparency.TRANSLUCENT);
        Graphics2D g2img = img.createGraphics();
            
        g2.setFont(font);
        g2.setRenderingHint(RenderingHints.KEY_RENDERING, 
                            RenderingHints.VALUE_RENDER_QUALITY);
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
                            RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, 
                            RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            
        TextLayout tl = new TextLayout 
            (label, font, g2.getFontRenderContext());
        Rectangle2D box = tl.getBounds();
        double bw = fillLabel ? (width-2*padx) : box.getWidth();
        double bh = box.getHeight();
            
        int arc = getRoundedDiameter (bh);
        int yoff = labelOffset == 0 ? offset : labelOffset;
        rect.setRoundRect
            ((width-(bw+2*padx))/2.,
             labelTop ? yoff : (height-(bh+2*pady+yoff)),
             bw+2*padx, bh+2*pady, arc, arc);
        g2.setPaint(textBackgroundColor);
        g2.fill(rect);
        g2.setPaint(Color.white);
        g2.drawString(label, (int)((width -box.getWidth()- padx)/2.+ 0.5),
                      (int)(rect.getY() + bh + pady + 0.5));
            
        g2img.dispose();
        g2.drawImage(img, 0, 0, img.getWidth(), img.getHeight(), null);
    }
        
    static int getRoundedDiameter (double controlHeight) {
        int roundedDiameter = (int) (controlHeight * .98 + 0.5);
        return roundedDiameter - roundedDiameter % 2;
    }
} // LabelRenderer
